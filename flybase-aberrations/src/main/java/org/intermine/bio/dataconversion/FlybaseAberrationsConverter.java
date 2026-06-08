package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2026 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.BufferedReader;
import java.io.Reader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loader for FlyBase aberrations (FBab) and balancers (FBba).
 * Inputs in src.data.dir:
 *   - fb_synonym_*.tsv                                  identity for FBab + FBba
 *   - aberration_experimental_gene_del_dup_data_*.tsv   gene del/dup wiring
 *   - fbba_to_fbab*.tsv                                 curated balancer composition
 *
 * Design notes (revised 2026-06-02 after AllianceMineDev integrate found empty
 * collection tables):
 *   - File processing order is non-deterministic. All three parsers must be
 *     order-tolerant: a delDup or curated-balancer row that arrives before the
 *     corresponding synonyms row must still create the stub Aberration/Balancer
 *     item. processSynonyms enriches existing items rather than overwriting.
 *   - Any .gz companion file is skipped — if both fb_*.tsv and fb_*.tsv.gz are
 *     in src.data.dir (operator left the gz around after gunzip), the
 *     decompressed .tsv is the source of truth.
 *
 * @author new_flymine session, 2026-06-01 (bug fix 2026-06-02)
 */
public class FlybaseAberrationsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseAberrationsConverter.class);

    private static final String DATASET_TITLE = "FlyBase aberrations and balancers";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    private static final Pattern ABERRATION_PREFIX =
        Pattern.compile("^(Df|Dp|In|T)\\(.*\\).*");

    private final Map<String, Item> aberrationsById = new HashMap<String, Item>();
    private final Map<String, Item> aberrationsBySymbol = new HashMap<String, Item>();
    private final Map<String, Item> balancersById   = new HashMap<String, Item>();
    private final Map<String, Item> genesByFbgn     = new HashMap<String, Item>();
    /** FBab -> set of distinct cytological band strings (e.g. "51A5", "26A2-26A5").
     *  The breakpoints TSV is per-(aberration, breakpoint feature) but the cyto_loc
     *  column repeats per-aberration content across all of an FBab's breakpoint rows;
     *  dedupe here to avoid producing multiple identical CytologicalBand items. */
    private final Map<String, Set<String>> bandsByFbab = new HashMap<String, Set<String>>();
    private Item organism;

    /**
     * Constructor
     * @param writer the ItemWriter
     * @param model the target data model
     */
    public FlybaseAberrationsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    /**
     * Dispatch by input filename. Skips any *.gz companions so a leftover
     * compressed file does not double-process content as garbage.
     * @param reader reader over the current input file
     */
    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-aberrations: skipping gz companion " + name);
            return;
        }
        LOG.info("flybase-aberrations: processing " + name
                 + " (state before: aberrations=" + aberrationsById.size()
                 + " balancers=" + balancersById.size()
                 + " genes=" + genesByFbgn.size() + ")");
        if (name.startsWith("fb_synonym")) {
            processSynonyms(reader);
        } else if (name.startsWith("aberration_experimental_gene_del_dup_data")) {
            processDelDup(reader);
        } else if (name.startsWith("fbba_to_fbab")) {
            processCuratedBalancers(reader);
        } else if (name.startsWith("aberration_cytological_breakpoints")) {
            processBreakpoints(reader);
        } else {
            LOG.info("flybase-aberrations: unrecognised filename " + name + " — skipping");
        }
        LOG.info("flybase-aberrations: after " + name
                 + ": aberrations=" + aberrationsById.size()
                 + " balancers=" + balancersById.size()
                 + " genes=" + genesByFbgn.size());
    }

    /**
     * Get the Aberration item for the given FBab; create a stub (primaryIdentifier
     * + organism only) if not yet seen. The stub is enriched later by
     * processSynonyms with symbol + aberrationType.
     */
    private Item getOrCreateAberration(String fbab) throws ObjectStoreException {
        Item ab = aberrationsById.get(fbab);
        if (ab == null) {
            ab = createItem("Aberration");
            ab.setAttribute("primaryIdentifier", fbab);
            ab.setReference("organism", getOrganism());
            aberrationsById.put(fbab, ab);
        }
        return ab;
    }

    /**
     * Get the Balancer item for the given FBba; create a stub if not yet seen.
     */
    private Item getOrCreateBalancer(String fbba) throws ObjectStoreException {
        Item ba = balancersById.get(fbba);
        if (ba == null) {
            ba = createItem("Balancer");
            ba.setAttribute("primaryIdentifier", fbba);
            ba.setReference("organism", getOrganism());
            balancersById.put(fbba, ba);
        }
        return ba;
    }

    /**
     * Get the Gene item for the given FBgn; create a stub if not yet seen.
     * The InterMine loader merges by primaryIdentifier (per
     * flybase-aberrations_keys.properties) so a stub created here will be
     * merged with the chado-db gene of the same FBgn at load time.
     */
    private Item getOrCreateGene(String fbgn) throws ObjectStoreException {
        Item g = genesByFbgn.get(fbgn);
        if (g == null) {
            g = createItem("Gene");
            g.setAttribute("primaryIdentifier", fbgn);
            g.setReference("organism", getOrganism());
            genesByFbgn.put(fbgn, g);
        }
        return g;
    }

    /**
     * Parse fb_synonym_*.tsv. Enriches (or creates) Aberration / Balancer items
     * with symbol + (for Aberration) aberrationType. Order-tolerant: if an
     * item already exists from a prior del/dup or curated-balancer pass, this
     * sets the missing attributes without losing already-attached collections.
     */
    void processSynonyms(Reader reader) throws Exception {
        BufferedReader br = new BufferedReader(reader);
        String line;
        int abCount = 0;
        int baCount = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < 3) {
                continue;
            }
            String fbId = cols[0];
            String symbol = cols[2];
            if (fbId.startsWith("FBab")) {
                Item ab = getOrCreateAberration(fbId);
                if (!symbol.isEmpty()) {
                    if (ab.getAttribute("symbol") == null) {
                        ab.setAttribute("symbol", symbol);
                    }
                    if (ab.getAttribute("aberrationType") == null) {
                        ab.setAttribute("aberrationType", aberrationTypeFromSymbol(symbol));
                    }
                    // Index by symbol for the curated-balancers symbol-ref lookup.
                    // Only first-seen wins on duplicate symbols (synonyms file has
                    // one row per FBab so this is normally a no-op).
                    aberrationsBySymbol.putIfAbsent(symbol, ab);
                } else if (ab.getAttribute("aberrationType") == null) {
                    ab.setAttribute("aberrationType", "other");
                }
                abCount++;
            } else if (fbId.startsWith("FBba")) {
                Item ba = getOrCreateBalancer(fbId);
                if (!symbol.isEmpty() && ba.getAttribute("symbol") == null) {
                    ba.setAttribute("symbol", symbol);
                }
                baCount++;
            }
        }
        LOG.info("flybase-aberrations: processSynonyms rows: FBab=" + abCount
                 + " FBba=" + baCount);
    }

    /**
     * Parse aberration_experimental_gene_del_dup_data_*.tsv. Skips "not deleted"
     * / "not duplicated" negative rows. Order-tolerant — creates stub Aberration
     * + Gene items if not yet seen; processSynonyms / chado-db will merge or
     * enrich them.
     */
    void processDelDup(Reader reader) throws Exception {
        BufferedReader br = new BufferedReader(reader);
        String line;
        int delCount = 0;
        int dupCount = 0;
        int skipped = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < 5) {
                continue;
            }
            String fbgn = cols[0];
            String type = cols[2].toLowerCase();
            String fbab = cols[3];
            if (type.startsWith("not ")) {
                skipped++;
                continue;
            }
            Item aberration = getOrCreateAberration(fbab);
            Item gene = getOrCreateGene(fbgn);
            if (type.contains("deleted") || type.contains("disrupted")) {
                aberration.addToCollection("deletedGenes", gene);
                delCount++;
            } else if (type.contains("duplicated")) {
                aberration.addToCollection("duplicatedGenes", gene);
                dupCount++;
            }
        }
        LOG.info("flybase-aberrations: processDelDup deletedRows=" + delCount
                 + " duplicatedRows=" + dupCount + " negativesSkipped=" + skipped);
    }

    /**
     * Parse the curated fbba_to_fbab*.tsv. Each row lists a balancer (FBba)
     * and the FBab inversions that compose it. The fbab_ids column accepts
     * EITHER FBab primary identifiers (FBab0000001, etc.) OR symbols of
     * constituent inversions (e.g. In(2LR)CyO). Symbol refs are resolved
     * against the aberrationsBySymbol map built from the fb_synonym pass —
     * curators don't need to look up FBab IDs by hand for the well-known
     * balancers whose inversion symbols are stable across FlyBase releases.
     *
     * Order-tolerant: creates stub Balancer + Aberration items as needed.
     * Symbol refs that don't resolve are skipped silently (logged at INFO).
     */
    void processCuratedBalancers(Reader reader) throws Exception {
        BufferedReader br = new BufferedReader(reader);
        String line;
        int balancerRows = 0;
        int composedAdds = 0;
        int symbolResolved = 0;
        int symbolUnresolved = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < 3) {
                continue;
            }
            String fbba = cols[0];
            String fbabRefs = cols[2];
            Item balancer = getOrCreateBalancer(fbba);
            balancerRows++;
            if (fbabRefs.isEmpty()) {
                continue;
            }
            for (String ref : fbabRefs.split("\\|")) {
                String trimmed = ref.trim();
                if (trimmed.isEmpty()) {
                    continue;
                }
                Item aberration;
                if (trimmed.startsWith("FBab")) {
                    // Direct FBab id: stub-create if not yet seen.
                    aberration = getOrCreateAberration(trimmed);
                } else {
                    // Symbol reference: must already be in the synonyms-built map.
                    aberration = aberrationsBySymbol.get(trimmed);
                    if (aberration == null) {
                        symbolUnresolved++;
                        LOG.info("flybase-aberrations: curated balancer "
                                 + fbba + " references unknown aberration symbol '"
                                 + trimmed + "' — skipped");
                        continue;
                    }
                    symbolResolved++;
                }
                balancer.addToCollection("composedOfAberrations", aberration);
                composedAdds++;
            }
        }
        LOG.info("flybase-aberrations: processCuratedBalancers rows=" + balancerRows
                 + " composedAdds=" + composedAdds
                 + " symbolResolved=" + symbolResolved
                 + " symbolUnresolved=" + symbolUnresolved);
    }

    /**
     * Parse the builder-side aberration_cytological_breakpoints.tsv dump
     * (the chado SQL extraction documented in
     * agr_intermine_builder/docs/FLYMINE_SIBLING_REPLY_2026_06_08.md §2).
     *
     * Columns (tab-separated, empty for NULL):
     *   1: fbab             (FBab uniquename)
     *   2: breakpoint_name  (chromosome_breakpoint feature, e.g. Df(2R)03072:bk1)
     *   3: genomic_loc      (text or empty, e.g. X_r6:7901330..7956278)
     *   4: chrom            (parsed from genomic_loc)
     *   5: fmin             (int)
     *   6: fmax             (int)
     *   7: cyto_loc         (per-FBab derived_attributed_breakpoint;
     *                        semicolon-separated for multi-breakpoint, e.g.
     *                        "51A5;51C1" or "25E1-25E2;26A2-26A5")
     *
     * The cyto_loc field is per-aberration (NOT per-breakpoint). Every row
     * for a given FBab carries the same cyto_loc text. To avoid producing
     * duplicate CytologicalBand items we dedupe (fbab, band) pairs via
     * bandsByFbab; the CytologicalBand items themselves are created in
     * close() after all rows are read.
     *
     * v1 scope: only cytological bands (col 7). Per-breakpoint genomic
     * coords (cols 3-6) are deferred — would require a separate
     * Chromosome+Location item chain per breakpoint feature, and the
     * coverage (9,918 of 68,591 rows = 14%) is sparse enough that a
     * dedicated v2 follow-up is cleaner than a partial v1.
     */
    void processBreakpoints(Reader reader) throws Exception {
        BufferedReader br = new BufferedReader(reader);
        String line;
        int rows = 0;
        int distinctFbabs = 0;
        int newBandPairs = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < 7) {
                continue;
            }
            String fbab = cols[0].trim();
            if (fbab.isEmpty() || !fbab.startsWith("FBab")) {
                continue;
            }
            rows++;
            String cytoLoc = cols[6].trim();
            if (cytoLoc.isEmpty()) {
                continue;
            }
            Set<String> bands = bandsByFbab.get(fbab);
            if (bands == null) {
                bands = new LinkedHashSet<String>();
                bandsByFbab.put(fbab, bands);
                distinctFbabs++;
            }
            for (String band : cytoLoc.split(";")) {
                String trimmed = band.trim();
                if (trimmed.isEmpty()) {
                    continue;
                }
                if (bands.add(trimmed)) {
                    newBandPairs++;
                }
            }
        }
        LOG.info("flybase-aberrations: processBreakpoints rows=" + rows
                 + " distinct FBabs with cyto_loc=" + distinctFbabs
                 + " distinct (fbab, band) pairs=" + newBandPairs);
    }

    /**
     * Map FlyBase aberration symbol prefix to a coarse aberrationType.
     * @param symbol the current_symbol field from fb_synonym (e.g. "Df(2R)min")
     * @return one of: deletion | duplication | inversion | translocation | other
     */
    static String aberrationTypeFromSymbol(String symbol) {
        if (symbol == null) {
            return "other";
        }
        Matcher m = ABERRATION_PREFIX.matcher(symbol);
        if (!m.matches()) {
            return "other";
        }
        switch (m.group(1)) {
            case "Df": return "deletion";
            case "Dp": return "duplication";
            case "In": return "inversion";
            case "T":  return "translocation";
            default:   return "other";
        }
    }

    Item getOrganism() throws ObjectStoreException {
        if (organism == null) {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", DMEL_TAXON);
            store(organism);
        }
        return organism;
    }

    /**
     * Deferred store: items are accumulated during process() across multiple
     * input files in non-deterministic order, mutated as more files are read,
     * then stored at the end so collection mutations from later files persist.
     * Store order: Aberration -> Balancer -> Gene.
     */
    @Override
    public void close() throws Exception {
        // Materialise CytologicalBand items from the deduped per-FBab
        // band set BEFORE storing aberrations, so the
        // cytologicalBreakpoints collection refs land on the persisted
        // Item. CytologicalBand extends Location; v1 only sets
        // cytologicalCoordinates (inherited start/end/strand/locatedOn
        // left null since chado's band-level granularity has no
        // molecular coordinates to map cleanly).
        int totalBands = 0;
        for (Map.Entry<String, Set<String>> e : bandsByFbab.entrySet()) {
            Item aberration = getOrCreateAberration(e.getKey());
            for (String band : e.getValue()) {
                Item cb = createItem("CytologicalBand");
                cb.setAttribute("cytologicalCoordinates", band);
                store(cb);
                aberration.addToCollection("cytologicalBreakpoints", cb);
                totalBands++;
            }
        }
        LOG.info("flybase-aberrations: close(): storing "
                 + aberrationsById.size() + " aberrations, "
                 + balancersById.size() + " balancers, "
                 + genesByFbgn.size() + " genes, "
                 + totalBands + " CytologicalBand items (across "
                 + bandsByFbab.size() + " FBabs)");
        for (Item ab : aberrationsById.values()) {
            store(ab);
        }
        for (Item ba : balancersById.values()) {
            store(ba);
        }
        for (Item g : genesByFbgn.values()) {
            store(g);
        }
        super.close();
    }
}
