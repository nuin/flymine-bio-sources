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
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loader for FlyBase transposable-element insertion sites (FBti).
 *
 * The FB2026_01 precomputed file ships 420,755 FBti records. Without this
 * source, Aim 2 cross-class queries ("stocks x transgenes x genes x
 * insertion sites") return empty and Aim 3 ID converter rejects FBti
 * inputs. The historical chado-db source doesn't surface FBti as a
 * first-class entity; this lean source closes the gap by parsing the
 * insertion_mapping_*.tsv file directly.
 *
 * Input filename: insertion_mapping_*.tsv (header lines start with `##`).
 *
 * Columns (tab-separated):
 *   col 1: FBti{nnnnnnn}
 *   col 2: symbol (e.g. P{lacW}gukh[k08210])
 *   col 3: genomic_location (e.g. 2L:7,238,500..7,238,500) — parsed into
 *          chromosome + start/end of a Location item
 *   col 4: range (e.g. 2L:7238500-7238500) — currently ignored (redundant
 *          with col 3 in most rows; can be added in v2)
 *   col 5: orientation (+ / - / blank) — Location.strand "1" / "-1"
 *   col 6: estimated cytogenetic location — currently ignored
 *   col 7: observed cytogenetic location — TransposableElementInsertionSite.cytogeneticLocation
 *
 * @author new_flymine session, 2026-06-04
 */
public class FlybaseInsertionsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseInsertionsConverter.class);

    private static final String DATASET_TITLE = "FlyBase transposable-element insertion sites";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    private static final int COL_FBTI = 0;
    private static final int COL_SYMBOL = 1;
    private static final int COL_GENOMIC_LOCATION = 2;
    private static final int COL_ORIENTATION = 4;
    private static final int COL_CYTO_OBSERVED = 6;
    private static final int MIN_COLS = 7;

    /** Matches "<chr>:<start>..<end>" with optional thousands separators. */
    private static final Pattern LOCATION_PATTERN =
        Pattern.compile("^([^:]+):([0-9,]+)\\.\\.([0-9,]+)$");

    private final Map<String, Item> insertionsByFbti = new HashMap<String, Item>();
    private final Map<String, Item> chromosomesByName = new HashMap<String, Item>();
    private Item organism;

    public FlybaseInsertionsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-insertions: skipping gz companion " + name);
            return;
        }
        if (!name.startsWith("insertion_mapping")) {
            LOG.info("flybase-insertions: not an insertion_mapping file ("
                     + name + "), skipping");
            return;
        }
        LOG.info("flybase-insertions: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int rows = 0;
        int withLocation = 0;
        int skipped = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < MIN_COLS) {
                skipped++;
                continue;
            }
            String fbti = cols[COL_FBTI].trim();
            if (fbti.isEmpty() || !fbti.startsWith("FBti")) {
                skipped++;
                continue;
            }
            String symbol = cols[COL_SYMBOL].trim();
            String genomic = cols[COL_GENOMIC_LOCATION].trim();
            String orientation = cols[COL_ORIENTATION].trim();
            String cytoObserved = cols[COL_CYTO_OBSERVED].trim();

            Item insertion = getOrCreateInsertion(fbti);
            if (!symbol.isEmpty()) {
                insertion.setAttribute("symbol", symbol);
            }
            if (!cytoObserved.isEmpty()) {
                insertion.setAttribute("cytogeneticLocation", cytoObserved);
            }
            if (!genomic.isEmpty()) {
                Item location = parseLocation(genomic, orientation, insertion);
                if (location != null) {
                    insertion.setReference("chromosomeLocation", location);
                    withLocation++;
                }
            }
            rows++;
        }
        LOG.info("flybase-insertions: rows=" + rows
                 + " withLocation=" + withLocation
                 + " skipped=" + skipped);
    }

    /**
     * Parse "2L:7,238,500..7,238,500" into a Location item; null if the
     * input doesn't match. Also creates and attaches a Chromosome item
     * with primaryIdentifier=2L.
     */
    private Item parseLocation(String genomic, String orientation, Item insertion)
            throws ObjectStoreException {
        Matcher m = LOCATION_PATTERN.matcher(genomic);
        if (!m.matches()) {
            return null;
        }
        String chrName = m.group(1);
        String startStr = m.group(2).replace(",", "");
        String endStr = m.group(3).replace(",", "");
        Item chr = getOrCreateChromosome(chrName);
        Item loc = createItem("Location");
        loc.setAttribute("start", startStr);
        loc.setAttribute("end", endStr);
        if ("+".equals(orientation)) {
            loc.setAttribute("strand", "1");
        } else if ("-".equals(orientation)) {
            loc.setAttribute("strand", "-1");
        }
        loc.setReference("locatedOn", chr);
        loc.setReference("feature", insertion);
        insertion.setReference("chromosome", chr);
        store(loc);
        return loc;
    }

    private Item getOrCreateInsertion(String fbti) throws ObjectStoreException {
        Item i = insertionsByFbti.get(fbti);
        if (i == null) {
            i = createItem("TransposableElementInsertionSite");
            i.setAttribute("primaryIdentifier", fbti);
            i.setReference("organism", getOrganism());
            insertionsByFbti.put(fbti, i);
        }
        return i;
    }

    private Item getOrCreateChromosome(String chrName) throws ObjectStoreException {
        Item c = chromosomesByName.get(chrName);
        if (c == null) {
            c = createItem("Chromosome");
            c.setAttribute("primaryIdentifier", chrName);
            c.setReference("organism", getOrganism());
            chromosomesByName.put(chrName, c);
        }
        return c;
    }

    private Item getOrganism() throws ObjectStoreException {
        if (organism == null) {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", DMEL_TAXON);
            store(organism);
        }
        return organism;
    }

    @Override
    public void close() throws Exception {
        LOG.info("flybase-insertions: close(): storing "
                 + insertionsByFbti.size() + " insertions, "
                 + chromosomesByName.size() + " chromosomes");
        for (Item i : insertionsByFbti.values()) {
            store(i);
        }
        for (Item c : chromosomesByName.values()) {
            store(c);
        }
        super.close();
    }
}
