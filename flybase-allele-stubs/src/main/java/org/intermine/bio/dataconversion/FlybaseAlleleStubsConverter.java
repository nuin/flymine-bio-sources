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

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loader for the canonical FlyBase FBal -> FBgn mapping
 * (precomputed file fbal_to_fbgn_*.tsv).
 *
 * The FB2026_01 precomputed file contains 305,332 FBal entries; chado's
 * Dmel-filtered feature table only carries 239,366. The ~66K delta is the
 * transgene-vector alleles (Ecol\lacZ, Scer\GAL4/UAS, Avic\GFP, Ppyr\luc,
 * Hsap\* fusions) which in chado are tagged with their source organism
 * rather than Dmel; chado-db-flybase-dmel's organism filter drops them.
 * For a fly geneticist these alleles are essential (every UAS driver
 * line, every reporter, every fluorescent fusion is one of these).
 *
 * This source loads ALL FBal rows from the precomputed file. The
 * InterMine loader merges by Allele.primaryIdentifier:
 *   - FBals already in chado-loaded data: no-op merge (same attrs)
 *   - FBals only in the precomputed file: fresh Allele stubs created
 *     with primaryIdentifier + symbol + gene ref + organism=Dmel
 *
 * Filename:
 *   fbal_to_fbgn_*.tsv
 *
 * Columns (tab-separated):
 *   col 1: FBal{nnnnnnn}                -> Allele.primaryIdentifier
 *   col 2: allele symbol (e.g. Ecol\lacZ[T125])  -> Allele.symbol
 *   col 3: FBgn{nnnnnnn}                -> Gene reference (merged by FBgn)
 *   col 4: gene symbol (ignored — chado/FBgn carries it)
 *
 * Header lines start with '#'. Empty / malformed rows are skipped.
 * Defence-in-depth: dispatches by filename so a misconfigured
 * src.data.dir.includes does not misparse a sibling file.
 *
 * @author new_flymine session, 2026-06-04
 */
public class FlybaseAlleleStubsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseAlleleStubsConverter.class);

    private static final String DATASET_TITLE = "FlyBase allele-gene mapping (canonical FBal set)";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    private static final int COL_FBAL = 0;
    private static final int COL_SYMBOL = 1;
    private static final int COL_FBGN = 2;
    private static final int MIN_COLS = 4;

    private final Map<String, Item> allelesByFbal = new HashMap<String, Item>();
    private final Map<String, Item> genesByFbgn = new HashMap<String, Item>();
    private Item organism;

    /**
     * Constructor.
     * @param writer the ItemWriter
     * @param model the target data model
     */
    public FlybaseAlleleStubsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-allele-stubs: skipping gz companion " + name);
            return;
        }
        if (!name.startsWith("fbal_to_fbgn")) {
            LOG.info("flybase-allele-stubs: not an fbal_to_fbgn file ("
                     + name + "), skipping");
            return;
        }
        LOG.info("flybase-allele-stubs: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int alleleRows = 0;
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
            String fbal = cols[COL_FBAL].trim();
            String symbol = cols[COL_SYMBOL].trim();
            String fbgn = cols[COL_FBGN].trim();
            if (fbal.isEmpty() || !fbal.startsWith("FBal")) {
                skipped++;
                continue;
            }
            Item allele = getOrCreateAllele(fbal);
            if (!symbol.isEmpty()) {
                allele.setAttribute("symbol", symbol);
            }
            if (!fbgn.isEmpty() && fbgn.startsWith("FBgn")) {
                Item gene = getOrCreateGene(fbgn);
                allele.setReference("gene", gene);
            }
            alleleRows++;
        }
        LOG.info("flybase-allele-stubs: alleleRows=" + alleleRows
                 + " skipped=" + skipped);
    }

    private Item getOrCreateAllele(String fbal) throws ObjectStoreException {
        Item a = allelesByFbal.get(fbal);
        if (a == null) {
            a = createItem("Allele");
            a.setAttribute("primaryIdentifier", fbal);
            a.setReference("organism", getOrganism());
            allelesByFbal.put(fbal, a);
        }
        return a;
    }

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
        LOG.info("flybase-allele-stubs: close(): storing "
                 + allelesByFbal.size() + " alleles, "
                 + genesByFbgn.size() + " gene stubs");
        for (Item a : allelesByFbal.values()) {
            store(a);
        }
        for (Item g : genesByFbgn.values()) {
            store(g);
        }
        super.close();
    }
}
