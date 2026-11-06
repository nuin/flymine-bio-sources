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
 * Loader for FlyBase precomputed gene summaries.
 *
 * Reads two files from src.data.dir (dispatched by filename):
 *   best_gene_summary_*.tsv       -> Gene.description
 *   automated_gene_summaries_*.tsv -> Gene.briefDescription
 *
 * best_gene_summary header:
 *   #FBgn_ID  Gene_Symbol  Summary_Source  Summary
 *   col 0 = FBgn (Gene.primaryIdentifier), col 3 = Summary (Gene.description)
 *
 * automated_gene_summaries header:
 *   #primary_FBgn  summary_text
 *   col 0 = FBgn, col 1 = summary_text (Gene.briefDescription)
 *
 * Gene rows already exist in the merged DB from chado-db-flybase-dmel; this
 * source enriches them by FBgn merge (Gene.key = primaryIdentifier). Blank
 * summaries are skipped. Genes are cached across both files and stored once
 * in close().
 *
 * @author flymine builder session, 2026-09
 */
public class FlybaseGeneDescriptionsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseGeneDescriptionsConverter.class);

    private static final String DATASET_TITLE = "FlyBase gene summaries";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    // best_gene_summary columns
    private static final int BEST_COL_FBGN = 0;
    private static final int BEST_COL_SUMMARY = 3;
    private static final int BEST_MIN_COLS = 4;

    // automated_gene_summaries columns
    private static final int AUTO_COL_FBGN = 0;
    private static final int AUTO_COL_SUMMARY = 1;
    private static final int AUTO_MIN_COLS = 2;

    private final Map<String, Item> genesByFbgn = new HashMap<String, Item>();
    private Item organism;

    /**
     * Constructor.
     * @param writer the ItemWriter
     * @param model the target data model
     */
    public FlybaseGeneDescriptionsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-gene-descriptions: skipping gz companion " + name);
            return;
        }
        if (name.startsWith("best_gene_summary")) {
            processBest(reader);
        } else if (name.startsWith("automated_gene_summaries")) {
            processAutomated(reader);
        } else {
            LOG.info("flybase-gene-descriptions: unrecognised file (" + name
                     + "), skipping");
        }
    }

    private void processBest(Reader reader) throws Exception {
        LOG.info("flybase-gene-descriptions: processing best_gene_summary");
        BufferedReader br = new BufferedReader(reader);
        String line;
        int rows = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < BEST_MIN_COLS) {
                continue;
            }
            String fbgn = cols[BEST_COL_FBGN].trim();
            String summary = cols[BEST_COL_SUMMARY].trim();
            if (fbgn.isEmpty() || summary.isEmpty()) {
                continue;
            }
            getOrCreateGene(fbgn).setAttribute("description", summary);
            rows++;
        }
        LOG.info("flybase-gene-descriptions: best_gene_summary rows=" + rows);
    }

    private void processAutomated(Reader reader) throws Exception {
        LOG.info("flybase-gene-descriptions: processing automated_gene_summaries");
        BufferedReader br = new BufferedReader(reader);
        String line;
        int rows = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < AUTO_MIN_COLS) {
                continue;
            }
            String fbgn = cols[AUTO_COL_FBGN].trim();
            String summary = cols[AUTO_COL_SUMMARY].trim();
            if (fbgn.isEmpty() || summary.isEmpty()) {
                continue;
            }
            getOrCreateGene(fbgn).setAttribute("briefDescription", summary);
            rows++;
        }
        LOG.info("flybase-gene-descriptions: automated_gene_summaries rows=" + rows);
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
        LOG.info("flybase-gene-descriptions: close(): storing "
                 + genesByFbgn.size() + " genes");
        for (Item g : genesByFbgn.values()) {
            store(g);
        }
        super.close();
    }
}
