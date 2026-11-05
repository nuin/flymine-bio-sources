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
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loader for FlyBase D. melanogaster - human orthologs.
 *
 * Reads one file from src.data.dir:
 *   dmel_human_orthologs_disease_*.tsv
 *
 * Header (leading "##" comment line):
 *   Dmel_gene_ID Dmel_gene_symbol Human_gene_HGNC_ID Human_gene_OMIM_ID
 *   Human_gene_symbol DIOPT_score OMIM_Phenotype_IDs OMIM_Phenotype_IDs[name]
 *
 *   col 0 = Dmel_gene_ID       (FBgn)  -> fly Gene.primaryIdentifier (merge key)
 *   col 2 = Human_gene_HGNC_ID (e.g. HGNC:8918, may be empty)
 *   col 4 = Human_gene_symbol  (e.g. PHEX)
 *
 * For each row: merges the Dmel Gene by FBgn (taxon 7227); creates/merges a
 * HUMAN Gene stub (taxon 9606) whose primaryIdentifier is the HGNC id when
 * present, else the human symbol, with symbol set; and creates a Homologue
 * (gene=fly, homologue=human, type="orthologue"). The Gene.homologues
 * collection is populated via the Homologue.gene reverse-reference.
 *
 * Creating cross-organism human Gene stubs is intended. Homologues are
 * deduped on (FBgn, humanIdentifier) since a fly gene can map to the same
 * human gene through multiple OMIM phenotype rows.
 *
 * @author flymine builder session, 2026-09
 */
public class FlybaseHumanOrthologsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseHumanOrthologsConverter.class);

    private static final String DATASET_TITLE =
        "FlyBase D. melanogaster-Human Orthologs";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";
    private static final String HUMAN_TAXON = "9606";

    private static final int COL_FBGN = 0;
    private static final int COL_HGNC = 2;
    private static final int COL_HUMAN_SYMBOL = 4;
    private static final int MIN_COLS = 5;

    private final Map<String, Item> flyGenesByFbgn = new HashMap<String, Item>();
    private final Map<String, Item> humanGenesById = new HashMap<String, Item>();
    private final Set<String> seenHomologues = new HashSet<String>();
    private Item dmelOrganism;
    private Item humanOrganism;

    /**
     * Constructor.
     * @param writer the ItemWriter
     * @param model the target data model
     */
    public FlybaseHumanOrthologsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-human-orthologs: skipping gz companion " + name);
            return;
        }
        LOG.info("flybase-human-orthologs: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int homologues = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < MIN_COLS) {
                continue;
            }
            String fbgn = cols[COL_FBGN].trim();
            String hgnc = cols[COL_HGNC].trim();
            String humanSymbol = cols[COL_HUMAN_SYMBOL].trim();
            String humanId = !hgnc.isEmpty() ? hgnc : humanSymbol;
            if (fbgn.isEmpty() || humanId.isEmpty()) {
                continue;
            }
            String dedupKey = fbgn + "\t" + humanId;
            if (!seenHomologues.add(dedupKey)) {
                continue;
            }
            Item flyGene = getOrCreateFlyGene(fbgn);
            Item humanGene = getOrCreateHumanGene(humanId, humanSymbol);

            Item homologue = createItem("Homologue");
            homologue.setAttribute("type", "orthologue");
            homologue.setReference("gene", flyGene);
            homologue.setReference("homologue", humanGene);
            store(homologue);
            homologues++;
        }
        LOG.info("flybase-human-orthologs: homologues=" + homologues);
    }

    private Item getOrCreateFlyGene(String fbgn) throws ObjectStoreException {
        Item g = flyGenesByFbgn.get(fbgn);
        if (g == null) {
            g = createItem("Gene");
            g.setAttribute("primaryIdentifier", fbgn);
            g.setReference("organism", getDmelOrganism());
            flyGenesByFbgn.put(fbgn, g);
        }
        return g;
    }

    private Item getOrCreateHumanGene(String humanId, String humanSymbol)
        throws ObjectStoreException {
        Item g = humanGenesById.get(humanId);
        if (g == null) {
            g = createItem("Gene");
            g.setAttribute("primaryIdentifier", humanId);
            if (!humanSymbol.isEmpty()) {
                g.setAttribute("symbol", humanSymbol);
            }
            g.setReference("organism", getHumanOrganism());
            humanGenesById.put(humanId, g);
        }
        return g;
    }

    private Item getDmelOrganism() throws ObjectStoreException {
        if (dmelOrganism == null) {
            dmelOrganism = createItem("Organism");
            dmelOrganism.setAttribute("taxonId", DMEL_TAXON);
            store(dmelOrganism);
        }
        return dmelOrganism;
    }

    private Item getHumanOrganism() throws ObjectStoreException {
        if (humanOrganism == null) {
            humanOrganism = createItem("Organism");
            humanOrganism.setAttribute("taxonId", HUMAN_TAXON);
            store(humanOrganism);
        }
        return humanOrganism;
    }

    @Override
    public void close() throws Exception {
        LOG.info("flybase-human-orthologs: close(): storing "
                 + flyGenesByFbgn.size() + " fly genes, "
                 + humanGenesById.size() + " human gene stubs");
        for (Item g : flyGenesByFbgn.values()) {
            store(g);
        }
        for (Item g : humanGenesById.values()) {
            store(g);
        }
        super.close();
    }
}
