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
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loader for FlyBase transgenic constructs (FBtp).
 *
 * 173,232 FBtp records in FB2026_01. Each FBtp is a molecular tool used
 * to make transgenic stocks (P{lacW}, P{UAS-GAL4}, P{Avic\GFP} etc.) and
 * is associated with: one or more genes it encodes (FBgn) and one or
 * more insertion events that use it (FBti).
 *
 * Without this source, Aim 2's "transgenes" entity is empty and Aim 4
 * BDSC stock pages can't link to construct descriptions.
 *
 * Input filename: transgenic_construct_descriptions_*.tsv (header lines
 * start with `##` or `#`).
 *
 * Columns (tab-separated):
 *   col 1: FBtp{nnnnnnn}
 *   col 2: construct symbol (e.g. P{lacW})
 *   col 3: encoded features (semicolon-space-separated text description)
 *   col 4: encoded FBgn list (semicolon-space-separated)
 *   col 5: related FBti symbols (semicolon-space-separated) — ignored,
 *          use col 6 ids
 *   col 6: related FBti IDs (semicolon-space-separated)
 *
 * @author new_flymine session, 2026-06-04
 */
public class FlybaseTransgenicConstructsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseTransgenicConstructsConverter.class);

    private static final String DATASET_TITLE = "FlyBase transgenic constructs";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    private static final int COL_FBTP = 0;
    private static final int COL_SYMBOL = 1;
    private static final int COL_ENCODED_FEATURES = 2;
    private static final int COL_ENCODED_FBGN = 3;
    private static final int COL_RELATED_FBTI = 5;
    private static final int MIN_COLS = 4;

    /** Semicolon possibly followed by whitespace. */
    private static final Pattern SEMI_SPLIT = Pattern.compile(";\\s*");

    private final Map<String, Item> constructsByFbtp = new HashMap<String, Item>();
    private final Map<String, Item> genesByFbgn = new HashMap<String, Item>();
    private final Map<String, Item> insertionsByFbti = new HashMap<String, Item>();
    private Item organism;

    public FlybaseTransgenicConstructsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-transgenic-constructs: skipping gz companion " + name);
            return;
        }
        if (!name.startsWith("transgenic_construct_descriptions")) {
            LOG.info("flybase-transgenic-constructs: not a "
                     + "transgenic_construct_descriptions file (" + name
                     + "), skipping");
            return;
        }
        LOG.info("flybase-transgenic-constructs: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int rows = 0;
        int geneRefs = 0;
        int insertRefs = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < MIN_COLS) {
                continue;
            }
            String fbtp = cols[COL_FBTP].trim();
            if (fbtp.isEmpty() || !fbtp.startsWith("FBtp")) {
                continue;
            }
            Item construct = getOrCreateConstruct(fbtp);
            if (cols.length > COL_SYMBOL) {
                String symbol = cols[COL_SYMBOL].trim();
                if (!symbol.isEmpty()) {
                    construct.setAttribute("symbol", symbol);
                }
            }
            if (cols.length > COL_ENCODED_FEATURES) {
                String features = cols[COL_ENCODED_FEATURES].trim();
                if (!features.isEmpty()) {
                    construct.setAttribute("description", features);
                }
            }
            if (cols.length > COL_ENCODED_FBGN) {
                for (String fbgn : SEMI_SPLIT.split(cols[COL_ENCODED_FBGN])) {
                    String trimmed = fbgn.trim();
                    if (trimmed.isEmpty() || !trimmed.startsWith("FBgn")) {
                        continue;
                    }
                    Item gene = getOrCreateGene(trimmed);
                    construct.addToCollection("encodedGenes", gene);
                    geneRefs++;
                }
            }
            if (cols.length > COL_RELATED_FBTI) {
                for (String fbti : SEMI_SPLIT.split(cols[COL_RELATED_FBTI])) {
                    String trimmed = fbti.trim();
                    if (trimmed.isEmpty() || !trimmed.startsWith("FBti")) {
                        continue;
                    }
                    Item insertion = getOrCreateInsertion(trimmed);
                    construct.addToCollection("relatedInsertions", insertion);
                    insertRefs++;
                }
            }
            rows++;
        }
        LOG.info("flybase-transgenic-constructs: rows=" + rows
                 + " geneRefs=" + geneRefs
                 + " insertRefs=" + insertRefs);
    }

    private Item getOrCreateConstruct(String fbtp) throws ObjectStoreException {
        Item c = constructsByFbtp.get(fbtp);
        if (c == null) {
            c = createItem("TransgenicConstruct");
            c.setAttribute("primaryIdentifier", fbtp);
            c.setReference("organism", getOrganism());
            constructsByFbtp.put(fbtp, c);
        }
        return c;
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
        LOG.info("flybase-transgenic-constructs: close(): storing "
                 + constructsByFbtp.size() + " constructs, "
                 + genesByFbgn.size() + " gene stubs, "
                 + insertionsByFbti.size() + " insertion stubs");
        for (Item c : constructsByFbtp.values()) {
            store(c);
        }
        for (Item g : genesByFbgn.values()) {
            store(g);
        }
        for (Item i : insertionsByFbti.values()) {
            store(i);
        }
        super.close();
    }
}
