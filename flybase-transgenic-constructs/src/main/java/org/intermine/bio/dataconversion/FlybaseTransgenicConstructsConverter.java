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
 * REWRITTEN 2026-06-05: the actual FB2026_01 file
 * (transgenic_construct_descriptions_*.tsv) is 17 columns and
 * allele-centric — one row per (component allele, transgenic construct)
 * pair. The round 3 converter assumed a 6-col construct-centric layout
 * and silently skipped every row (defence-in-depth filter
 * name.startsWith("FBtp") on col 1, which actually holds an FBal symbol).
 *
 * Input filename: transgenic_construct_descriptions_*.tsv (header lines
 * start with `##` or `#`).
 *
 * Columns (tab-separated, 1-indexed):
 *    1 Component Allele (symbol)            — text
 *    2 Component Allele (id)                — FBal
 *    3 Transgenic Construct (symbol)        — text
 *    4 Transgenic Construct (id)            — FBtp  [group key]
 *    5 Transgenic Product class (term)      — text -> productClass
 *    6 Transgenic Product class (id)        — FBcv (currently ignored)
 *    7 Regulatory region (symbol)           — text -> regulatoryRegion
 *    8 Regulatory region (id)               — FBgn (currently ignored)
 *    9 Encoded product/tool (symbol)        — text -> encodedProduct
 *   10 Encoded product/tool (id)            — FBgn (currently ignored)
 *   11 Tagged with (symbol)                 — text -> taggedWith
 *   12 Tagged with (id)                     — FBgn (currently ignored)
 *   13 Also carries (symbol)                — text -> alsoCarries
 *   14 Also carries (id)                    — FBgn list (semicolon-separated)
 *   15 Description (text)                   — text -> description
 *   16 Description (supporting reference)   — FBrf (currently ignored)
 *   17 Stocks (number)                      — int -> stockCount
 *
 * Grouping strategy: a TransgenicConstruct Item is created on first
 * sight of an FBtp (col 4). Scalar attributes use first-non-empty-wins
 * across the rows for that FBtp. The componentAlleles collection
 * accumulates one ref per row's FBal (col 2). The alsoCarriesGenes
 * collection accumulates semicolon-split FBgn refs from col 14.
 *
 * @author new_flymine session, 2026-06-04 (rewrite 2026-06-05)
 */
public class FlybaseTransgenicConstructsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseTransgenicConstructsConverter.class);

    private static final String DATASET_TITLE = "FlyBase transgenic constructs";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    // 0-indexed
    private static final int COL_ALLELE_FBAL    = 1;
    private static final int COL_CONSTRUCT_SYM  = 2;
    private static final int COL_CONSTRUCT_FBTP = 3;
    private static final int COL_PRODUCT_CLASS  = 4;
    private static final int COL_REG_REGION_SYM = 6;
    private static final int COL_ENCODED_SYM    = 8;
    private static final int COL_TAGGED_WITH    = 10;
    private static final int COL_ALSO_CARRIES   = 12;
    private static final int COL_ALSO_CARRIES_FBGN = 13;
    private static final int COL_DESCRIPTION    = 14;
    private static final int COL_STOCK_COUNT    = 16;
    private static final int MIN_COLS = 4;   // need at least up to FBtp

    /** Semicolon possibly followed by whitespace. */
    private static final Pattern SEMI_SPLIT = Pattern.compile(";\\s*");

    private final Map<String, Item> constructsByFbtp = new HashMap<String, Item>();
    private final Map<String, Item> allelesByFbal = new HashMap<String, Item>();
    private final Map<String, Item> genesByFbgn = new HashMap<String, Item>();
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
            LOG.info("flybase-transgenic-constructs: not a transgenic_construct_descriptions"
                     + " file (" + name + "), skipping");
            return;
        }
        LOG.info("flybase-transgenic-constructs: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int rows = 0;
        int alleleAdds = 0;
        int alsoCarriesAdds = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < MIN_COLS) {
                continue;
            }
            String fbtp = cols[COL_CONSTRUCT_FBTP].trim();
            if (fbtp.isEmpty() || !fbtp.startsWith("FBtp")) {
                continue;
            }
            Item construct = getOrCreateConstruct(fbtp);

            // construct symbol — set on first non-empty
            if (cols.length > COL_CONSTRUCT_SYM) {
                setAttributeIfMissing(construct, "symbol", cols[COL_CONSTRUCT_SYM]);
            }
            // scalar attrs, first-non-empty-wins
            if (cols.length > COL_PRODUCT_CLASS) {
                setAttributeIfMissing(construct, "productClass", cols[COL_PRODUCT_CLASS]);
            }
            if (cols.length > COL_REG_REGION_SYM) {
                setAttributeIfMissing(construct, "regulatoryRegion", cols[COL_REG_REGION_SYM]);
            }
            if (cols.length > COL_ENCODED_SYM) {
                setAttributeIfMissing(construct, "encodedProduct", cols[COL_ENCODED_SYM]);
            }
            if (cols.length > COL_TAGGED_WITH) {
                setAttributeIfMissing(construct, "taggedWith", cols[COL_TAGGED_WITH]);
            }
            if (cols.length > COL_ALSO_CARRIES) {
                setAttributeIfMissing(construct, "alsoCarries", cols[COL_ALSO_CARRIES]);
            }
            if (cols.length > COL_DESCRIPTION) {
                setAttributeIfMissing(construct, "description", cols[COL_DESCRIPTION]);
            }
            if (cols.length > COL_STOCK_COUNT) {
                String stockCount = cols[COL_STOCK_COUNT].trim();
                if (!stockCount.isEmpty() && construct.getAttribute("stockCount") == null) {
                    try {
                        Integer.parseInt(stockCount);   // validate
                        construct.setAttribute("stockCount", stockCount);
                    } catch (NumberFormatException nfe) {
                        // ignore non-numeric stockCount
                    }
                }
            }

            // componentAlleles: one Allele ref per row's FBal (col 2)
            if (cols.length > COL_ALLELE_FBAL) {
                String fbal = cols[COL_ALLELE_FBAL].trim();
                if (!fbal.isEmpty() && fbal.startsWith("FBal")) {
                    Item allele = getOrCreateAllele(fbal);
                    construct.addToCollection("componentAlleles", allele);
                    alleleAdds++;
                }
            }

            // alsoCarriesGenes: semicolon-split FBgn list from col 14
            if (cols.length > COL_ALSO_CARRIES_FBGN) {
                String fbgnField = cols[COL_ALSO_CARRIES_FBGN];
                if (!fbgnField.isEmpty()) {
                    for (String fbgn : SEMI_SPLIT.split(fbgnField)) {
                        String trimmed = fbgn.trim();
                        if (trimmed.isEmpty() || !trimmed.startsWith("FBgn")) {
                            continue;
                        }
                        Item gene = getOrCreateGene(trimmed);
                        construct.addToCollection("alsoCarriesGenes", gene);
                        alsoCarriesAdds++;
                    }
                }
            }

            rows++;
        }
        LOG.info("flybase-transgenic-constructs: rows=" + rows
                 + " distinct FBtp=" + constructsByFbtp.size()
                 + " componentAllele adds=" + alleleAdds
                 + " alsoCarriesGenes adds=" + alsoCarriesAdds);
    }

    private void setAttributeIfMissing(Item item, String attrName, String value) {
        if (value == null) {
            return;
        }
        String trimmed = value.trim();
        if (trimmed.isEmpty()) {
            return;
        }
        if (item.getAttribute(attrName) != null) {
            return;
        }
        item.setAttribute(attrName, trimmed);
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
        LOG.info("flybase-transgenic-constructs: close(): storing "
                 + constructsByFbtp.size() + " constructs, "
                 + allelesByFbal.size() + " allele stubs, "
                 + genesByFbgn.size() + " gene stubs");
        for (Item c : constructsByFbtp.values()) {
            store(c);
        }
        for (Item a : allelesByFbal.values()) {
            store(a);
        }
        for (Item g : genesByFbgn.values()) {
            store(g);
        }
        super.close();
    }
}
