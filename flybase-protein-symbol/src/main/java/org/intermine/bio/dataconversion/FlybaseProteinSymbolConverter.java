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
 * Loader for FlyBase polypeptide (protein) symbols.
 *
 * Reads one file from src.data.dir:
 *   fbgn_fbtr_fbpp_expanded_*.tsv
 *
 * Header (leading "## " comment line):
 *   organism gene_type gene_ID gene_symbol gene_fullname annotation_ID
 *   transcript_type transcript_ID transcript_symbol polypeptide_ID
 *   polypeptide_symbol
 *
 *   col  9 = polypeptide_ID     (FBpp) -> Protein.secondaryIdentifier (merge key)
 *   col 10 = polypeptide_symbol         -> Protein.symbol
 *
 * chado-db-flybase-dmel already set Protein.secondaryIdentifier = FBpp
 * uniquename, so this source merges by secondaryIdentifier (see
 * keys.properties) and only adds the symbol. Rows without an FBpp
 * (ncRNA / tRNA genes) are skipped. Each FBpp repeats across transcript
 * rows sharing the polypeptide; the first occurrence wins.
 *
 * The file uses DOS line endings (CRLF); a trailing carriage return on the
 * last column is stripped.
 *
 * @author flymine builder session, 2026-09
 */
public class FlybaseProteinSymbolConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseProteinSymbolConverter.class);

    private static final String DATASET_TITLE = "FlyBase protein symbols";
    private static final String DATA_SOURCE_NAME = "FlyBase";

    private static final int COL_FBPP = 9;
    private static final int COL_SYMBOL = 10;
    private static final int MIN_COLS = 11;

    private final Map<String, Item> proteinsByFbpp = new HashMap<String, Item>();

    /**
     * Constructor.
     * @param writer the ItemWriter
     * @param model the target data model
     */
    public FlybaseProteinSymbolConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-protein-symbol: skipping gz companion " + name);
            return;
        }
        LOG.info("flybase-protein-symbol: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int rows = 0;
        while ((line = br.readLine()) != null) {
            // strip any trailing CR from CRLF line endings
            if (line.endsWith("\r")) {
                line = line.substring(0, line.length() - 1);
            }
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < MIN_COLS) {
                continue;
            }
            String fbpp = cols[COL_FBPP].trim();
            String symbol = cols[COL_SYMBOL].trim();
            if (fbpp.isEmpty() || symbol.isEmpty()) {
                continue;
            }
            if (proteinsByFbpp.containsKey(fbpp)) {
                continue;
            }
            Item protein = createItem("Protein");
            protein.setAttribute("secondaryIdentifier", fbpp);
            protein.setAttribute("symbol", symbol);
            proteinsByFbpp.put(fbpp, protein);
            rows++;
        }
        LOG.info("flybase-protein-symbol: distinct FBpp with symbol=" + rows);
    }

    @Override
    public void close() throws Exception {
        LOG.info("flybase-protein-symbol: close(): storing "
                 + proteinsByFbpp.size() + " proteins");
        for (Item p : proteinsByFbpp.values()) {
            store(p);
        }
        super.close();
    }
}
