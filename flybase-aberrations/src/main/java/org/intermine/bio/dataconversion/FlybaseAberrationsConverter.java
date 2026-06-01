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

import java.io.Reader;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
 * @author new_flymine session, 2026-06-01
 */
public class FlybaseAberrationsConverter extends BioFileConverter
{
    private static final String DATASET_TITLE = "FlyBase aberrations and balancers";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    private static final Pattern ABERRATION_PREFIX =
        Pattern.compile("^(Df|Dp|In|T)\\(.*\\).*");

    private final Map<String, Item> aberrationsById = new HashMap<String, Item>();
    private final Map<String, Item> balancersById   = new HashMap<String, Item>();
    private final Map<String, Item> genesByFbgn     = new HashMap<String, Item>();
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
     * Dispatch by input filename.
     * @param reader reader over the current input file
     */
    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.startsWith("fb_synonym")) {
            processSynonyms(reader);
        } else if (name.startsWith("aberration_experimental_gene_del_dup_data")) {
            processDelDup(reader);
        } else if (name.startsWith("fbba_to_fbab")) {
            processCuratedBalancers(reader);
        }
        // unknown files are silently skipped
    }

    void processSynonyms(Reader reader) throws Exception {
        throw new UnsupportedOperationException("not yet implemented");
    }

    void processDelDup(Reader reader) throws Exception {
        throw new UnsupportedOperationException("not yet implemented");
    }

    void processCuratedBalancers(Reader reader) throws Exception {
        throw new UnsupportedOperationException("not yet implemented");
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
}
