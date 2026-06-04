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
 * Loader for FlyBase classical + insertion allele descriptions (v1).
 *
 * Reads exactly one file from src.data.dir:
 *   dmel_classical_and_insertion_allele_descriptions_*.tsv  (21 columns)
 *
 * The other files that sit in the same FB precomputed alleles/ subdir
 * (allele_genetic_interactions, fbal_to_fbgn, genotype_phenotype_data,
 * split_system_combinations) have different column counts; the source's
 * src.data.dir.includes property in project.xml narrows iteration to just
 * the descriptions file. As a defence-in-depth, this converter also
 * dispatches by filename so a misconfigured includes filter does not
 * misparse a sibling.
 *
 * v1 scope: extract description text + supporting publication FBrf refs
 * and merge into existing Allele items (chado-db-flybase-dmel already
 * created the Allele rows). All other 21-column fields (insertion details,
 * encoded product, tagged with, "also carries", alleleClass multi-values)
 * are deferred to v2.
 *
 * Columns used:
 *   col 2 (index 1): Allele FBal id    -> Allele.primaryIdentifier
 *   col 19 (index 18): Description text -> Allele.description
 *   col 20 (index 19): pipe-separated FBrf ids -> Allele.descriptionReferences
 *
 * Skips rows with empty description (most rows: builder reported only a
 * minority of FBal entries carry a description in FB2026_01).
 *
 * @author new_flymine session, 2026-06-03
 */
public class FlybaseAlleleDescriptionsConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseAlleleDescriptionsConverter.class);

    private static final String DATASET_TITLE = "FlyBase classical and insertion allele descriptions";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    private static final int COL_FBAL = 1;
    private static final int COL_DESCRIPTION = 18;
    private static final int COL_REFERENCES = 19;
    private static final int MIN_COLS = 20;

    private final Map<String, Item> allelesByFbal = new HashMap<String, Item>();
    private final Map<String, Item> publicationsByFbrf = new HashMap<String, Item>();
    private Item organism;

    /**
     * Constructor.
     * @param writer the ItemWriter
     * @param model the target data model
     */
    public FlybaseAlleleDescriptionsConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        if (name.endsWith(".gz")) {
            LOG.info("flybase-allele-descriptions: skipping gz companion " + name);
            return;
        }
        if (!name.startsWith("dmel_classical_and_insertion_allele_descriptions")) {
            LOG.info("flybase-allele-descriptions: not a descriptions file ("
                     + name + "), skipping");
            return;
        }
        LOG.info("flybase-allele-descriptions: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int rowsTotal = 0;
        int rowsWithDescription = 0;
        int refsTotal = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < MIN_COLS) {
                continue;
            }
            rowsTotal++;
            String fbal = cols[COL_FBAL];
            String description = cols[COL_DESCRIPTION];
            if (description.isEmpty()) {
                continue;
            }
            rowsWithDescription++;
            Item allele = getOrCreateAllele(fbal);
            allele.setAttribute("description", description);

            if (cols.length > COL_REFERENCES) {
                String refsField = cols[COL_REFERENCES];
                if (!refsField.isEmpty()) {
                    for (String fbrf : refsField.split("\\|")) {
                        String trimmed = fbrf.trim();
                        if (trimmed.isEmpty() || !trimmed.startsWith("FBrf")) {
                            continue;
                        }
                        Item pub = getOrCreatePublication(trimmed);
                        allele.addToCollection("descriptionReferences", pub);
                        refsTotal++;
                    }
                }
            }
        }
        LOG.info("flybase-allele-descriptions: rowsTotal=" + rowsTotal
                 + " rowsWithDescription=" + rowsWithDescription
                 + " refsTotal=" + refsTotal);
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

    private Item getOrCreatePublication(String fbrf) throws ObjectStoreException {
        Item p = publicationsByFbrf.get(fbrf);
        if (p == null) {
            p = createItem("Publication");
            p.setAttribute("primaryIdentifier", fbrf);
            publicationsByFbrf.put(fbrf, p);
        }
        return p;
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
        LOG.info("flybase-allele-descriptions: close(): storing "
                 + allelesByFbal.size() + " alleles, "
                 + publicationsByFbrf.size() + " publications");
        for (Item a : allelesByFbal.values()) {
            store(a);
        }
        for (Item p : publicationsByFbrf.values()) {
            store(p);
        }
        super.close();
    }
}
