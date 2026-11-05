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
 * Loader for FlyBase gene external cross-references.
 *
 * Reads one staged 3-column TSV from src.data.dir:
 *   gene_dbxrefs.tsv
 *
 *   col 0 = FBgn        -> Gene.primaryIdentifier (merge key)
 *   col 1 = db_name     -> DataSource.name
 *   col 2 = accession   -> CrossReference.identifier
 *
 * The orchestrator generates this file by extracting FlyBase gene
 * feature_dbxref from chado with a db allow-list; this converter only reads
 * the 3-column format (## / blank lines skipped).
 *
 * For each row it merges the Gene by FBgn (Gene.key = primaryIdentifier),
 * gets-or-creates a DataSource keyed by db_name, and creates a CrossReference
 * (identifier=accession, source=DataSource) attached to the Gene via
 * CrossReference.subject (reverse-reference Gene.crossReferences). Rows are
 * deduped on (FBgn, db_name, accession).
 *
 * @author flymine builder session, 2026-09
 */
public class FlybaseGeneCrossReferencesConverter extends BioFileConverter
{
    private static final Logger LOG =
        Logger.getLogger(FlybaseGeneCrossReferencesConverter.class);

    private static final String DATASET_TITLE = "FlyBase gene cross references";
    private static final String DATA_SOURCE_NAME = "FlyBase";
    private static final String DMEL_TAXON = "7227";

    private static final int COL_FBGN = 0;
    private static final int COL_DB = 1;
    private static final int COL_ACCESSION = 2;
    private static final int MIN_COLS = 3;

    private final Map<String, Item> genesByFbgn = new HashMap<String, Item>();
    private final Map<String, Item> dataSourcesByName = new HashMap<String, Item>();
    private final Set<String> seen = new HashSet<String>();
    private Item organism;

    /**
     * Constructor.
     * @param writer the ItemWriter
     * @param model the target data model
     */
    public FlybaseGeneCrossReferencesConverter(ItemWriter writer, Model model) {
        super(writer, model, DATA_SOURCE_NAME, DATASET_TITLE);
    }

    @Override
    public void process(Reader reader) throws Exception {
        String name = getCurrentFile().getName();
        LOG.info("flybase-gene-crossreferences: processing " + name);

        BufferedReader br = new BufferedReader(reader);
        String line;
        int xrefs = 0;
        while ((line = br.readLine()) != null) {
            if (line.isEmpty() || line.startsWith("#")) {
                continue;
            }
            String[] cols = line.split("\\t", -1);
            if (cols.length < MIN_COLS) {
                continue;
            }
            String fbgn = cols[COL_FBGN].trim();
            String db = cols[COL_DB].trim();
            String accession = cols[COL_ACCESSION].trim();
            if (fbgn.isEmpty() || db.isEmpty() || accession.isEmpty()) {
                continue;
            }
            String dedupKey = fbgn + "\t" + db + "\t" + accession;
            if (!seen.add(dedupKey)) {
                continue;
            }
            Item gene = getOrCreateGene(fbgn);
            Item source = getOrCreateDataSource(db);
            Item xref = createItem("CrossReference");
            xref.setAttribute("identifier", accession);
            xref.setReference("source", source);
            xref.setReference("subject", gene);
            store(xref);
            xrefs++;
        }
        LOG.info("flybase-gene-crossreferences: crossReferences=" + xrefs);
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

    private Item getOrCreateDataSource(String dbName) throws ObjectStoreException {
        Item ds = dataSourcesByName.get(dbName);
        if (ds == null) {
            ds = createItem("DataSource");
            ds.setAttribute("name", dbName);
            store(ds);
            dataSourcesByName.put(dbName, ds);
        }
        return ds;
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
        LOG.info("flybase-gene-crossreferences: close(): storing "
                 + genesByFbgn.size() + " genes, "
                 + dataSourcesByName.size() + " data sources");
        for (Item g : genesByFbgn.values()) {
            store(g);
        }
        super.close();
    }
}
