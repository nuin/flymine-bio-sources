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

import java.io.File;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.LinkedHashMap;

import org.intermine.dataconversion.ItemsTestCase;
import org.intermine.dataconversion.MockItemWriter;
import org.intermine.metadata.Model;
import org.intermine.model.fulldata.Item;

/**
 * Unit tests for FlybaseTransgenicConstructsConverter.
 *
 * Rewritten 2026-06-05 to match the actual 17-col allele-centric file
 * format. Fixture has:
 *   - 2 rows for FBtp0000001 (P{lacW}) — components lacZ[T125] + white[mini]
 *   - 1 row for FBtp0000002 (P{Gal4-sd}) — component GAL4[sd-SG29.1]
 *   - 1 noise row (col 4 not an FBtp)
 *   - 1 malformed row (too few cols)
 * Expected: 2 distinct TransgenicConstruct items.
 */
public class FlybaseTransgenicConstructsConverterTest extends ItemsTestCase
{
    public FlybaseTransgenicConstructsConverterTest(String name) {
        super(name);
    }

    /** 4 well-formed FBtp rows aggregate to 2 distinct TransgenicConstruct items. */
    public void testConstructsCreatedAndGroupedByFBtp() throws Exception {
        MockItemWriter writer = new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("transgenic_construct_descriptions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        long constructs = items.stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName()))
            .count();
        assertEquals(2, constructs);
    }

    /**
     * FBtp0000001 has 2 component-allele rows -> componentAlleles
     * collection size 2; FBtp0000002 has 1 -> size 1.
     */
    public void testComponentAllelesCollection() throws Exception {
        MockItemWriter writer = new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("transgenic_construct_descriptions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        Item lacW = items.stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBtp0000001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBtp0000001 not stored"));
        long componentCount = lacW.getCollections().stream()
            .filter(col -> "componentAlleles".equals(col.getName()))
            .mapToLong(col -> col.getRefIds().size()).sum();
        assertEquals(2, componentCount);
    }

    /** FBtp0000001 row 1 col 14: "FBgn0003996; FBgn0003997" -> 2 alsoCarriesGenes. */
    public void testAlsoCarriesGenesParsed() throws Exception {
        MockItemWriter writer = new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("transgenic_construct_descriptions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        Item lacW = items.stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBtp0000001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBtp0000001 not stored"));
        long alsoCarriesCount = lacW.getCollections().stream()
            .filter(col -> "alsoCarriesGenes".equals(col.getName()))
            .mapToLong(col -> col.getRefIds().size()).sum();
        assertEquals(2, alsoCarriesCount);
    }

    /**
     * Scalar attrs use first-non-empty-wins. FBtp0000001 row 1 has
     * description="Classic lacZ reporter construct"; row 2 has it blank.
     * The Item must carry row 1's value (not be overwritten or cleared).
     */
    public void testDescriptionWinsFirstNonEmpty() throws Exception {
        MockItemWriter writer = new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("transgenic_construct_descriptions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        Item lacW = items.stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBtp0000001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBtp0000001 not stored"));
        boolean hasDesc = lacW.getAttributes().stream().anyMatch(a ->
            "description".equals(a.getName())
            && "Classic lacZ reporter construct".equals(a.getValue()));
        assertTrue("FBtp0000001 must keep row-1's description across rows", hasDesc);

        // stockCount also set from row 1 (=5)
        boolean hasStockCount = lacW.getAttributes().stream().anyMatch(a ->
            "stockCount".equals(a.getName()) && "5".equals(a.getValue()));
        assertTrue("FBtp0000001 stockCount=5", hasStockCount);
    }

    /** Header rows (## / #) and malformed rows skipped. */
    public void testHeaderAndMalformedRowsSkipped() throws Exception {
        MockItemWriter writer = new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("transgenic_construct_descriptions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();
        // 3 well-formed FBtp rows + 1 noise (col 4 not FBtp) + 1 malformed
        // -> 2 distinct TransgenicConstruct items, 3 distinct FBal items
        Collection<Item> items = writer.getItems();
        long constructs = items.stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName())).count();
        long alleles = items.stream()
            .filter(i -> "Allele".equals(i.getClassName())).count();
        assertEquals(2, constructs);
        // FBal0000099, FBal0000100, FBal0000200 — FBal0000300 not picked up
        // because its row's col 4 is "FBgn0040000" not an FBtp.
        assertEquals(3, alleles);
    }

    /** Filename dispatch defence-in-depth. */
    public void testNonConstructsFileSkipped() throws Exception {
        MockItemWriter writer = new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("fbal_to_fbgn_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();
        long constructs = writer.getItems().stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName())).count();
        assertEquals(0, constructs);
    }
}
