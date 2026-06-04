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
 */
public class FlybaseTransgenicConstructsConverterTest extends ItemsTestCase
{
    public FlybaseTransgenicConstructsConverterTest(String name) {
        super(name);
    }

    /**
     * 4 FBtp data rows in fixture all stored (lenient cols-check).
     * Header lines skipped.
     */
    public void testConstructsCreatedAndHeaderRowsSkipped() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
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
        // FBtp0000003 has only 3 cols -> >= MIN_COLS (4) is false -> skipped
        // FBtp0000004 has 2 cols -> skipped
        // FBtp0000001 + FBtp0000002 stored
        assertEquals(2, constructs);
    }

    /**
     * Semicolon-space-separated FBgn list explodes into multiple Gene refs;
     * each Gene merged-by-FBgn.
     */
    public void testEncodedGenesParsed() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("transgenic_construct_descriptions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        // FBtp0000001 has 2 FBgns -> 2 encodedGenes
        Item construct1 = items.stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBtp0000001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBtp0000001 not stored"));
        long encoded = construct1.getCollections().stream()
            .filter(col -> "encodedGenes".equals(col.getName()))
            .mapToLong(col -> col.getRefIds().size()).sum();
        assertEquals(2, encoded);
    }

    /**
     * Semicolon-space-separated FBti list explodes into multiple
     * Insertion refs.
     */
    public void testRelatedInsertionsParsed() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("transgenic_construct_descriptions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        // FBtp0000001 has 2 FBtis -> 2 relatedInsertions
        Item construct1 = items.stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBtp0000001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBtp0000001 not stored"));
        long related = construct1.getCollections().stream()
            .filter(col -> "relatedInsertions".equals(col.getName()))
            .mapToLong(col -> col.getRefIds().size()).sum();
        assertEquals(2, related);
    }

    /**
     * Filename dispatch defence-in-depth: a non-transgenic_construct file
     * is skipped.
     */
    public void testNonConstructsFileSkipped() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseTransgenicConstructsConverter c =
            new FlybaseTransgenicConstructsConverter(
                writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("fbal_to_fbgn_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/transgenic_construct_descriptions_test.tsv")));
        c.close();
        long constructs = writer.getItems().stream()
            .filter(i -> "TransgenicConstruct".equals(i.getClassName()))
            .count();
        assertEquals(0, constructs);
    }
}
