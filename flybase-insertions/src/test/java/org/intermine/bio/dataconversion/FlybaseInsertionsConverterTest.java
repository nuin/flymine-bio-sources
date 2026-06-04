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
 * Unit tests for FlybaseInsertionsConverter.
 */
public class FlybaseInsertionsConverterTest extends ItemsTestCase
{
    public FlybaseInsertionsConverterTest(String name) {
        super(name);
    }

    /**
     * Fixture has 4 data rows: FBti0000001/2 have full location, FBti0000003
     * has only an FBti id, FBti0000004 has location but no strand. All four
     * are valid Insertion items. Header rows (##) skipped.
     */
    public void testInsertionsCreatedAndHeaderRowsSkipped() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseInsertionsConverter c =
            new FlybaseInsertionsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("insertion_mapping_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/insertion_mapping_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        long inserts = items.stream()
            .filter(i -> "TransposableElementInsertionSite".equals(i.getClassName()))
            .count();
        // FBti0000003 has only 1 col -> below MIN_COLS, skipped
        // FBti0000004 has 5 cols -> also below MIN_COLS (need 7), skipped
        // FBti0000001 + FBti0000002 have full 7 cols -> stored
        assertEquals(2, inserts);
    }

    /**
     * Genomic location parser: "2L:7,238,500..7,238,500" yields
     * chromosome=2L, start=end=7238500.
     */
    public void testChromosomeLocationParsed() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseInsertionsConverter c =
            new FlybaseInsertionsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("insertion_mapping_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/insertion_mapping_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        long locations = items.stream()
            .filter(i -> "Location".equals(i.getClassName())).count();
        assertEquals(2, locations);

        Item firstLoc = items.stream()
            .filter(i -> "Location".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "start".equals(a.getName()) && "7238500".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("expected Location with start=7238500"));
        boolean hasEnd = firstLoc.getAttributes().stream().anyMatch(a ->
            "end".equals(a.getName()) && "7238500".equals(a.getValue()));
        assertTrue("expected matching end=7238500", hasEnd);
    }

    /**
     * Orientation column -> Location.strand: "+" yields "1", "-" yields "-1".
     */
    public void testOrientationStrand() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseInsertionsConverter c =
            new FlybaseInsertionsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("insertion_mapping_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/insertion_mapping_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        // FBti0000001 (+) -> strand="1"; FBti0000002 (-) -> strand="-1"
        long plus = items.stream()
            .filter(i -> "Location".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "strand".equals(a.getName()) && "1".equals(a.getValue())))
            .count();
        long minus = items.stream()
            .filter(i -> "Location".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "strand".equals(a.getName()) && "-1".equals(a.getValue())))
            .count();
        assertEquals(1, plus);
        assertEquals(1, minus);
    }

    /**
     * Filename dispatch defence-in-depth: a non-insertion_mapping file must
     * be skipped without parsing.
     */
    public void testNonInsertionsFileSkipped() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseInsertionsConverter c =
            new FlybaseInsertionsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("fbal_to_fbgn_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/insertion_mapping_test.tsv")));
        c.close();
        long inserts = writer.getItems().stream()
            .filter(i -> "TransposableElementInsertionSite".equals(i.getClassName()))
            .count();
        assertEquals(0, inserts);
    }
}
