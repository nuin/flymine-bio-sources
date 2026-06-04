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
 * Unit tests for FlybaseAlleleStubsConverter.
 */
public class FlybaseAlleleStubsConverterTest extends ItemsTestCase
{
    public FlybaseAlleleStubsConverterTest(String name) {
        super(name);
    }

    /**
     * 4 valid FBal rows in the fixture, 1 header block, 1 malformed (single
     * token), 1 with empty symbol -> 4 Allele items, 4 Gene refs (by FBgn).
     * The InterMine loader's merge-by-primaryIdentifier handles the
     * "already in chado" vs "new stub" distinction; the converter just
     * emits items for every valid row.
     */
    public void testStubsCreatedFromTsv() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAlleleStubsConverter c =
            new FlybaseAlleleStubsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("fbal_to_fbgn_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fbal_to_fbgn_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        long alleles = items.stream()
            .filter(i -> "Allele".equals(i.getClassName())).count();
        long genes = items.stream()
            .filter(i -> "Gene".equals(i.getClassName())).count();
        // 4 valid FBal rows; 4 distinct FBgns
        assertEquals(4, alleles);
        assertEquals(4, genes);
    }

    /**
     * The FBal with empty symbol should still produce an Allele item
     * (primaryIdentifier + gene ref + organism), just without symbol set.
     */
    public void testEmptySymbolStillCreatesAllele() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAlleleStubsConverter c =
            new FlybaseAlleleStubsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("fbal_to_fbgn_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fbal_to_fbgn_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        Item noSymbolAllele = items.stream()
            .filter(i -> "Allele".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBal0000300".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBal0000300 not stored"));
        boolean hasSymbolAttr = noSymbolAllele.getAttributes().stream()
            .anyMatch(a -> "symbol".equals(a.getName()));
        assertFalse("FBal0000300 has empty symbol col; symbol must not be set",
                    hasSymbolAttr);
    }

    /**
     * Gene reference matches the FBgn column.
     */
    public void testGeneRefByFBgn() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAlleleStubsConverter c =
            new FlybaseAlleleStubsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("fbal_to_fbgn_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fbal_to_fbgn_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        // The Ecol\lacZ allele FBal0000099 has FBgn0014447 in col 3
        Item lacZ = items.stream()
            .filter(i -> "Allele".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBal0000099".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBal0000099 not stored"));
        String geneRefId = lacZ.getReferences().stream()
            .filter(r -> "gene".equals(r.getName()))
            .map(r -> r.getRefId())
            .findFirst()
            .orElseThrow(() -> new AssertionError("lacZ has no gene ref"));
        // The refId points to the Gene item with primaryIdentifier=FBgn0014447
        Item gene = items.stream()
            .filter(i -> "Gene".equals(i.getClassName()))
            .filter(i -> geneRefId.equals(i.getIdentifier()))
            .findFirst()
            .orElseThrow(() -> new AssertionError("gene refId doesn't match any Gene"));
        boolean hasFbgn = gene.getAttributes().stream().anyMatch(a ->
            "primaryIdentifier".equals(a.getName())
            && "FBgn0014447".equals(a.getValue()));
        assertTrue("gene ref must point at FBgn0014447", hasFbgn);
    }

    /**
     * Filename dispatch defence-in-depth: a non-fbal_to_fbgn file (e.g.
     * the descriptions TSV) must be skipped without parsing.
     */
    public void testNonStubsFileSkipped() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAlleleStubsConverter c =
            new FlybaseAlleleStubsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("dmel_classical_and_insertion_allele_descriptions_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fbal_to_fbgn_test.tsv")));
        c.close();
        long alleles = writer.getItems().stream()
            .filter(i -> "Allele".equals(i.getClassName())).count();
        assertEquals(0, alleles);
    }
}
