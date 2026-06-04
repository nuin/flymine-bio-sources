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
import java.io.Reader;
import java.util.Collection;
import java.util.LinkedHashMap;

import org.intermine.dataconversion.ItemsTestCase;
import org.intermine.dataconversion.MockItemWriter;
import org.intermine.metadata.Model;
import org.intermine.model.fulldata.Item;

/**
 * Unit tests for FlybaseAlleleDescriptionsConverter.
 */
public class FlybaseAlleleDescriptionsConverterTest extends ItemsTestCase
{
    public FlybaseAlleleDescriptionsConverterTest(String name) {
        super(name);
    }

    /**
     * Description present + FBrf refs wired; empty-description rows skipped;
     * pipe-separated FBrf list explodes into multiple Publication items.
     */
    public void testDescriptionsAndReferencesWired() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAlleleDescriptionsConverter c =
            new FlybaseAlleleDescriptionsConverter(writer, Model.getInstanceByName("genomic"));
        Reader r = new InputStreamReader(getClass().getResourceAsStream(
            "/dmel_classical_and_insertion_allele_descriptions_test.tsv"));
        c.setCurrentFile(new File(
            "dmel_classical_and_insertion_allele_descriptions_test.tsv"));
        c.process(r);
        c.close();

        Collection<Item> items = writer.getItems();

        // 3 alleles with non-empty descriptions: FBal0000008, FBal0000099, FBal0000100
        // (FBal0000001 has empty description -> skipped, no Allele created)
        long alleleCount = items.stream()
            .filter(i -> "Allele".equals(i.getClassName())).count();
        assertEquals(3, alleleCount);

        // FBal0000099 has 3 FBrf refs (Multi-reference description)
        Item multiRefAllele = items.stream()
            .filter(i -> "Allele".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBal0000099".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBal0000099 not stored"));
        long refCount = multiRefAllele.getCollections().stream()
            .filter(c2 -> "descriptionReferences".equals(c2.getName()))
            .mapToLong(c2 -> c2.getRefIds().size()).sum();
        assertEquals(3, refCount);

        // FBal0000100 has description but empty references -> 0 refs
        Item noRefAllele = items.stream()
            .filter(i -> "Allele".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBal0000100".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBal0000100 not stored"));
        long emptyRefCount = noRefAllele.getCollections().stream()
            .filter(c2 -> "descriptionReferences".equals(c2.getName()))
            .mapToLong(c2 -> c2.getRefIds().size()).sum();
        assertEquals(0, emptyRefCount);

        // 4 total Publication items: FBrf0000001 + FBrf0000010 + FBrf0000020 + FBrf0000030
        long pubCount = items.stream()
            .filter(i -> "Publication".equals(i.getClassName())).count();
        assertEquals(4, pubCount);
    }

    /**
     * A non-descriptions file (e.g. allele_genetic_interactions) must be
     * skipped without parsing — defence-in-depth for misconfigured
     * src.data.dir.includes.
     */
    public void testNonDescriptionsFileSkipped() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAlleleDescriptionsConverter c =
            new FlybaseAlleleDescriptionsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("allele_genetic_interactions_test.tsv"));
        c.process(new InputStreamReader(getClass().getResourceAsStream(
            "/dmel_classical_and_insertion_allele_descriptions_test.tsv")));
        c.close();
        // converter saw a non-descriptions filename: should bail out, no items stored
        long alleles = writer.getItems().stream()
            .filter(i -> "Allele".equals(i.getClassName())).count();
        assertEquals(0, alleles);
    }
}
