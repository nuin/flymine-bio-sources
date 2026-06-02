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
 * Unit tests for FlybaseAberrationsConverter.
 * Tests for the three input branches accumulate as Tasks 7-10 land.
 */
public class FlybaseAberrationsConverterTest extends ItemsTestCase
{
    public FlybaseAberrationsConverterTest(String name) {
        super(name);
    }

    /** Task 7 — static prefix parser. */
    public void testAberrationTypeFromSymbol() {
        assertEquals("deletion",
            FlybaseAberrationsConverter.aberrationTypeFromSymbol("Df(2R)min"));
        assertEquals("duplication",
            FlybaseAberrationsConverter.aberrationTypeFromSymbol("Dp(1;2)test"));
        assertEquals("inversion",
            FlybaseAberrationsConverter.aberrationTypeFromSymbol("In(2LR)a[M60]"));
        assertEquals("translocation",
            FlybaseAberrationsConverter.aberrationTypeFromSymbol("T(2;3)foo"));
        assertEquals("other",
            FlybaseAberrationsConverter.aberrationTypeFromSymbol("Mystery(test)"));
        assertEquals("other",
            FlybaseAberrationsConverter.aberrationTypeFromSymbol(null));
    }

    /** Task 8 — synonyms file creates Aberration + Balancer items. */
    public void testSynonymsCreatesAberrationsAndBalancers() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAberrationsConverter c =
            new FlybaseAberrationsConverter(writer, Model.getInstanceByName("genomic"));
        Reader r = new InputStreamReader(
            getClass().getResourceAsStream("/fb_synonym_test.tsv"));
        c.setCurrentFile(new File("fb_synonym_test.tsv"));
        c.process(r);
        c.close();

        Collection<Item> items = writer.getItems();
        long aberrations = items.stream()
            .filter(i -> "Aberration".equals(i.getClassName())).count();
        long balancers = items.stream()
            .filter(i -> "Balancer".equals(i.getClassName())).count();
        assertEquals(5, aberrations);
        assertEquals(2, balancers);
    }

    /**
     * Task 9 — del/dup file wires deletedGenes / duplicatedGenes;
     * skips "not deleted" / "not duplicated" negative rows.
     */
    public void testDelDupSkipsNegativesAndWiresGenes() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAberrationsConverter c =
            new FlybaseAberrationsConverter(writer, Model.getInstanceByName("genomic"));
        // first pass: synonyms (creates the Aberration items keyed by FBab)
        c.setCurrentFile(new File("fb_synonym_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fb_synonym_test.tsv")));
        // second pass: del/dup
        c.setCurrentFile(new File("aberration_del_dup_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/aberration_del_dup_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        // Find FBab0001001 — should have 2 deletedGenes (rows 1,2) and 0 duplicated
        Item ab1 = items.stream()
            .filter(i -> "Aberration".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBab0001001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBab0001001 not stored"));

        long deletedCount = ab1.getCollections().stream()
            .filter(c2 -> "deletedGenes".equals(c2.getName()))
            .mapToLong(c2 -> c2.getRefIds().size()).sum();
        assertEquals(2, deletedCount);

        long duplicatedCount = ab1.getCollections().stream()
            .filter(c2 -> "duplicatedGenes".equals(c2.getName()))
            .mapToLong(c2 -> c2.getRefIds().size()).sum();
        assertEquals(0, duplicatedCount);

        // The 2 "not deleted"/"not duplicated" rows should NOT create their
        // gene items (geneD, geneE). Total Gene items = 3 (geneA, B, C).
        long genes = items.stream()
            .filter(i -> "Gene".equals(i.getClassName())).count();
        assertEquals(3, genes);
    }

    /**
     * Bugfix 2026-06-02 — file iteration order is non-deterministic in
     * BioFileConverter. Run delDup BEFORE synonyms and confirm the resulting
     * stored items still have populated deletedGenes / duplicatedGenes
     * collections. (Symptom that motivated the fix: AllianceMineDev integrate
     * produced 23,870 Aberrations + 642 Balancers but 0 rows in
     * aberrationdeletedgenes/duplicatedgenes/balancercomposedofaberrations.)
     */
    public void testDelDupBeforeSynonymsStillWiresCollections() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAberrationsConverter c =
            new FlybaseAberrationsConverter(writer, Model.getInstanceByName("genomic"));
        // INTENTIONALLY OUT OF ORDER: del/dup before synonyms
        c.setCurrentFile(new File("aberration_del_dup_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/aberration_del_dup_test.tsv")));
        c.setCurrentFile(new File("fb_synonym_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fb_synonym_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        Item ab1 = items.stream()
            .filter(i -> "Aberration".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBab0001001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBab0001001 not stored"));

        long deletedCount = ab1.getCollections().stream()
            .filter(c2 -> "deletedGenes".equals(c2.getName()))
            .mapToLong(c2 -> c2.getRefIds().size()).sum();
        assertEquals("out-of-order: FBab0001001 should still have 2 deletedGenes",
                     2, deletedCount);

        // Also confirm the enriched FBab carries symbol + aberrationType from synonyms
        boolean hasSymbol = ab1.getAttributes().stream().anyMatch(a ->
            "symbol".equals(a.getName()) && "Df(2R)test1".equals(a.getValue()));
        boolean hasType = ab1.getAttributes().stream().anyMatch(a ->
            "aberrationType".equals(a.getName()) && "deletion".equals(a.getValue()));
        assertTrue("synonyms must enrich (not replace) the stub Aberration: symbol", hasSymbol);
        assertTrue("synonyms must enrich (not replace) the stub Aberration: type", hasType);
    }

    /**
     * Bugfix 2026-06-02 — companion .tsv.gz file must not be re-processed.
     * If both fb_synonym_test.tsv and fb_synonym_test.tsv.gz are passed in,
     * the .gz must be skipped to avoid garbage-parsing or duplicate stubs.
     */
    public void testGzCompanionIsSkipped() throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAberrationsConverter c =
            new FlybaseAberrationsConverter(writer, Model.getInstanceByName("genomic"));
        // .tsv first (normal processing), then a fictional .tsv.gz with the same name —
        // should be skipped without exception even though we pass a non-gzip Reader.
        c.setCurrentFile(new File("fb_synonym_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fb_synonym_test.tsv")));
        c.setCurrentFile(new File("fb_synonym_test.tsv.gz"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fb_synonym_test.tsv")));  // reader is real but file is .gz
        c.close();

        Collection<Item> items = writer.getItems();
        long aberrations = items.stream()
            .filter(i -> "Aberration".equals(i.getClassName())).count();
        // 5 from .tsv; .tsv.gz skipped → still 5
        assertEquals(5, aberrations);
    }

    /**
     * Task 10 — curated fbba_to_fbab table wires composedOfAberrations;
     * balancers with empty composition still exist as Balancer items.
     */
    public void testCuratedBalancersComposeAberrationsAndTolerateEmpty()
            throws Exception {
        MockItemWriter writer =
            new MockItemWriter(new LinkedHashMap<String, Item>());
        FlybaseAberrationsConverter c =
            new FlybaseAberrationsConverter(writer, Model.getInstanceByName("genomic"));
        c.setCurrentFile(new File("fb_synonym_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fb_synonym_test.tsv")));
        c.setCurrentFile(new File("fbba_to_fbab_test.tsv"));
        c.process(new InputStreamReader(
            getClass().getResourceAsStream("/fbba_to_fbab_test.tsv")));
        c.close();

        Collection<Item> items = writer.getItems();
        // FBba0001001 has 2 fbab_ids -> expect 2 in composedOfAberrations
        Item bal1 = items.stream()
            .filter(i -> "Balancer".equals(i.getClassName()))
            .filter(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBba0001001".equals(a.getValue())))
            .findFirst()
            .orElseThrow(() -> new AssertionError("FBba0001001 not stored"));
        long composedCount = bal1.getCollections().stream()
            .filter(col -> "composedOfAberrations".equals(col.getName()))
            .mapToLong(col -> col.getRefIds().size()).sum();
        assertEquals(2, composedCount);

        // FBba0001002 has empty fbab_ids — Balancer must still exist
        boolean bal2Exists = items.stream()
            .filter(i -> "Balancer".equals(i.getClassName()))
            .anyMatch(i -> i.getAttributes().stream().anyMatch(a ->
                "primaryIdentifier".equals(a.getName())
                && "FBba0001002".equals(a.getValue())));
        assertTrue("FBba0001002 (empty composition) should still be stored",
                   bal2Exists);
    }
}
