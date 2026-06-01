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
}
