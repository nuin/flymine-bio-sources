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

import junit.framework.TestCase;

/**
 * Unit tests for FlybaseAberrationsConverter.
 * Expanded across tasks 7-10; this initial test exercises the static
 * symbol-prefix parser.
 */
public class FlybaseAberrationsConverterTest extends TestCase
{
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
}
