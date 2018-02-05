/* Copyright (C) 2009 Rajarshi Guha
 *               2009,2011 Egon Willighagen <egonw@users.sf.net>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 * Contact: Rajarshi Guha <rajarshi@users.sourceforge.net>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All I ask is that proper credit is given for my work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.fingerprint;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.fingerprint.ICountFingerprint;

import static org.hamcrest.CoreMatchers.is;

/**
 * @cdk.module test-fingerprint
 */
public class AtomPairs2DFingerprintTest extends AbstractFingerprinterTest {

    SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());

    @Test
    public void testFingerprint() throws Exception {
    	/*
    	 * We are going to test hexane. Hexane is a good test b/c it has 10 carbons.
    	 * Since the max distance for this fingerprint is 10, the final C-C fingerprint slot
    	 * at distance 10 should return false, while all the other C-C fingerprint slots
    	 * should return true.
    	 */
    	IFingerprinter printer = new AtomPairs2DFingerprinter();
        IAtomContainer mol1 = parser.parseSmiles("cccccccccc");
        BitSetFingerprint bsfp = (BitSetFingerprint) printer.getBitFingerprint(mol1);
        Assert.assertEquals(9, bsfp.cardinality());
        Assert.assertEquals(true, bsfp.get(0));		//Distance 1
        Assert.assertEquals(true, bsfp.get(78));	//Distance 2
        Assert.assertEquals(true, bsfp.get(156));	//Distance 3
        Assert.assertEquals(true, bsfp.get(234));	//Distance 4
        Assert.assertEquals(true, bsfp.get(312));	//Distance 5
        Assert.assertEquals(true, bsfp.get(390));	//Distance 6
        Assert.assertEquals(true, bsfp.get(468));	//Distance 7
        Assert.assertEquals(true, bsfp.get(546));	//Distance 8
        Assert.assertEquals(true, bsfp.get(624));	//Distance 9
        Assert.assertEquals(false, bsfp.get(702)); 	//Distance 10
    }
    
    @Test
    public void testGetCountFingerprint() throws Exception {
    	IFingerprinter printer = new AtomPairs2DFingerprinter();
        IAtomContainer mol1 = parser.parseSmiles("cccccccccc");
        ICountFingerprint icfp = printer.getCountFingerprint(mol1);
        Assert.assertEquals(780, icfp.numOfPopulatedbins());
        Assert.assertEquals(780, icfp.size());
        
    }
    
    @Test
    public void testGetRawFingerprint() throws Exception {
    	IFingerprinter printer = new AtomPairs2DFingerprinter();
    }
}
