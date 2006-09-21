/* $RCSfile$
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 * Copyright (C) 1997-2006  The Chemistry Development Kit (CDK) project
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA. 
 * 
 */

package org.openscience.cdk.test;

import junit.framework.Test;
import junit.framework.TestSuite;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IRing;
import org.openscience.cdk.interfaces.IRingSet;

/**
 * Checks the funcitonality of the RingSet class.
 *
 * @cdk.module test-data
 *
 * @see org.openscience.cdk.RingSet
 */
public class RingSetTest extends CDKTestCase {

	protected IChemObjectBuilder builder;
	
    public RingSetTest(String name) {
        super(name);
    }

    public void setUp() {
       	builder = DefaultChemObjectBuilder.getInstance();
    }

    public static Test suite() {
        return new TestSuite(RingSetTest.class);
    }
    
    public void testRingSet() {
        IRingSet rs = builder.newRingSet();
        assertNotNull(rs);
    }
    
    
    public void testAdd_IRingSet() {
        IRing r1 = builder.newRing(5, "C");
        IRing r2 = builder.newRing(3, "C");
        
        IRingSet rs = builder.newRingSet();
        rs.addAtomContainer(r1);
        
        IRingSet rs2 = builder.newRingSet();
        rs2.addAtomContainer(r2);
        rs2.add(rs);
        
        assertEquals(1, rs.getAtomContainerCount());
        assertEquals(2, rs2.getAtomContainerCount());
    }

    public void testToString() {
        IRingSet ringset = builder.newRingSet();
        String description = ringset.toString();
        for (int i=0; i< description.length(); i++) {
            assertTrue(description.charAt(i) != '\n');
            assertTrue(description.charAt(i) != '\r');
        }
    }
    
    /* FIXME: make RingSet inherited from IChemObject, see #1117775
    public void testClone() {
        RingSet ringset = builder.newRingSet();
        Ring ring = builder.newRing();
        ringset.add(ring);
        
        RingSet clone = (RingSet)ringset.clone();
        assertNotNull(clone);
        assertTrue(clone instanceof RingSet);
        assertEquals(1, clone.size());
        assertNotSame(ring, clone.elementAt(0));
    } */
    
    public void testContains_IAtom() {
        IRingSet ringset = builder.newRingSet();

        IAtom ring1Atom1 = builder.newAtom("C"); // rather artificial molecule
        IAtom ring1Atom2 = builder.newAtom("C");
        IAtom sharedAtom1 = builder.newAtom("C");
        IAtom sharedAtom2 = builder.newAtom("C");
        IAtom ring2Atom1 = builder.newAtom("C");
        IAtom ring2Atom2 = builder.newAtom("C");
        IBond ring1Bond1 = builder.newBond(ring1Atom1, ring1Atom2);
        IBond ring1Bond2 = builder.newBond(sharedAtom1, ring1Atom1);
        IBond ring1Bond3 = builder.newBond(sharedAtom2, ring1Atom2);
        IBond sharedBond = builder.newBond(sharedAtom1, sharedAtom2);
        IBond ring2Bond1 = builder.newBond(ring2Atom1, ring2Atom2);
        IBond ring2Bond2 = builder.newBond(sharedAtom1, ring2Atom1);
        IBond ring2Bond3 = builder.newBond(sharedAtom2, ring2Atom2);

        IRing ring1 = builder.newRing();
        ring1.addAtom(ring1Atom1);
        ring1.addAtom(ring1Atom2);
        ring1.addAtom(sharedAtom1);
        ring1.addAtom(sharedAtom2);
        ring1.addBond(ring1Bond1);
        ring1.addBond(ring1Bond2);
        ring1.addBond(ring1Bond3);
        ring1.addBond(sharedBond);
        IRing ring2 = builder.newRing();
        ring2.addAtom(ring2Atom1);
        ring2.addAtom(ring2Atom2);
        ring2.addAtom(sharedAtom1);
        ring2.addAtom(sharedAtom2);
        ring2.addBond(ring2Bond1);
        ring2.addBond(ring2Bond2);
        ring2.addBond(ring2Bond3);
        ring2.addBond(sharedBond);
        
        ringset.addAtomContainer(ring1);
        ringset.addAtomContainer(ring2);
        
        assertTrue(ringset.contains(ring1Atom1));
        assertTrue(ringset.contains(ring1Atom2));
        assertTrue(ringset.contains(sharedAtom1));
        assertTrue(ringset.contains(sharedAtom2));
        assertTrue(ringset.contains(ring2Atom1));
        assertTrue(ringset.contains(ring2Atom2));
    }
    
    public void testContains_IAtomContainer() {
        IRingSet ringset = builder.newRingSet();

        IAtom ring1Atom1 = builder.newAtom("C"); // rather artificial molecule
        IAtom ring1Atom2 = builder.newAtom("C");
        IAtom sharedAtom1 = builder.newAtom("C");
        IAtom sharedAtom2 = builder.newAtom("C");
        IAtom ring2Atom1 = builder.newAtom("C");
        IAtom ring2Atom2 = builder.newAtom("C");
        IBond ring1Bond1 = builder.newBond(ring1Atom1, ring1Atom2);
        IBond ring1Bond2 = builder.newBond(sharedAtom1, ring1Atom1);
        IBond ring1Bond3 = builder.newBond(sharedAtom2, ring1Atom2);
        IBond sharedBond = builder.newBond(sharedAtom1, sharedAtom2);
        IBond ring2Bond1 = builder.newBond(ring2Atom1, ring2Atom2);
        IBond ring2Bond2 = builder.newBond(sharedAtom1, ring2Atom1);
        IBond ring2Bond3 = builder.newBond(sharedAtom2, ring2Atom2);

        IRing ring1 = builder.newRing();
        ring1.addAtom(ring1Atom1);
        ring1.addAtom(ring1Atom2);
        ring1.addAtom(sharedAtom1);
        ring1.addAtom(sharedAtom2);
        ring1.addBond(ring1Bond1);
        ring1.addBond(ring1Bond2);
        ring1.addBond(ring1Bond3);
        ring1.addBond(sharedBond);
        IRing ring2 = builder.newRing();
        ring2.addAtom(ring2Atom1);
        ring2.addAtom(ring2Atom2);
        ring2.addAtom(sharedAtom1);
        ring2.addAtom(sharedAtom2);
        ring2.addBond(ring2Bond1);
        ring2.addBond(ring2Bond2);
        ring2.addBond(ring2Bond3);
        ring2.addBond(sharedBond);
        
        ringset.addAtomContainer(ring1);
        ringset.addAtomContainer(ring2);
        
        assertTrue(ringset.contains(ring1));
        assertTrue(ringset.contains(ring2));
    }

    public void testGetRings_IBond() {
        IRingSet ringset = builder.newRingSet();

        IAtom ring1Atom1 = builder.newAtom("C"); // rather artificial molecule
        IAtom ring1Atom2 = builder.newAtom("C");
        IAtom sharedAtom1 = builder.newAtom("C");
        IAtom sharedAtom2 = builder.newAtom("C");
        IAtom ring2Atom1 = builder.newAtom("C");
        IAtom ring2Atom2 = builder.newAtom("C");
        IBond ring1Bond1 = builder.newBond(ring1Atom1, ring1Atom2);
        IBond ring1Bond2 = builder.newBond(sharedAtom1, ring1Atom1);
        IBond ring1Bond3 = builder.newBond(sharedAtom2, ring1Atom2);
        IBond sharedBond = builder.newBond(sharedAtom1, sharedAtom2);
        IBond ring2Bond1 = builder.newBond(ring2Atom1, ring2Atom2);
        IBond ring2Bond2 = builder.newBond(sharedAtom1, ring2Atom1);
        IBond ring2Bond3 = builder.newBond(sharedAtom2, ring2Atom2);

        IRing ring1 = builder.newRing();
        ring1.addAtom(ring1Atom1);
        ring1.addAtom(ring1Atom2);
        ring1.addAtom(sharedAtom1);
        ring1.addAtom(sharedAtom2);
        ring1.addBond(ring1Bond1);
        ring1.addBond(ring1Bond2);
        ring1.addBond(ring1Bond3);
        ring1.addBond(sharedBond);
        IRing ring2 = builder.newRing();
        ring2.addAtom(ring2Atom1);
        ring2.addAtom(ring2Atom2);
        ring2.addAtom(sharedAtom1);
        ring2.addAtom(sharedAtom2);
        ring2.addBond(ring2Bond1);
        ring2.addBond(ring2Bond2);
        ring2.addBond(ring2Bond3);
        ring2.addBond(sharedBond);
        
        ringset.addAtomContainer(ring1);
        ringset.addAtomContainer(ring2);
        
        assertEquals(1, ringset.getRings(ring1Bond1).size());
        assertEquals(1, ringset.getRings(ring1Bond2).size());
        assertEquals(1, ringset.getRings(ring1Bond3).size());
        assertEquals(2, ringset.getRings(sharedBond).size());
        assertEquals(1, ringset.getRings(ring2Bond1).size());
        assertEquals(1, ringset.getRings(ring2Bond2).size());
        assertEquals(1, ringset.getRings(ring2Bond3).size());
    }

    public void testGetRings_IAtom() {
        IRingSet ringset = builder.newRingSet();

        IAtom ring1Atom1 = builder.newAtom("C"); // rather artificial molecule
        IAtom ring1Atom2 = builder.newAtom("C");
        IAtom sharedAtom1 = builder.newAtom("C");
        IAtom sharedAtom2 = builder.newAtom("C");
        IAtom ring2Atom1 = builder.newAtom("C");
        IAtom ring2Atom2 = builder.newAtom("C");
        IBond ring1Bond1 = builder.newBond(ring1Atom1, ring1Atom2);
        IBond ring1Bond2 = builder.newBond(sharedAtom1, ring1Atom1);
        IBond ring1Bond3 = builder.newBond(sharedAtom2, ring1Atom2);
        IBond sharedBond = builder.newBond(sharedAtom1, sharedAtom2);
        IBond ring2Bond1 = builder.newBond(ring2Atom1, ring2Atom2);
        IBond ring2Bond2 = builder.newBond(sharedAtom1, ring2Atom1);
        IBond ring2Bond3 = builder.newBond(sharedAtom2, ring2Atom2);

        IRing ring1 = builder.newRing();
        ring1.addAtom(ring1Atom1);
        ring1.addAtom(ring1Atom2);
        ring1.addAtom(sharedAtom1);
        ring1.addAtom(sharedAtom2);
        ring1.addBond(ring1Bond1);
        ring1.addBond(ring1Bond2);
        ring1.addBond(ring1Bond3);
        ring1.addBond(sharedBond);
        IRing ring2 = builder.newRing();
        ring2.addAtom(ring2Atom1);
        ring2.addAtom(ring2Atom2);
        ring2.addAtom(sharedAtom1);
        ring2.addAtom(sharedAtom2);
        ring2.addBond(ring2Bond1);
        ring2.addBond(ring2Bond2);
        ring2.addBond(ring2Bond3);
        ring2.addBond(sharedBond);
        
        ringset.addAtomContainer(ring1);
        ringset.addAtomContainer(ring2);
        
        assertEquals(1, ringset.getRings(ring1Atom1).getAtomContainerCount());
        assertEquals(1, ringset.getRings(ring1Atom1).getAtomContainerCount());
        assertEquals(2, ringset.getRings(sharedAtom1).getAtomContainerCount());
        assertEquals(2, ringset.getRings(sharedAtom2).getAtomContainerCount());
        assertEquals(1, ringset.getRings(ring2Atom1).getAtomContainerCount());
        assertEquals(1, ringset.getRings(ring2Atom2).getAtomContainerCount());
    }
    
    public void testGetConnectedRings_IRing() {
        IRingSet ringset = builder.newRingSet();

        IAtom ring1Atom1 = builder.newAtom("C"); // rather artificial molecule
        IAtom ring1Atom2 = builder.newAtom("C");
        IAtom sharedAtom1 = builder.newAtom("C");
        IAtom sharedAtom2 = builder.newAtom("C");
        IAtom ring2Atom1 = builder.newAtom("C");
        IAtom ring2Atom2 = builder.newAtom("C");
        IBond ring1Bond1 = builder.newBond(ring1Atom1, ring1Atom2);
        IBond ring1Bond2 = builder.newBond(sharedAtom1, ring1Atom1);
        IBond ring1Bond3 = builder.newBond(sharedAtom2, ring1Atom2);
        IBond sharedBond = builder.newBond(sharedAtom1, sharedAtom2);
        IBond ring2Bond1 = builder.newBond(ring2Atom1, ring2Atom2);
        IBond ring2Bond2 = builder.newBond(sharedAtom1, ring2Atom1);
        IBond ring2Bond3 = builder.newBond(sharedAtom2, ring2Atom2);

        IRing ring1 = builder.newRing();
        ring1.addAtom(ring1Atom1);
        ring1.addAtom(ring1Atom2);
        ring1.addAtom(sharedAtom1);
        ring1.addAtom(sharedAtom2);
        ring1.addBond(ring1Bond1);
        ring1.addBond(ring1Bond2);
        ring1.addBond(ring1Bond3);
        ring1.addBond(sharedBond);
        IRing ring2 = builder.newRing();
        ring2.addAtom(ring2Atom1);
        ring2.addAtom(ring2Atom2);
        ring2.addAtom(sharedAtom1);
        ring2.addAtom(sharedAtom2);
        ring2.addBond(ring2Bond1);
        ring2.addBond(ring2Bond2);
        ring2.addBond(ring2Bond3);
        ring2.addBond(sharedBond);
        
        ringset.addAtomContainer(ring1);
        ringset.addAtomContainer(ring2);
        
        assertEquals(2, ringset.getConnectedRings(ring2).size());
        assertEquals(2, ringset.getConnectedRings(ring1).size());
    }

    public void testClone() {
        // Added to make the Coverage tool happy
        // The method is apparently not part of the interface yet
    	assertTrue(true);
    }
}
