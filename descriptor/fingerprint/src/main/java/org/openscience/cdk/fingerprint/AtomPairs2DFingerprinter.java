/*
 * Note: this code was ported from Yap Chun Wei's PaDEL source code, which can be found here:
 * http://www.yapcwsoft.com/dd/padeldescriptor/
 */



package org.openscience.cdk.fingerprint;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Map;
import java.util.TreeMap;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.matrix.TopologicalMatrix;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;



public class AtomPairs2DFingerprinter extends AbstractFingerprinter implements IFingerprinter {
	
	private static final int maxDistance = 10;
    public String[] names;
    private static final String[] atypes = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "X"};
    private ArrayList<Integer[]> atypesInt;
    
    private int[][] distance;
    private int[] counts;
    private Map<Integer, Integer> map = new TreeMap<Integer, Integer>();

    private class Pair {
        public Integer[] first;
        public Integer[] second;

        public Pair() {
        
        }

        public Pair(Integer[] first, Integer[] second) {
            this.first = first;
            this.second = second;
        }
    }
    
    private Pair[] atomPairs;

    public AtomPairs2DFingerprinter() {
        atypesInt = new ArrayList<Integer[]>();
        for (int i=0; i<atypes.length-1; ++i){
            atypesInt.add(new Integer[1]);
        }
        
        atypesInt.add(new Integer[4]);
        atypesInt.get(0)[0] = 6; // C
        atypesInt.get(1)[0] = 7; // N
        atypesInt.get(2)[0] = 8; // O
        atypesInt.get(3)[0] = 16; // S
        atypesInt.get(4)[0] = 15; // P
        atypesInt.get(5)[0] = 9; // F
        atypesInt.get(6)[0] = 17; // Cl
        atypesInt.get(7)[0] = 35; // Br
        atypesInt.get(8)[0] = 53; // I
        atypesInt.get(9)[0] = 5; // B
        atypesInt.get(10)[0] = 14; // Si
        atypesInt.get(11)[0] = 9; // X
        atypesInt.get(11)[1] = 17; // X
        atypesInt.get(11)[2] = 35; // X
        atypesInt.get(11)[3] = 53; // X

        atomPairs = new Pair[atypes.length*atypes.length/2 + atypes.length/2];
        int index = 0;
        for (int i=0; i<atypes.length; ++i) {
            for (int j=i; j<atypes.length; ++j){
                atomPairs[index++] = new Pair(atypesInt.get(i), atypesInt.get(j));
            }
        }

        index = 0;
        names = new String[maxDistance*atomPairs.length];
        for (int d=1; d<=maxDistance; ++d) {
            for (int i=0; i<atypes.length; ++i) {
            	for (int j=i; j<atypes.length; ++j) {
                    names[index++] = "AP2D" + d + "_" + atypes[i] + "_" + atypes[j];
                }
            }
        }
    }

    @Override
    public int getSize() {
        return maxDistance*atomPairs.length;
    }

    private boolean isInArray(Integer[] array, int val) {
        for (int i=0; i<array.length; ++i){
            if (array[i]==val) return true;
        }
        return false;
    }

	@Override
	public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
		IAtomContainer local = AtomContainerManipulator.removeHydrogens(container);
        int natom = local.getAtomCount();
        int[][] distance = TopologicalMatrix.getMatrix(local);
        BitSet fp = new BitSet(maxDistance*atomPairs.length);
        for (int i=0; i<natom; ++i){
            int a1 = local.getAtom(i).getAtomicNumber();
            for (int j=i+1; j<natom; ++j){
                int a2 = local.getAtom(j).getAtomicNumber();
                if (distance[i][j]<=maxDistance) {
                    for (int a=0; a<atomPairs.length; ++a){
                        if ((isInArray(atomPairs[a].first, a1) && isInArray(atomPairs[a].second, a2)) ||
                            (isInArray(atomPairs[a].first, a2) && isInArray(atomPairs[a].second, a1))) {
                            fp.set((distance[i][j]-1)*atomPairs.length + a, true);
                        }
                    }
                }
            }
        }
        return new BitSetFingerprint(fp);
	}

	/**
     * Invalid: it is not appropriate to convert the integer hash codes into strings.
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer mol) throws CDKException {
        throw new UnsupportedOperationException();
    }

	@Override
	public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {
		IAtomContainer local = AtomContainerManipulator.removeHydrogens(container);
        int natom = local.getAtomCount();
        distance = TopologicalMatrix.getMatrix(local);
        counts = new int[maxDistance * atomPairs.length];
        BitSet fp = new BitSet(maxDistance*atomPairs.length);
        for (int i=0; i<natom; ++i){
            int a1 = local.getAtom(i).getAtomicNumber();
            for (int j=i+1; j<natom; ++j){
                int a2 = local.getAtom(j).getAtomicNumber();
                if (distance[i][j]<=maxDistance) {
                    for (int a=0; a<atomPairs.length; ++a){
                        if ((isInArray(atomPairs[a].first, a1) && isInArray(atomPairs[a].second, a2)) ||
                            (isInArray(atomPairs[a].first, a2) && isInArray(atomPairs[a].second, a1))) {
                            counts[(distance[i][j]-1)*atomPairs.length + a]++;
                            fp.set((distance[i][j]-1)*atomPairs.length + a, true);
                        }
                    }
                }
            }
        }
        final BitSet fp2 = fp;
        return new ICountFingerprint() {

            @Override
            public long size() {
                return counts.length;
            }

            @Override
            public int numOfPopulatedbins() {
                return counts.length;
            }

            @Override
            public int getCount(int index){
            	return counts[index];
            }

            @Override
            public int getHash(int index) {
                if(fp2.get(index)){
                	return index;
                }
                else{
                	return 0;
                }
            }

            @Override
            public void merge(ICountFingerprint fp) {}

            @Override
            public void setBehaveAsBitFingerprint(boolean behaveAsBitFingerprint) {}

            @Override
            public boolean hasHash(int hash) {
                return fp2.get(hash);
            }

            @Override
            public int getCountForHash(int hash) {
            	if(counts.length < hash){
            		return counts[hash];
            	}
            	else{
            		return 0;
            	}
            }
        };
	}


}
