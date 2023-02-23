/*
 * The MIT License
 *
 * Copyright 2023 Rami Manaf Abdullah.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.jme.forcefield.mmff;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * inner class for handling MMFF94 parameter files
 * @author Rami Manaf Abdullah
 */
class MMFF94Parameters {

    final Hashtable<String, Integer> typeDefinition = new Hashtable<>();
    final Hashtable<Integer, Integer[]> atomTypeEquivalenceParameters = new Hashtable<>();
    //first static number is the atom which connected by other atoms in the nested map
    //and afstatic fected by the float bond charge when connected to them
    final List<Float[]> bondChargeIncrements = new ArrayList<>();
    final List<StretchParameters> bondStretchParameters = new ArrayList<>();
    final List<StretchEmpiricalParameters> bondStretchEmpiricalParameterses = new ArrayList<>();
    List<VdwParameters> vdwParameters;
    final List<Float[]> angleBendingParameters = new ArrayList<>();
    final List<Float[]> stretchBendParameters = new ArrayList<>();
    final List<Float[]> outOfPlaneParameters = new ArrayList<>();
    final List<Float[]> outOfPlaneParametersStatic = new ArrayList<>();
    final List<Float[]> torsionParameters = new ArrayList<>();
    final List<Float[]> torsionParametersStatic = new ArrayList<>();
    final List<GeometricParameters> geometricProperties = new ArrayList<>();
    final List<Float[]> partialChargeIncrement = new ArrayList<>();
    final List<Float[]> stretchBendEmpiricalParamters = new ArrayList<>();
    final float[][] angleBendingEmpiricalParamters = new float[][]{
        {1, 1.395f, -1},
        {5, -1, 0.704f},
        {6, 2.494f, 1.016f},
        {7, 2.711f, 1.113f},
        {8, 3.045f, 1.337f},
        {9, 2.847f, -1},
        {14, 2.350f, 0.811f},
        {15, 2.350f, 1.068f},
        {16, 2.980f, 1.249f},
        {17, 2.909f, 1.078f},
        {33, -1, 0.825f},
        {35, 3.017f, -1},
        {53, 3.086f, -1}};

    final float[][] torsionEmpiricalParamters = new float[][]{
        {6, 2, 2.12f},
        {7, 2, 1.5f},
        {8, 2, .2f},
        {14, 1.25f, 1.22f},
        {15, 1.25f, 2.4f},
        {16, 1.25f, .49f}
    };

    final Float[][] covalentRadiiPaulingElectronegativities = new Float[][]{
        //    {atomicNum, covRad, pauEle}, 
        {1f, 0.33f, 2.20f},
        {3f, 1.34f, 0.97f},
        {6f, 0.77f, 2.50f},
        {7f, 0.73f, 3.07f},
        {8f, 0.72f, 3.50f},
        {9f, 0.74f, 4.12f},
        {11f, 1.54f, 1.01f},
        {12f, 1.30f, 1.23f},
        {14f, 1.15f, 1.74f},
        {15f, 1.09f, 2.06f},
        {16f, 1.03f, 2.44f},
        {17f, 1.01f, 2.83f},
        {19f, 1.96f, 0.91f},
        {20f, 1.74f, 1.04f},
        {29f, 1.38f, 1.75f},
        {30f, 1.31f, 1.66f},
        {35f, 1.15f, 2.74f},
        {53f, 1.33f, 2.21f}};

    final Float[][] herschbachLaurie = new Float[][]{
        //    Parameters for Badger's Rule
        //    i j a_ij d_ij dp_ij
        {0f, 0f, 1.26f, 0.025f, 0.025f},
        {0f, 1f, 1.66f, 0.30f, 0.36f},
        {0f, 2f, 1.84f, 0.38f, 0.58f},
        {0f, 3f, 1.98f, 0.49f, 0.65f},
        {0f, 4f, 2.03f, 0.51f, 0.80f},
        {0f, 5f, 2.03f, 0.25f, 0.81f},
        {0f, 30f, 1.85f, 0.15f, 0.53f},
        {0f, 40f, 1.84f, 0.61f, 0.61f},
        {0f, 50f, 1.78f, 0.97f, 0.62f},
        {1f, 1f, 1.91f, 0.68f, 0.68f},
        {1f, 2f, 2.28f, 0.74f, 0.92f},
        {1f, 3f, 2.35f, 0.85f, 1.02f},
        {1f, 4f, 2.33f, 0.68f, 1.12f},
        {1f, 5f, 2.50f, 0.97f, 1.22f},
        {1f, 30f, 2.08f, 1.14f, 0.97f},
        {1f, 40f, 2.34f, 1.17f, 1.08f},
        {2f, 2f, 2.41f, 1.18f, 1.18f},
        {2f, 3f, 2.52f, 1.02f, 1.28f},
        {2f, 4f, 2.61f, 1.28f, 1.40f},
        {2f, 5f, 2.60f, 0.84f, 1.24f},
        {3f, 3f, 2.58f, 1.41f, 1.35f},
        {3f, 4f, 2.66f, 0.86f, 1.48f},
        {3f, 5f, 2.75f, 1.14f, 1.55f},
        {4f, 4f, 2.85f, 1.62f, 1.62f},
        {4f, 5f, 2.76f, 1.25f, 1.51f}};

    private static MMFF94Parameters INSTANCE;

    private MMFF94Parameters() throws IOException {
        loadParameters();
    }

    static MMFF94Parameters getInstance() throws IOException {
        if (INSTANCE == null) {
            INSTANCE = new MMFF94Parameters();
        }
        return INSTANCE;
    }

    private void loadParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("/org/openscience/cdk/forcefield/mmff/mmff-symb-mapping.tsv"), StandardCharsets.UTF_8))) {
            String line;
            reader.readLine();
            while ((line = reader.readLine()) != null) {
                String[] parameters = line.split("\t");
                typeDefinition.put(parameters[0], Integer.parseInt(parameters[1]));
            };
        }
        parseVdwParameters();
        parseElectrostaticParameters();
        parseBondStretchingParameters();
        parseBondStretchingEmpiricalParameters();
        parseAngleBendingParameters();
        parseStretchBendParameters();
        parseOutOfPlaneParameters(true);
        parseOutOfPlaneParameters(false);
        parseTorsionParameters(true);
        parseTorsionParameters(false);
        parseAtomTypeEquivalenceParameters();
        parseGeometricProperties();
        parsePartialBondChargeIncrement();
        parseStretchBendParameters();
        parseEmpiricalStretchBendParameters();
        parseMMFFFORMCHG(getClass().getResourceAsStream("/org/openscience/cdk/forcefield/mmff/MMFFFORMCHG.PAR"), fCharges);
    }
    Map<String, BigDecimal> fCharges = new HashMap<>();

    private static void parseMMFFFORMCHG(InputStream in, Map<String, BigDecimal> fcharges) throws IOException {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
            String line = br.readLine(); // header
            while ((line = br.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                final String[] cols = line.split("\\s+");
                fcharges.put(cols[0], new BigDecimal(cols[1]));
            }
        }
    }

    private void parseVdwParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-vdw.par"), StandardCharsets.UTF_8))) {
            //TODO replace streams with loops for performance
            vdwParameters = reader.lines().filter((t) -> !t.startsWith("*")).map((t) -> {
                String[] parameter = t.split("\t");
                return new VdwParameters(Integer.parseInt(parameter[0]), Float.parseFloat(parameter[1]), Float.parseFloat(parameter[2]), Float.parseFloat(parameter[3]), Float.parseFloat(parameter[4]), parameter[5]);
            }).collect(Collectors.toList());
        }
    }

    private void parseElectrostaticParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-bond-charge-increment.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                float type1 = Float.parseFloat(parameter[1]);
                float type2 = Float.parseFloat(parameter[2]);
                float bondCharge = Float.parseFloat(parameter[3]);
                bondChargeIncrements.add(new Float[]{Float.parseFloat(parameter[0]), type1, type2, -bondCharge});
            };
        }
    }

    private void parseBondStretchingParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-bond-stretch.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                int bondType = Integer.parseInt(parameter[0]);
                int type1 = Integer.parseInt(parameter[1]);
                int type2 = Integer.parseInt(parameter[2]);
                float kb = Float.parseFloat(parameter[3]);
                float r0 = Float.parseFloat(parameter[4]);
                bondStretchParameters.add(new StretchParameters(bondType, type1, type2, kb, r0));
            }
        }
    }

    private void parseBondStretchingEmpiricalParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-bond-stretch-empirical.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                int type1 = Integer.parseInt(parameter[0]);
                int type2 = Integer.parseInt(parameter[1]);
                float r0 = Float.parseFloat(parameter[2]);
                float kb = Float.parseFloat(parameter[3]);
                bondStretchEmpiricalParameterses.add(new StretchEmpiricalParameters(type1, type2, kb, r0));
            }
        }
    }

    private void parseAngleBendingParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-angle-bending.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                angleBendingParameters.add(new Float[]{Float.parseFloat(parameter[0]), Float.parseFloat(parameter[1]), Float.parseFloat(parameter[2]), Float.parseFloat(parameter[3]), Float.parseFloat(parameter[4]), Float.parseFloat(parameter[5])});
            }
        }
    }

    private void parseStretchBendParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-stretch-bend.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                stretchBendParameters.add(new Float[]{Float.parseFloat(parameter[0]), Float.parseFloat(parameter[1]), Float.parseFloat(parameter[2]), Float.parseFloat(parameter[3]), Float.parseFloat(parameter[4]), Float.parseFloat(parameter[5])});
            }
        }
    }

    private void parseOutOfPlaneParameters(boolean mmff94s) throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(mmff94s ? "MMFF94s-out-of-plane.par" : "MMFF94-out-of-plane.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                List<Float[]> oopParameters = mmff94s ? outOfPlaneParametersStatic : outOfPlaneParameters;
                oopParameters.add(new Float[]{Float.parseFloat(parameter[0]), Float.parseFloat(parameter[1]), Float.parseFloat(parameter[2]), Float.parseFloat(parameter[3]), Float.parseFloat(parameter[4])});
            }
        }
    }

    private void parseTorsionParameters(boolean mmff94s) throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(mmff94s ? "MMFF94s-torsion.par" : "MMFF94-torsion.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                List<Float[]> torsionParameters = mmff94s ? this.torsionParametersStatic : this.torsionParameters;
                torsionParameters.add(new Float[]{Float.parseFloat(parameter[0]), Float.parseFloat(parameter[1]), Float.parseFloat(parameter[2]), Float.parseFloat(parameter[3]), Float.parseFloat(parameter[4]), Float.parseFloat(parameter[5]), Float.parseFloat(parameter[6]), Float.parseFloat(parameter[7])});
            }
        }
    }

    private void parseAtomTypeEquivalenceParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-atom-type-equivalence.par"), StandardCharsets.UTF_8))) {
            String line;
            int n = 0;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                atomTypeEquivalenceParameters.put(Integer.valueOf(parameter[1]), new Integer[]{Integer.parseInt(parameter[2]), Integer.parseInt(parameter[3]), Integer.parseInt(parameter[4]), Integer.parseInt(parameter[5])});
                n++;
            }
        }
    }

    private void parseGeometricProperties() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-geometric-properties.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                geometricProperties.add(new GeometricParameters(Integer.parseInt(parameter[0]), Integer.parseInt(parameter[1]), Integer.parseInt(parameter[2]), Integer.parseInt(parameter[3]), parseBoolean(parameter[4]), Integer.parseInt(parameter[5]), parseBoolean(parameter[6]), parseBoolean(parameter[7]), parseBoolean(parameter[8])));
            }
        }
    }

    private void parsePartialBondChargeIncrement() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-partial-charge-increment.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                partialChargeIncrement.add(new Float[]{Float.valueOf(parameter[0]), Float.valueOf(parameter[1]), Float.valueOf(parameter[2])});
            }
        }
    }

    private void parseEmpiricalStretchBendParameters() throws IOException {
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream("MMFF94-stretch-bend-empirical.par"), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                String[] parameter = line.split("\t");
                stretchBendEmpiricalParamters.add(new Float[]{Float.parseFloat(parameter[0]), Float.parseFloat(parameter[1]), Float.parseFloat(parameter[2]), Float.parseFloat(parameter[3]), Float.parseFloat(parameter[4])});
            }
        }
    }

    private static boolean parseBoolean(String number) {
        return "1".equals(number);
    }

    static class VdwParameters {

        static final byte NONE = 0;
        static final byte DONER = 1;
        static final byte ACCEPTOR = 2;

        final int type;
        final float alpha, N, A, G;
        final byte DA;

        public VdwParameters(int type, float alpha, float N, float A, float G, String DA) {
            this.type = type;
            this.alpha = alpha;
            this.N = N;
            this.A = A;
            this.G = G;
            this.DA = DA.equals("D") ? DONER : (DA.equals("A") ? ACCEPTOR : NONE);
        }
    }

    static class GeometricParameters {

        /**
         * MMFF atom type
         */
        final int atomType;
        /**
         * atomic number
         */
        final int atomNumber;
        /**
         * the mandatory number of bonded neighbors
         */
        final int crd;
        /**
         * total number of bonds made to that atom type
         */
        final int valence;
        /**
         * pi lone pair capable of participating in resonance interactions with,
         * say, an adjacent multiple bond
         */
        final boolean piLonePair;
        /**
         * specifies cases in which double (2) or triple (3) bonds are expected
         * to be made to an atom having the listed atom type
         */
        final int mltb;
        final boolean aromatic;
        final boolean ideallyLinear;
        /**
         * designates atom types that can form either a multiple bond or a
         * delocalized single bond with an sp2- or sp-hybridized neighbor. An
         * example of such a single bond is the bond connecting the two central
         * carbons of butadiene
         */
        final boolean singleBondBetweenMultipleBonds;

        public GeometricParameters(int atomType, int atomNumber, int crd, int valence, boolean piLonePair, int mltb, boolean aromatic, boolean ideallyLinear, boolean singleBondBetweenMultipleBonds) {
            this.atomType = atomType;
            this.atomNumber = atomNumber;
            this.crd = crd;
            this.valence = valence;
            this.piLonePair = piLonePair;
            this.mltb = mltb;
            this.aromatic = aromatic;
            this.ideallyLinear = ideallyLinear;
            this.singleBondBetweenMultipleBonds = singleBondBetweenMultipleBonds;
        }
    }

    static class StretchParameters {

        final int bondType, i, j;
        final float kb, r0;

        public StretchParameters(int bondType, int i, int j, float kb, float r0) {
            this.bondType = bondType;
            this.i = i;
            this.j = j;
            this.kb = kb;
            this.r0 = r0;
        }
    }

    static class StretchEmpiricalParameters {

        final int i, j;
        final float kb, r0;

        public StretchEmpiricalParameters(int i, int j, float kb, float r0) {
            this.i = i;
            this.j = j;
            this.kb = kb;
            this.r0 = r0;
        }
    }

}
