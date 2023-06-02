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

import org.jme.forcefield.ForceField;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.forcefield.mmff.Mmff;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.jme.forcefield.mmff.MMFF94Parameters.VdwParameters;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.RingSearch;

/**
 * Merck Molecular Force Field (MMFF94) implementation. The implementation is
 * based on MMFF seven papers, JMol, and RDKit implementations and include both
 * MMFF94 and MMFF94s. All the calculated energies are in kcal/mol unit.
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94 implements ForceField {

    static final String MMFF94_TYPE = "mmff.type";
    static final String MMFF94_RINGS = "mmff.rings";
    static final String MMFF94_RINGS_ISOLATED = "mmff.rings.isolated";
    static final String MMFF94_PARAMETER_STRETCH = "mmff.parameters.stretch";
    static final String MMFF94_PARAMETER_GEOMETRIC_PROPERTIES = "mmff.parameters.gemetricProperties";
    static final String MMFF94_FORMAL_CHARGE = "mmff.formalCharge";
    static final Logger LOGGER = Logger.getLogger(MMFF94.class.getName());
    private final Mmff mmff = new Mmff();
    private MMFF94Parameters mmffParameters;
    private boolean mmff94s;
    private boolean debug = false;

    /**
     * initialize MMFF94 and load its parameters
     *
     * @param mmff94s true for mmff94s and false for mmff94
     */
    public MMFF94(boolean mmff94s) {
        this.mmff94s = mmff94s;
        try {
            this.mmffParameters = MMFF94Parameters.getInstance();
        } catch (IOException ex) {
            LOGGER.log(Level.SEVERE, "Could not load MMFF94 parameter files", ex);
            throw new RuntimeException(ex);
        }
    }

    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    /**
     * returns if debug mode is on
     *
     * @return
     */
    public boolean isDebug() {
        return debug;
    }

    /**
     * preliminary assignation of MMFF94 parameters for the molecule. This is
     * required before calculating any energy is possible.
     *
     * @param atomContainer
     */
    public void assignParameters(IAtomContainer atomContainer) {
        if (!mmff.assignAtomTypes(atomContainer)) {
            LOGGER.warning("Not all atom types are assigned");
        }
        for (IAtom atom : atomContainer.atoms()) {
            int atomType = mmffParameters.typeDefinition.get(atom.getAtomTypeName());
            atom.setProperty(MMFF94_TYPE, atomType);
            atom.setFlag(CDKConstants.ISAROMATIC, false);
            //there is a 4 number gap in mmff atom type number (83-86)
            atom.setProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES, mmffParameters.geometricProperties.get(atomType - (atomType >= 87 ? 5 : 1)));
        }
        for (IBond bond : atomContainer.bonds()) {
            if (bond.getProperty("mmff.arom") != null && ((boolean) bond.getProperty("mmff.arom")) == true) {
                bond.setFlag(CDKConstants.ISAROMATIC, true);
                bond.getBegin().setFlag(CDKConstants.ISAROMATIC, true);
                bond.getEnd().setFlag(CDKConstants.ISAROMATIC, true);
            } else {
                bond.setFlag(CDKConstants.ISAROMATIC, false);
            }
        }
        RingSearch ringSearch = new RingSearch(atomContainer);
        try {
            atomContainer.setProperty(MMFF94_RINGS, new AllRingsFinder().findAllRings(atomContainer));
        } catch (CDKException ex) {
            throw new RuntimeException(ex);
        }
        atomContainer.setProperty(MMFF94_RINGS_ISOLATED, ringSearch.isolated());
        //after assigning atom and bond types and found rings we can now assign parameters

        calculateFormalCharge(atomContainer);
        calculatePartialCharge(atomContainer);
        if (debug) {
            for (IAtom atom : atomContainer.atoms()) {
                System.out.println(atom.getAtomTypeName() + ":\t" + atom.getCharge());
            }
        }
        for (IBond bond : atomContainer.bonds()) {
            IAtom iAtom, jAtom;
            if (bond.getBegin().<Integer>getProperty(MMFF94_TYPE) <= bond.getEnd().<Integer>getProperty(MMFF94_TYPE)) {
                iAtom = bond.getBegin();
                jAtom = bond.getEnd();
            } else {
                iAtom = bond.getEnd();
                jAtom = bond.getBegin();
            }
            int i = iAtom.getProperty(MMFF94_TYPE);
            int j = jAtom.getProperty(MMFF94_TYPE);
            int bondType = findBondType(iAtom, jAtom);
            boolean found = false;
            for (MMFF94Parameters.StretchParameters parameters : mmffParameters.bondStretchParameters) {
                if (parameters.i == i && parameters.j == j && parameters.bondType == bondType) {
                    bond.setProperty(MMFF94_PARAMETER_STRETCH, parameters);
                    found = true;
                    break;
                }
            }
            if (!found) {
                bond.setProperty(MMFF94_PARAMETER_STRETCH, calculateEmpiricalStretchParameters(iAtom, jAtom));
            }
        }
    }

    private MMFF94Parameters.StretchParameters calculateEmpiricalStretchParameters(IAtom iAtom, IAtom jAtom) {
        int atomicNum1, atomicNum2;
        if (iAtom.getAtomicNumber() <= jAtom.getAtomicNumber()) {
            atomicNum1 = iAtom.getAtomicNumber();
            atomicNum2 = jAtom.getAtomicNumber();
        } else {
            atomicNum1 = jAtom.getAtomicNumber();
            atomicNum2 = iAtom.getAtomicNumber();
        }
        MMFF94Parameters.StretchEmpiricalParameters mmffBndkParams = null;
        for (MMFF94Parameters.StretchEmpiricalParameters parameters : mmffParameters.bondStretchEmpiricalParameterses) {
            if (parameters.i == atomicNum1 && parameters.j == atomicNum2) {
                mmffBndkParams = parameters;
                break;
            }
        }
        Float[][] radiiElectronegativityParameters = new Float[2][3];
        int found = 0;
        for (Float[] parameters : mmffParameters.covalentRadiiPaulingElectronegativities) {
            if (parameters[0] == atomicNum1) {
                radiiElectronegativityParameters[0] = parameters;
                found++;
            }
            if (parameters[0] == atomicNum2) {
                radiiElectronegativityParameters[1] = parameters;
                found++;
            }
            if (found == 2) {
                break;
            }
        }
        if (found != 2) {
            throw new RuntimeException(String.format("could not find empirical paramters for %s-%s bond", iAtom.getAtomTypeName(), jAtom.getAtomTypeName()));
        }
        double c = (((atomicNum1 == 1) || (atomicNum2 == 1)) ? 0.050 : 0.085);
        double n = 1.4;
        double delta = 0.0;
        float[] r0_i = new float[2];

        for (int i = 0; i < 2; i++) {
            r0_i[i] = radiiElectronegativityParameters[i][1];
        }
        float r0 = (float) (r0_i[0] + r0_i[1] - c * Math.pow(Math.abs(radiiElectronegativityParameters[0][2] - radiiElectronegativityParameters[1][2]), n) - delta);
        float kb;
        if (mmffBndkParams != null) {
            kb = (float) (mmffBndkParams.kb * Math.pow(mmffBndkParams.r0 / r0, 6));
        } else {
            int row1 = getHerschbachLauriePeriodicTableRow(iAtom);
            int row2 = getHerschbachLauriePeriodicTableRow(jAtom);
            if (row1 > row2) {
                int temp = row2;
                row2 = row1;
                row1 = temp;
            }
            Float[] HerschbachLaurieParameters = null;
            for (Float[] parameters : mmffParameters.herschbachLaurie) {
                if (row1 == parameters[0] && row2 == parameters[1]) {
                    HerschbachLaurieParameters = parameters;
                    break;
                }
            }
            Objects.requireNonNull(HerschbachLaurieParameters, String.format("could not find Herschbach-Laurie parameters for %s-%s bond", iAtom.getAtomTypeName(), jAtom.getAtomTypeName()));
            kb = (float) Math.pow(10.0, -(r0 - HerschbachLaurieParameters[2]) / HerschbachLaurieParameters[3]);
        }

        return new MMFF94Parameters.StretchParameters(findBondType(iAtom, jAtom), iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), kb, r0);
    }
    
    private void checkParapetersAssigned(IAtom iAtom){
        if (iAtom.getProperty(MMFF94_TYPE) == null) {
            throw new RuntimeException("parameters need to be assigned first");
        }
    }

    /**
     * calculates the potential energy for the molecule in kcal/mol. This
     * includes all energy terms in MMFF94
     *
     * @param atomContainer
     * @return
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        if (!atomContainer.isEmpty()) {
            checkParapetersAssigned(atomContainer.getAtom(0));
        }
        ArrayList<String>[] arr = new ArrayList[7];
        if (debug) {
            for (int i = 0; i < 7; i++) {
                arr[i] = new ArrayList<>();
            }
        }
        double energy = 0;
        double totalStretch = 0;
        double totalVdw = 0;
        double totalElectrostatic = 0;
        double totalAngleBend = 0;
        double totalStretchBend = 0;
        double totalOOP = 0;
        double totalTorsion = 0;
        for (IBond bond : atomContainer.bonds()) {
            double x = calculateBondStretchingEnergy(bond);
            energy += x;
            totalStretch += x;
            if (debug) {
                arr[0].add(String.format("Stretch:\t%s-%s\t%f", bond.getBegin().getAtomTypeName(), bond.getEnd().getAtomTypeName(), x));
            }
        }
        for (IAtom atom : atomContainer.atoms()) {
            //non bonded interaction
            for (IAtom atom2 : atomContainer.atoms()) {
                if (atom != atom2 && atom.getIndex() < atom2.getIndex() && isSeparatedByNOrMoreBonds(atom, atom2, 3)) {
                    double x = calculateVdwEnergy(atom, atom2);
                    energy += x;
                    totalVdw += x;
                    if (debug) {
                        arr[5].add(String.format("VDW:\t%s-%s\t%f", atom.getAtomTypeName(), atom2.getAtomTypeName(), x));
                    }
                    double y = calculateElectrostaticInteractionEnergy(atom, atom2, 1);
                    energy += y;
                    totalElectrostatic += y;
                    if (debug) {
                        arr[6].add(String.format("Electro:\t%s-%s\t%f", atom.getAtomTypeName(), atom2.getAtomTypeName(), y));
                    }
                }
            }
            //bonded interaction
            if (atom.getBondCount() >= 2) {
                for (IAtom atom2 : atomContainer.getConnectedAtomsList(atom)) {
                    for (IAtom atom3 : atomContainer.getConnectedAtomsList(atom)) {
                        if (atom3 == atom2) {
                            continue;
                        }
                        int j = atom.getProperty(MMFF94_TYPE);
                        int i = atom2.getProperty(MMFF94_TYPE);
                        int k = atom3.getProperty(MMFF94_TYPE);
                        if (i < k || (i == k && atom2.getIndex() < atom3.getIndex())) {
                            double x = calculateAngleBendingEnergy(atom2, atom, atom3);
                            energy += x;
                            totalAngleBend += x;
                            if (debug) {
                                arr[1].add(String.format("Bending:\t%d-%d-%d\t%.3f", atom2.getProperty(MMFF94_TYPE), atom.getProperty(MMFF94_TYPE), atom3.getProperty(MMFF94_TYPE), x));
                            }
                            double y = calculateStretchBendEnergy(atom2, atom, atom3);
                            energy += y;
                            totalStretchBend += y;
                            if (debug) {
                                arr[2].add(String.format("Str-Bnd:\t%d-%d-%d\t%.3f", atom2.getProperty(MMFF94_TYPE), atom.getProperty(MMFF94_TYPE), atom3.getProperty(MMFF94_TYPE), y));
                            }
                        }
                        if (atom.getBondCount() >= 3) {
                            for (IAtom atom4 : atomContainer.getConnectedAtomsList(atom)) {
                                if (atom4 == atom2 || atom4 == atom3) {
                                    continue;
                                }
                                if (i < k || (k == i && atom2.getIndex() < atom3.getIndex())) {
                                    double z = calculateOutOfPlaneEnergy(atom2, atom, atom3, atom4);
                                    energy += z;
                                    totalOOP += z;
                                    if (debug) {
                                        arr[3].add(String.format("OOP:\t%d-%d-%d-%d\t%.3f", atom2.getProperty(MMFF94_TYPE), atom.getProperty(MMFF94_TYPE), atom3.getProperty(MMFF94_TYPE), atom4.getProperty(MMFF94_TYPE), z));
                                    }
                                }
                            }
                        }
                        //torsion
                        if (!atom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear && !atom3.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear && atom3.getBondCount() >= 2) {
                            for (IAtom atom4 : atomContainer.getConnectedAtomsList(atom3)) {
                                if (atom4 == atom || atom4 == atom2) {
                                    continue;
                                }
                                int l = atom4.getProperty(MMFF94_TYPE);
                                //prevent torsion energy from calculating twice by forcing torsion paramter file order
                                if (j < k || (j == k && i < l) || (j == k && i == l && atom.getIndex() < atom3.getIndex())) {
                                    double x = calculateTorsionEnergy(atom2, atom, atom3, atom4);
                                    energy += x;
                                    totalTorsion += x;
                                    if (debug) {
                                        arr[4].add(String.format("Tor:\t%d-%d-%d-%d\t%f", atom2.getProperty(MMFF94_TYPE), atom.getProperty(MMFF94_TYPE), atom3.getProperty(MMFF94_TYPE), atom4.getProperty(MMFF94_TYPE), x));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (debug) {
            arr[0].add("total stretch:\t" + totalStretch);
            arr[1].add("total angle bend:\t" + totalAngleBend);
            arr[2].add("total str-bend:\t" + totalStretchBend);
            arr[3].add("total OOP:\t" + totalOOP);
            arr[4].add("total torsion:\t" + totalTorsion);
            arr[5].add("total Vdw:\t" + totalVdw);
            arr[6].add("total Electro:\t" + totalElectrostatic);
            for (ArrayList<String> arrayList : arr) {
                arrayList.forEach(System.out::println);
            }
        }
        return energy;
    }

    /**
     * calculate the bond stretching energy in the provided bond
     *
     * @param bond
     * @return
     */
    public double calculateBondStretchingEnergy(IBond bond) {
        IAtom iAtom = bond.getBegin();
        checkParapetersAssigned(iAtom);
        IAtom jAtom = bond.getEnd();
        MMFF94Parameters.StretchParameters parameters = bond.getProperty(MMFF94_PARAMETER_STRETCH);
        double deltaR = iAtom.getPoint3d().distance(jAtom.getPoint3d()) - parameters.r0;
        return .5 * 143.9325 * parameters.kb * deltaR * deltaR * (1 - 2 * deltaR + (7d / 12d) * 4 * deltaR * deltaR);
    }

    /**
     * calculate Van Der Waals energy between two atoms. Generally
     * intermolecular interactions are calculated between non-bonded atoms or
     * atoms separated by three bonds (two atoms at least)
     *
     * @param iAtom
     * @param jAtom
     * @return
     */
    public double calculateVdwEnergy(IAtom iAtom, IAtom jAtom) {
        checkParapetersAssigned(iAtom);
        checkParapetersAssigned(jAtom);
        int i = iAtom.getProperty(MMFF94_TYPE);
        int j = jAtom.getProperty(MMFF94_TYPE);
        //there is a 4 number gap in mmff atom type number (83-86)
        VdwParameters iParameters = mmffParameters.vdwParameters.get(i - (i >= 87 ? 5 : 1));
        VdwParameters jParameters = mmffParameters.vdwParameters.get(j - (j >= 87 ? 5 : 1));
        double RStarIJ;
        if (iAtom.getAtomTypeName().equals(jAtom.getAtomTypeName())) {
            //like pairs
            RStarIJ = iParameters.A * Math.pow(iParameters.alpha, 0.25d);
        } else {
            //unlike pairs
            double RStarII = iParameters.A * Math.pow(iParameters.alpha, 0.25d);
            double RStarJJ = jParameters.A * Math.pow(jParameters.alpha, 0.25d);

            double B = (iParameters.DA == VdwParameters.DONER || jParameters.DA == VdwParameters.DONER) ? 0 : 0.2;
            double beta = 12;
            double gamma = (RStarII - RStarJJ) / (RStarII + RStarJJ);
            RStarIJ = 0.5 * (RStarII + RStarJJ) * (1 + B * (1 - Math.exp(-beta * Math.pow(gamma, 2))));
        }

        double epsilonIJ = 181.16 * iParameters.G * jParameters.G * iParameters.alpha * jParameters.alpha;
        epsilonIJ /= (Math.sqrt(iParameters.alpha / iParameters.N) + Math.sqrt(jParameters.alpha / jParameters.N));
        epsilonIJ /= Math.pow(RStarIJ, 6);

        if (iParameters.DA + jParameters.DA == 3) {
            //one is doner and the other is acceptor
            //we must calculate epsilon with the unscaled RStarIJ
            RStarIJ *= 0.8;
            epsilonIJ *= 0.5;
        }
        double distance = iAtom.getPoint3d().distance(jAtom.getPoint3d());
        double vdwEnergy = epsilonIJ * Math.pow((1.07 * RStarIJ) / (distance + 0.07 * RStarIJ), 7);
        vdwEnergy = vdwEnergy * (((1.12 * Math.pow(RStarIJ, 7)) / (Math.pow(distance, 7) + 0.12 * Math.pow(RStarIJ, 7))) - 2);
        return vdwEnergy;
    }

    /**
     * calculates electrostatic interaction energy between the two atoms.
     * Generally intermolecular interactions are calculated between non-bonded
     * atoms or atoms separated by three bonds (two atoms at least)
     *
     * @param iAtom
     * @param jAtom
     * @param dielectricConstant dielectric constant of the solvent. Default
     * value is 1
     * @return
     */
    public double calculateElectrostaticInteractionEnergy(IAtom iAtom, IAtom jAtom, double dielectricConstant) {
        checkParapetersAssigned(iAtom);
        checkParapetersAssigned(jAtom);
        double qI = iAtom.getCharge();
        double qJ = jAtom.getCharge();
        double electrostaticEnergy = 332.0716 * qI * qJ / (dielectricConstant * (iAtom.getPoint3d().distance(jAtom.getPoint3d()) + 0.05));
        if (iAtom.getContainer() == jAtom.getContainer() && isSeparatedByNBonds(iAtom, jAtom, 3)) {
            return electrostaticEnergy * 0.75;
        } else {
            return electrostaticEnergy;
        }
    }

    /**
     *
     * @param iAtom
     * @return
     */
    private void calculateFormalCharge(IAtomContainer container) {
        boolean[] conjugates = new boolean[container.getAtomCount()];
        for (IAtom iAtom : container.atoms()) {
            iAtom.setProperty(MMFF94_FORMAL_CHARGE, 0d);
            ATOM_TYPE_SWITCH:
            switch ((int) iAtom.getProperty(MMFF94_TYPE)) {
                case 32:
                case 72: {
                    int secondDegreeAmine = 0;
                    int terminalOxygenOrSulfer = 0;
                    for (IAtom neighbor : iAtom.getContainer().getConnectedAtomsList(iAtom)) {
                        for (IAtom neighbor2 : iAtom.getContainer().getConnectedAtomsList(neighbor)) {
                            if (neighbor2.getAtomicNumber() == 7 && neighbor2.getBondCount() == 2 && !neighbor2.isAromatic()) {
                                secondDegreeAmine++;
                            } else if ((neighbor2.getAtomicNumber() == 8 || neighbor2.getAtomicNumber() == 16) && neighbor2.getBondCount() == 1) {
                                terminalOxygenOrSulfer++;
                            }
                        }
                        if (neighbor.getAtomicNumber() == 16 && terminalOxygenOrSulfer == 2 && secondDegreeAmine == 1) {
                            //sulfonamide
                            secondDegreeAmine = 0;
                        }
                        if (neighbor.getAtomicNumber() == 6 && terminalOxygenOrSulfer != 0) {
                            if (terminalOxygenOrSulfer == 1) {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, -1d);
                            } else {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, -(terminalOxygenOrSulfer - 1) / (double) terminalOxygenOrSulfer);
                            }
                            break ATOM_TYPE_SWITCH;
                        }
                        int neighborType = neighbor.getProperty(MMFF94_TYPE);
                        if (neighborType == 45 && terminalOxygenOrSulfer == 3) {
                            iAtom.setProperty(MMFF94_FORMAL_CHARGE, -1d / 3d);
                            break ATOM_TYPE_SWITCH;
                        }
                        if (neighborType == 25 && terminalOxygenOrSulfer != 0) {
                            if (terminalOxygenOrSulfer == 1) {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, 0d);
                            } else {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, -(terminalOxygenOrSulfer - 1) / (double) terminalOxygenOrSulfer);
                            }
                            break ATOM_TYPE_SWITCH;
                        }
                        if (neighborType == 18 && terminalOxygenOrSulfer != 0) {
                            if (secondDegreeAmine + terminalOxygenOrSulfer == 2) {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, 0d);
                            } else {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, -(secondDegreeAmine + terminalOxygenOrSulfer - 2) / (double) terminalOxygenOrSulfer);
                            }
                            break ATOM_TYPE_SWITCH;
                        }
                        if (neighborType == 73 && terminalOxygenOrSulfer != 0) {
                            if (terminalOxygenOrSulfer == 1) {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, 0d);
                            } else {
                                iAtom.setProperty(MMFF94_FORMAL_CHARGE, -(terminalOxygenOrSulfer - 1) / (double) terminalOxygenOrSulfer);
                            }
                            break ATOM_TYPE_SWITCH;
                        }
                        if (neighborType == 77 && terminalOxygenOrSulfer != 0) {
                            iAtom.setProperty(MMFF94_FORMAL_CHARGE, -1d / (double) terminalOxygenOrSulfer);
                            break ATOM_TYPE_SWITCH;
                        }
                    }
                    break;
                }
                case 76: {
                    IRingSet rings = ((IRingSet) iAtom.getContainer().getProperty(MMFF94_RINGS)).getRings(iAtom);
                    if (!rings.isEmpty()) {
                        //TODO what if the atom was shared between two rings?
                        int nitrogens = 0;
                        for (IAtom atomInRing : rings.getAtomContainer(0).atoms()) {
                            if (atomInRing.getProperty(MMFF94_TYPE).equals(76)) {
                                nitrogens++;
                            }
                        }
                        iAtom.setProperty(MMFF94_FORMAL_CHARGE, -1d / (double) nitrogens);
                    }
                    break;
                }
                case 55:
                case 56:
                case 81: {
                    Arrays.fill(conjugates, false);
                    conjugates[iAtom.getIndex()] = true;
                    int conjugatesCount = 1;
                    int oldConjugates = 0;
                    int formalCharge = iAtom.getFormalCharge();
                    while (conjugatesCount > oldConjugates) {
                        oldConjugates = conjugatesCount;
                        for (int i = 0; i < container.getAtomCount(); i++) {
                            if (!conjugates[i]) {
                                continue;
                            }
                            for (IAtom neighbor : container.getConnectedAtomsList(container.getAtom(i))) {
                                int type = neighbor.getProperty(MMFF94_TYPE);
                                if ((type != 57) && (type != 80)) {
                                    continue;
                                }
                                for (IAtom neighbor2 : container.getConnectedAtomsList(neighbor)) {
                                    int type2 = neighbor2.getProperty(MMFF94_TYPE);
                                    if (type2 != 55 && type2 != 56 && type2 != 81) {
                                        continue;
                                    }
                                    if (!conjugates[neighbor2.getIndex()]) {
                                        conjugates[neighbor2.getIndex()] = true;
                                        formalCharge += neighbor2.getFormalCharge();
                                        conjugatesCount++;
                                    }
                                }
                            }
                        }
                    }
                    iAtom.setProperty(MMFF94_FORMAL_CHARGE, formalCharge / ((double) conjugatesCount));
                    break;
                }
                case 61: {
                    for (IAtom neighbor : iAtom.getContainer().getConnectedAtomsList(iAtom)) {
                        if (neighbor.getProperty(MMFF94_TYPE).equals(42)) {
                            iAtom.setProperty(MMFF94_FORMAL_CHARGE, 1d);
                            break ATOM_TYPE_SWITCH;
                        }
                    }
                    break;
                }
                case 34:
                case 49:
                case 51:
                case 54:
                case 58:
                case 92:
                case 93:
                case 94:
                case 97: {
                    iAtom.setProperty(MMFF94_FORMAL_CHARGE, 1d);
                    break;
                }
                case 87:
                case 95:
                case 96:
                case 98:
                case 99: {
                    iAtom.setProperty(MMFF94_FORMAL_CHARGE, 2d);
                    break;
                }
                case 88: {
                    iAtom.setProperty(MMFF94_FORMAL_CHARGE, 3d);
                    break;
                }
                case 35:
                case 62:
                case 89:
                case 90:
                case 91: {
                    iAtom.setProperty(MMFF94_FORMAL_CHARGE, -1d);
                    break;
                }

            }
        }
    }

    private void calculatePartialCharge(IAtomContainer container) {
        for (IAtom iAtom : container.atoms()) {
            double q0 = iAtom.getProperty(MMFF94_FORMAL_CHARGE);
            Float[] parameters = mmffParameters.partialChargeIncrement.get(((int) iAtom.getProperty(MMFF94_TYPE)) - 1);
            if (parameters[2] == 0d) {
                for (IAtom neighbor : container.getConnectedAtomsList(iAtom)) {
                    double k0 = neighbor.getProperty(MMFF94_FORMAL_CHARGE);
                    if (k0 < 0d) {
                        q0 += k0 / (2d * neighbor.getBondCount());
                    }
                }
            }
            if (iAtom.getProperty(MMFF94_TYPE).equals(62)) {
                for (IAtom neighbor : container.getConnectedAtomsList(iAtom)) {
                    double k0 = neighbor.getProperty(MMFF94_FORMAL_CHARGE);
                    if (k0 > 0d) {
                        q0 -= k0 / 2d;
                    }
                }
            }
            boolean foundParamter;
            int i, j;
            double qI = 0, qKSum = 0;
            for (IAtom atom : iAtom.getContainer().getConnectedAtomsList(iAtom)) {
                qKSum += (double) atom.getProperty(MMFF94_FORMAL_CHARGE);
                int bondType = findBondType(iAtom, atom);
                i = Math.min(iAtom.getProperty(MMFF94_TYPE), atom.getProperty(MMFF94_TYPE));
                j = Math.max(iAtom.getProperty(MMFF94_TYPE), atom.getProperty(MMFF94_TYPE));
                boolean reversed = (int) iAtom.getProperty(MMFF94_TYPE) > (int) atom.getProperty(MMFF94_TYPE);
                foundParamter = false;
                int index = binarySearch(mmffParameters.bondChargeIncrements, 1, i);
                if (index != -1) {
                    for (int z = index; z < mmffParameters.bondChargeIncrements.size(); z++) {
                        Float[] bondChargeIncrement = mmffParameters.bondChargeIncrements.get(z);
                        if (bondChargeIncrement[1] == i && bondChargeIncrement[2] <= j) {
                            if (bondChargeIncrement[2] == j && bondType == bondChargeIncrement[0]) {
                                qI += (bondChargeIncrement[3] * (reversed ? -1 : 1));
                                foundParamter = true;
                                break;
                            }
                        } else {
                            break;
                        }
                    }
                }
                if (!foundParamter) {
                    qI += (mmffParameters.partialChargeIncrement.get(((int) iAtom.getProperty(MMFF94_TYPE)) - 1)[1] - mmffParameters.partialChargeIncrement.get(((int) atom.getProperty(MMFF94_TYPE)) - 1)[1]);
                }
            }
            qI += (1 - iAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).crd * parameters[2]) * q0 + qKSum * parameters[2];
            iAtom.setCharge((double) qI);
        }
    }

    /**
     * calculates the angle bending energy for the angle that is formed from the
     * three atoms.
     *
     * @param iAtom the atom at one side
     * @param jAtom the atom in the middle
     * @param kAtom the atom at the other side
     * @return energy
     */
    public double calculateAngleBendingEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom) {
        if(iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom)==null){
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k");
        }
        checkParapetersAssigned(iAtom);
        Float[] parameters = findAngleBendingParameters(iAtom, jAtom, kAtom);
        double angle = calculateAngleBetween(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d());
        double deltaAngle = angle - parameters[5];
        if (jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear) {
            return 143.9325 * parameters[4] * (1 + Math.cos(Math.toRadians(angle)));
        } else {
            return Math.toRadians(Math.toRadians(143.9325)) * parameters[4] * .5 * deltaAngle * deltaAngle * (1 - 0.006981317 * deltaAngle);
        }
    }

    /**
     * calculates stretch bend energy that couples i-j and j-k bonds stretching
     * to the angle formed by i-j-k. This energy favors bonds stretching when
     * the angle decrease.
     *
     * @param iAtom first atom
     * @param jAtom middle atom
     * @param kAtom end atom
     * @return
     */
    public double calculateStretchBendEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom) {
        if(iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom)==null){
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k");
        }
        checkParapetersAssigned(iAtom);
        if (jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear) {
            return 0;
        }
        Float[] angleBendingParameter = findAngleBendingParameters(iAtom, jAtom, kAtom);
        double deltaAngle = calculateAngleBetween(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d()) - angleBendingParameter[5];
        MMFF94Parameters.StretchParameters ijStretch = iAtom.getContainer().getBond(iAtom, jAtom).<MMFF94Parameters.StretchParameters>getProperty(MMFF94_PARAMETER_STRETCH);
        MMFF94Parameters.StretchParameters jkStretch = iAtom.getContainer().getBond(jAtom, kAtom).<MMFF94Parameters.StretchParameters>getProperty(MMFF94_PARAMETER_STRETCH);
        double deltaRij = iAtom.getPoint3d().distance(jAtom.getPoint3d()) - ijStretch.r0;
        double deltaRjk = jAtom.getPoint3d().distance(kAtom.getPoint3d()) - jkStretch.r0;
        Float[] stretchBendParameter = findStretchBendParameters(iAtom, jAtom, kAtom);
        return 2.51210 * (stretchBendParameter[4] * deltaRij + stretchBendParameter[5] * deltaRjk) * deltaAngle;
    }

    private Vector3d ijVector = new Vector3d(), kjVector = new Vector3d(), ljVector = new Vector3d(), perpendicular = new Vector3d();

    /**
     * calculates out of plane energy for the atoms
     *
     * @param iAtom first atom forming the plane
     * @param jAtom the middle atom forming the plane
     * @param kAtom the last atom forming the plane
     * @param lAtom out of the plane atom
     * @return
     */
    public double calculateOutOfPlaneEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
        if(jAtom.getBond(iAtom) == null || jAtom.getBond(kAtom)==null && jAtom.getBond(lAtom)==null){
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k/l where l is bonded to j");
        }
        checkParapetersAssigned(iAtom);
        int[] ikl = new int[3];
        ikl[0] = iAtom.getProperty(MMFF94_TYPE);
        int j = jAtom.getProperty(MMFF94_TYPE);
        ikl[1] = kAtom.getProperty(MMFF94_TYPE);
        ikl[2] = lAtom.getProperty(MMFF94_TYPE);
        int equivalentIndex = 0;
        List<Float[]> outOfPlaneParameters = mmff94s ? mmffParameters.outOfPlaneParametersStatic : mmffParameters.outOfPlaneParameters;
        while (true) {
            Arrays.sort(ikl);
            int index = binarySearch(outOfPlaneParameters, 1, j);
            if (index != -1) {
                for (int z = index; z < outOfPlaneParameters.size(); z++) {
                    Float[] parameters = outOfPlaneParameters.get(z);
                    if (parameters[1] == j) {
                        if (((ikl[0] == parameters[0] && ikl[1] == parameters[2] && ikl[2] == parameters[3]))) {
                            ijVector.sub(iAtom.getPoint3d(), jAtom.getPoint3d());
                            ijVector.normalize();
                            kjVector.sub(kAtom.getPoint3d(), jAtom.getPoint3d());
                            kjVector.normalize();
                            perpendicular.cross(ijVector, kjVector);
                            ljVector.sub(lAtom.getPoint3d(), jAtom.getPoint3d());
                            ljVector.normalize();
                            double angle = Math.asin(perpendicular.dot(ljVector) / Math.sin(calculateAngleBetween(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d()) * Math.PI / 180)) * 180 / Math.PI;
                            return Math.toRadians(Math.toRadians(143.9325)) * .5 * parameters[4] * angle * angle;
                        }
                    }
                }
            }
            if (equivalentIndex > 3) {
                //no parameters where found. We can omit oop
                return 0;
            }
            ikl[0] = findAtomTypeEquivalence(iAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            ikl[1] = findAtomTypeEquivalence(kAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            ikl[2] = findAtomTypeEquivalence(lAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            equivalentIndex++;
        }
    }

    /**
     * calculates the torsion energy for the torsion angle formed by i-j-k-l
     * atoms
     *
     * @param iAtom
     * @param jAtom
     * @param kAtom
     * @param lAtom
     * @return
     */
    public double calculateTorsionEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
        if(iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom)==null && kAtom.getBond(lAtom)==null){
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k-l");
        }
        checkParapetersAssigned(iAtom);
        if ((Integer) kAtom.getProperty(MMFF94_TYPE) < (Integer) jAtom.getProperty(MMFF94_TYPE)
                || (jAtom.getProperty(MMFF94_TYPE).equals((Integer) kAtom.getProperty(MMFF94_TYPE)) && (Integer) iAtom.getProperty(MMFF94_TYPE) > (Integer) lAtom.getProperty(MMFF94_TYPE))) {
            IAtom temp = jAtom;
            jAtom = kAtom;
            kAtom = temp;
            temp = iAtom;
            iAtom = lAtom;
            lAtom = temp;
        }
        Vector3d vectorIJ = new Vector3d(iAtom.getPoint3d().x - jAtom.getPoint3d().x, iAtom.getPoint3d().y - jAtom.getPoint3d().y, iAtom.getPoint3d().z - jAtom.getPoint3d().z);
        Vector3d vectorKJ = new Vector3d(kAtom.getPoint3d().x - jAtom.getPoint3d().x, kAtom.getPoint3d().y - jAtom.getPoint3d().y, kAtom.getPoint3d().z - jAtom.getPoint3d().z);
        Vector3d vectorJK = new Vector3d(vectorKJ);
        vectorJK.negate();
        Vector3d vectorLK = new Vector3d(lAtom.getPoint3d().x - kAtom.getPoint3d().x, lAtom.getPoint3d().y - kAtom.getPoint3d().y, lAtom.getPoint3d().z - kAtom.getPoint3d().z);
        double torsionAngle = calculateTorsionAngle(vectorIJ, vectorKJ, vectorJK, vectorLK);//in radians
        Float[] parameters = findTorsionParameters(iAtom, jAtom, kAtom, lAtom);
        return 0.5 * (parameters[5] * (1 + Math.cos(torsionAngle)) + parameters[6] * (1 - Math.cos(2 * torsionAngle)) + parameters[7] * (1 + Math.cos(3 * torsionAngle)));
    }

    /**
     * atoms must be sorted
     *
     * @param iAtom
     * @param jAtom
     * @param kAtom
     * @param lAtom
     * @return
     */
    private Float[] findTorsionParameters(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
        int i = iAtom.getProperty(MMFF94_TYPE);
        int j = jAtom.getProperty(MMFF94_TYPE);
        int k = kAtom.getProperty(MMFF94_TYPE);
        int l = lAtom.getProperty(MMFF94_TYPE);
        int ijBondType = findBondType(iAtom, jAtom);
        int jkBondType = findBondType(jAtom, kAtom);
        int klBondType = findBondType(kAtom, lAtom);
        int torsionType = jkBondType;
        if (jkBondType == 0 && jAtom.getBond(kAtom).getOrder().equals(IBond.Order.SINGLE) && (ijBondType == 1 || klBondType == 1)) {
            torsionType = 2;
        }
        if (iAtom.isInRing() && jAtom.isInRing() && kAtom.isInRing() && lAtom.isInRing()) {
            IRingSet rings = iAtom.getContainer().getProperty(MMFF94_RINGS);
            IRingSet irings = rings.getRings(iAtom);
            for (IAtomContainer ring : irings.atomContainers()) {
                if (ring.contains(jAtom) && ring.contains(kAtom) && ring.contains(lAtom)) {
                    if (ring.getAtomCount() == 4) {
                        torsionType = 4;
                        break;
                    } else if (ring.getAtomCount() == 5 && (i == 1 || j == 1 || k == 1 || l == 1)) {
                        torsionType = 5;
                        break;
                    }
                }
            }
        }
        int equivalentIndex = 0, iIndex = 0, lIndex = 0;
        List<Float[]> torsionParameters = mmff94s ? mmffParameters.torsionParametersStatic : mmffParameters.torsionParameters;
        while (true) {
            int index = binarySearch(torsionParameters, 2, j);
            if (index != -1) {
                for (int z = index; z < torsionParameters.size(); z++) {
                    Float[] parameters = torsionParameters.get(z);
                    if (parameters[0] == torsionType && parameters[2] == j && parameters[3] == k) {
                        if (i == parameters[1] && l == parameters[4]) {
                            return parameters;
                        }
                    }
                }
            }
            if (equivalentIndex > 3) {
                //empirical rules
                Float[] paramters = new Float[]{(float) torsionType,
                    ((Integer) iAtom.getProperty(MMFF94_TYPE)).floatValue(),
                    ((Integer) jAtom.getProperty(MMFF94_TYPE)).floatValue(),
                    ((Integer) kAtom.getProperty(MMFF94_TYPE)).floatValue(),
                    ((Integer) lAtom.getProperty(MMFF94_TYPE)).floatValue(), 0f, 0f, 0f};
                float[] jEmpiricalParamters = findTorsionEmpiricalParamters(jAtom);
                float[] kEmpiricalParamters = findTorsionEmpiricalParamters(kAtom);
                MMFF94Parameters.GeometricParameters jGeometricProperties = jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES);
                MMFF94Parameters.GeometricParameters kGeometricProperties = kAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES);
                if (jGeometricProperties.ideallyLinear || kGeometricProperties.ideallyLinear) {
                    //rule (a)
                    paramters[5] = 0f;
                    paramters[6] = 0f;
                    paramters[7] = 0f;
                    return paramters;
                } else if (jGeometricProperties.aromatic && kGeometricProperties.aromatic && ((IRingSet) jAtom.getContainer().getProperty(MMFF94_RINGS)).getRings(jAtom).contains(kAtom)) {
                    //rule (b)
                    float pi = .3f, beta = 6;
                    if (!jGeometricProperties.piLonePair && !kGeometricProperties.piLonePair) {
                        pi = .5f;
                    }
                    if ((jGeometricProperties.valence == 3 && kGeometricProperties.valence == 4)
                            || (jGeometricProperties.valence == 4 && kGeometricProperties.valence == 3)) {
                        beta = 3;
                    }
                    paramters[6] = (float) (beta * pi * Math.sqrt(jEmpiricalParamters[1] * kEmpiricalParamters[1]));
                } else if (jAtom.getBond(kAtom).getOrder().equals(IBond.Order.DOUBLE)) {
                    //rule (c)
                    float pi = .4f;
                    if (jGeometricProperties.mltb == 2 && kGeometricProperties.mltb == 2) {
                        pi = 1;
                    }
                    paramters[6] = (float) (6 * pi * Math.sqrt(jEmpiricalParamters[1] * kEmpiricalParamters[1]));
                } else if (jGeometricProperties.crd == 4 && kGeometricProperties.crd == 4) {
                    //rule (d)
                    paramters[7] = (float) (Math.sqrt(jEmpiricalParamters[2] * kEmpiricalParamters[2]) / ((jGeometricProperties.crd - 1) * (kGeometricProperties.crd - 1)));
                } else if (jGeometricProperties.crd == 4 && kGeometricProperties.crd != 4 || jGeometricProperties.crd != 4 && kGeometricProperties.crd == 4) {
                    //rule (e/f)
                    int nonTetraCoordinate = jGeometricProperties.crd == 4 ? k : j;
                    MMFF94Parameters.GeometricParameters nonTetraCoordinateGeometricProperties = nonTetraCoordinate == j ? jGeometricProperties : kGeometricProperties;
                    if (nonTetraCoordinateGeometricProperties.crd == 3
                            && ((nonTetraCoordinateGeometricProperties.valence == 4
                            || nonTetraCoordinateGeometricProperties.valence == 34)
                            || nonTetraCoordinateGeometricProperties.mltb != 0)) {
                        paramters[5] = 0f;
                        paramters[6] = 0f;
                        paramters[7] = 0f;
                        return paramters;
                    } else if (nonTetraCoordinateGeometricProperties.crd == 2
                            && (nonTetraCoordinateGeometricProperties.valence == 3
                            || nonTetraCoordinateGeometricProperties.mltb != 0)) {
                        paramters[5] = 0f;
                        paramters[6] = 0f;
                        paramters[7] = 0f;
                        return paramters;
                    } else {
                        paramters[7] = (float) (Math.sqrt(jEmpiricalParamters[2] * kEmpiricalParamters[2]) / ((jGeometricProperties.crd - 1) * (kGeometricProperties.crd - 1)));
                    }
                } else if ((jAtom.getBond(kAtom).getOrder().equals(IBond.Order.SINGLE) && jGeometricProperties.mltb != 0 && kGeometricProperties.mltb != 0)
                        || (jGeometricProperties.mltb != 0 && kGeometricProperties.piLonePair)
                        || (jGeometricProperties.piLonePair && kGeometricProperties.mltb != 0)) {
                    //rule (g)
                    //case 1
                    if (jGeometricProperties.piLonePair && kGeometricProperties.piLonePair) {
                        paramters[5] = 0f;
                        paramters[6] = 0f;
                        paramters[7] = 0f;
                        return paramters;
                    } else if (jGeometricProperties.piLonePair && kGeometricProperties.mltb != 0) {
                        //case 2
                        float pi = .15f;
                        if (jGeometricProperties.mltb == 1) {
                            pi = .5f;
                        } else if (getPeriodicTableRow(jAtom) == 1 || getPeriodicTableRow(kAtom) == 1) {
                            //beolong to row 2 in the periodic table
                            pi = .3f;
                        }
                        paramters[6] = (float) (6 * pi * Math.sqrt(jEmpiricalParamters[1] * kEmpiricalParamters[1]));
                    } else if (kGeometricProperties.piLonePair && jGeometricProperties.mltb != 0) {
                        //case 3
                        float pi = .15f;
                        if (kGeometricProperties.mltb == 1) {
                            pi = .5f;
                        } else if (getPeriodicTableRow(jAtom) == 1 || getPeriodicTableRow(kAtom) == 1) {
                            //beolong to row 2 in the periodic table
                            pi = .3f;
                        }
                        paramters[6] = (float) (6 * pi * Math.sqrt(jEmpiricalParamters[1] * kEmpiricalParamters[1]));
                    } else if ((jGeometricProperties.mltb == 1 || kGeometricProperties.mltb == 1) && (jAtom.getAtomicNumber() != 6 || kAtom.getAtomicNumber() != 6)) {
                        //case 4
                        paramters[6] = (float) (6 * .4 * Math.sqrt(jEmpiricalParamters[1] * kEmpiricalParamters[1]));
                    } else {
                        //case 5
                        paramters[6] = (float) (6 * .15 * Math.sqrt(jEmpiricalParamters[1] * kEmpiricalParamters[1]));
                    }
                } else if ((jAtom.getAtomicNumber() == 8 || jAtom.getAtomicNumber() == 16)
                        && (kAtom.getAtomicNumber() == 8 || kAtom.getAtomicNumber() == 16)) {
                    //rule (h)
                    float WW = jAtom.getAtomicNumber() == 8 ? 2 : 8;
                    WW *= kAtom.getAtomicNumber() == 8 ? 2 : 8;
                    paramters[6] = (float) -Math.sqrt(WW);
                } else {
                    paramters[7] = (float) Math.sqrt(jEmpiricalParamters[2] * kEmpiricalParamters[2]) / ((jGeometricProperties.crd - 1) * (kGeometricProperties.crd - 1));
                }
                return paramters;
            }
            if (equivalentIndex == 3 && torsionType == 5) {
                if (jkBondType == 1) {
                    torsionType = 1;
                    equivalentIndex = 0;
                    i = iAtom.getProperty(MMFF94_TYPE);
                    k = kAtom.getProperty(MMFF94_TYPE);
                    continue;
                } else if (jkBondType == 0 && (ijBondType == 1 || klBondType == 1)) {
                    torsionType = 2;
                    equivalentIndex = 0;
                    i = iAtom.getProperty(MMFF94_TYPE);
                    k = kAtom.getProperty(MMFF94_TYPE);
                    continue;
                }
            }
            iIndex = equivalentIndex;
            lIndex = equivalentIndex;
            if (equivalentIndex == 1) {
                iIndex = 1;
                lIndex = 3;
            } else if (equivalentIndex == 2) {
                iIndex = 3;
                lIndex = 1;
            }
            i = findAtomTypeEquivalence(iAtom.getProperty(MMFF94_TYPE), iIndex);
            l = findAtomTypeEquivalence(lAtom.getProperty(MMFF94_TYPE), lIndex);
            equivalentIndex++;
        }
    }

    private float[] findTorsionEmpiricalParamters(IAtom iAtom) {
        for (float[] torsionEmpiricalParamter : mmffParameters.torsionEmpiricalParamters) {
            if (((int) torsionEmpiricalParamter[0]) == iAtom.getAtomicNumber()) {
                return torsionEmpiricalParamter;
            }
        }
        throw new RuntimeException("no empirical torsion paramters was found for " + iAtom.getAtomTypeName());
    }

    private Vector3d n1 = new Vector3d(), n2 = new Vector3d();

    private float calculateTorsionAngle(Vector3d ijVector, Vector3d kjVector, Vector3d jkVector, Vector3d klVector) {
        n1.cross(ijVector, kjVector);
        n2.cross(jkVector, klVector);
        return (float) Math.acos(clampToOne(n1.dot(n2) / (n1.length() * n2.length())));
    }

    private double clampToOne(double n) {
        return n > 1 ? 1 : (n < -1 ? -1 : n);
    }

    private int findBondType(IAtom iAtom, IAtom jAtom) {
        IBond bond = iAtom.getContainer().getBond(iAtom, jAtom);
        Boolean aromaticBond = bond.<Boolean>getProperty("mmff.arom");
        if (bond.getOrder().equals(IBond.Order.SINGLE) && (aromaticBond == null || aromaticBond != true)) {
            MMFF94Parameters.GeometricParameters iAtomGeometricProperties = iAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES);
            MMFF94Parameters.GeometricParameters jAtomGeometricProperties = jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES);
            if (iAtomGeometricProperties.singleBondBetweenMultipleBonds
                    && jAtomGeometricProperties.singleBondBetweenMultipleBonds) {
                return 1;
            } else if (iAtomGeometricProperties.aromatic
                    && jAtomGeometricProperties.aromatic) {
                return 1;
            }
        }
        return 0;
    }

    private int findAngleType(IAtom iAtom, IAtom jAtom, IAtom kAtom) {
        int bondType = findBondType(iAtom, jAtom) + findBondType(jAtom, kAtom);
        if (iAtom.isInRing() && jAtom.isInRing() && kAtom.isInRing()) {
            IRingSet rings = iAtom.getContainer().getProperty(MMFF94_RINGS);
            IRingSet iRings = rings.getRings(iAtom);
            for (int ring = 0; ring < iRings.getAtomContainerCount(); ring++) {
                if (iRings.getAtomContainer(ring).contains(jAtom) && iRings.getAtomContainer(ring).contains(kAtom)) {
                    if (iRings.getAtomContainer(ring).getAtomCount() == 3) {
                        return bondType == 0 ? 3 : (3 + bondType + 1);
                    } else if (iRings.getAtomContainer(ring).getAtomCount() == 4) {
                        return bondType == 0 ? 4 : (4 + bondType + 2);
                    }
                }
            }
        }
        return bondType;
    }

    private int findAtomTypeEquivalence(int typeName, int index) {
        Integer[] equivalentTypes = mmffParameters.atomTypeEquivalenceParameters.get(typeName - (typeName >= 87 ? 5 : 1));
        if (equivalentTypes == null) {
            return 0;
        } else if (index >= 4) {
            throw new ArrayIndexOutOfBoundsException();
        }
        return equivalentTypes[index];
    }

    private double calculateAngleBetween(Point3d i, Point3d j, Point3d k) {
        double ux = i.x - j.x;
        double uy = i.y - j.y;
        double uz = i.z - j.z;
        double vx = k.x - j.x;
        double vy = k.y - j.y;
        double vz = k.z - j.z;
        return Math.acos(clampToOne((ux * vx + uy * vy + uz * vz) / (i.distance(j) * k.distance(j)))) * 180 / Math.PI;
    }

    private Float[] findAngleBendingParameters(IAtom iAtom, IAtom jAtom, IAtom kAtom) {
        int i = iAtom.getProperty(MMFF94_TYPE);
        int j = jAtom.getProperty(MMFF94_TYPE);
        int k = kAtom.getProperty(MMFF94_TYPE);
        int angleType = findAngleType(iAtom, jAtom, kAtom);
        int tempI, tempK, equivalentIndex = 0;
        while (true) {
            tempI = i;
            tempK = k;
            i = Math.min(tempI, tempK);
            k = Math.max(tempI, tempK);
            int index = binarySearch(mmffParameters.angleBendingParameters, 2, j);
            if (index != -1) {
                for (int z = index; z < mmffParameters.angleBendingParameters.size(); z++) {
                    Float[] parameters = mmffParameters.angleBendingParameters.get(z);
                    if (parameters[0] == angleType && parameters[1] == i && parameters[2] == j && parameters[3] == k) {
                        if (parameters[4] == .0) {
                            //use empirical parameters
                            return calculateEmpiricalAngleBendingParameters(iAtom, jAtom, kAtom, angleType, parameters[5]);
                        } else {
                            return parameters;
                        }
                    }
                }
            }
            if (equivalentIndex > 3) {
                //no paramters available. we will use the empirical rule
                if (debug) {
                    System.out.println("Warning: Empirical Angle Bending");
                }
                return calculateEmpiricalAngleBendingParameters(iAtom, jAtom, kAtom, angleType, null);
            }
            i = findAtomTypeEquivalence(iAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            k = findAtomTypeEquivalence(kAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            equivalentIndex++;
        }
    }

    private Float[] calculateEmpiricalAngleBendingParameters(IAtom iAtom, IAtom jAtom, IAtom kAtom, int angleType, Float idealAngle) {
        Float[] parameters = new Float[6];
        parameters[0] = (float) angleType;
        parameters[1] = ((Integer) iAtom.getProperty(MMFF94_TYPE)).floatValue();
        parameters[2] = ((Integer) jAtom.getProperty(MMFF94_TYPE)).floatValue();
        parameters[3] = ((Integer) kAtom.getProperty(MMFF94_TYPE)).floatValue();
        int bondCount = jAtom.getBondCount();
        parameters[5] = 120f;
        if (bondCount == 4) {
            parameters[5] = 109.45f;
        } else if (bondCount == 3) {
            MMFF94Parameters.GeometricParameters jGeometricProperties = jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES);
            if (jGeometricProperties.valence == 3 && jGeometricProperties.mltb == 0) {
                if (jAtom.getAtomicNumber() == 7) {
                    //nitrogen
                    parameters[5] = 107f;
                } else {
                    parameters[5] = 92f;
                }
            }
        } else if (bondCount == 2) {
            if (jAtom.getAtomicNumber() == 8) {
                //oxygen
                parameters[5] = 105f;
            } else if (jAtom.getAtomicNumber() > 10) {
                parameters[5] = 95f;
            } else if (jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear) {
                parameters[5] = 180f;
            }
        }
        float beta = 1.75f;
        if (iAtom.isInRing() && jAtom.isInRing() && kAtom.isInRing()) {
            IRingSet rings = iAtom.getContainer().getProperty(MMFF94_RINGS);
            IRingSet iRings = rings.getRings(iAtom);
            for (int ring = 0; ring < iRings.getAtomContainerCount(); ring++) {
                if (iRings.getAtomContainer(ring).contains(jAtom) && iRings.getAtomContainer(ring).contains(kAtom)) {
                    if (iRings.getAtomContainer(ring).getAtomCount() == 3) {
                        parameters[5] = 60f;
                        beta *= .05f;
                    } else if (iRings.getAtomContainer(ring).getAtomCount() == 4) {
                        parameters[5] = 90f;
                        beta *= .85f;
                    }
                    break;
                }
            }
        }
        if (idealAngle != null) {
            parameters[5] = idealAngle;
        }
        float rIJ = iAtom.getContainer().getBond(iAtom, jAtom).<MMFF94Parameters.StretchParameters>getProperty(MMFF94_PARAMETER_STRETCH).r0;
        float rJK = iAtom.getContainer().getBond(jAtom, kAtom).<MMFF94Parameters.StretchParameters>getProperty(MMFF94_PARAMETER_STRETCH).r0;
        float D = (float) (Math.pow(rIJ - rJK, 2) / Math.pow(rIJ + rJK, 2));
        float Zi = 0, Cj = 0, Zk = 0;
        int found = 0;
        for (int l = 0; l < mmffParameters.angleBendingEmpiricalParamters.length; l++) {
            float[] angleBendingEmpiricalParamter = mmffParameters.angleBendingEmpiricalParamters[l];
            if (angleBendingEmpiricalParamter[0] == iAtom.getAtomicNumber()) {
                Zi = angleBendingEmpiricalParamter[1];
                found++;
            }
            if (angleBendingEmpiricalParamter[0] == jAtom.getAtomicNumber()) {
                Cj = angleBendingEmpiricalParamter[2];
                found++;
            }
            if (angleBendingEmpiricalParamter[0] == kAtom.getAtomicNumber()) {
                Zk = angleBendingEmpiricalParamter[1];
                found++;
            }
            if (found == 3) {
                break;
            }
        }
        if (found != 3) {
            throw new RuntimeException(String.format("could not generate empirical angle bending parameters for atoms %s-%s-%s", iAtom.getAtomTypeName(), jAtom.getAtomTypeName(), kAtom.getAtomTypeName()));
        }
        parameters[4] = (float) (beta * Zi * Cj * Zk / ((rIJ + rJK) * Math.pow(Math.toRadians(parameters[5]), 2) * Math.exp(2 * D)));
        return parameters;
    }

    private Float[] findStretchBendParameters(IAtom iAtom, IAtom jAtom, IAtom kAtom) {
        int i = iAtom.getProperty(MMFF94_TYPE);
        int j = jAtom.getProperty(MMFF94_TYPE);
        int k = kAtom.getProperty(MMFF94_TYPE);
        int stretchBendType = 0;
        int angleType = findAngleType(iAtom, jAtom, kAtom);
        int ijBondType = findBondType(iAtom, jAtom);
        int jkBondType = findBondType(jAtom, kAtom);
        boolean swap = false;
        if (i == k && ijBondType < jkBondType) {
            int temp = jkBondType;
            jkBondType = ijBondType;
            ijBondType = temp;
            swap = true;
        }
        switch (angleType) {
            case 1: {
                stretchBendType = ((ijBondType != 0 || (ijBondType == jkBondType)) ? 1 : 2);
                break;
            }

            case 2: {
                stretchBendType = 3;
                break;
            }
            case 4: {
                stretchBendType = 4;
                break;
            }

            case 3: {
                stretchBendType = 5;
                break;
            }

            case 5: {
                stretchBendType = ((ijBondType != 0 || (ijBondType == jkBondType)) ? 6 : 7);
                break;
            }

            case 6: {
                stretchBendType = 8;
                break;
            }

            case 7: {
                stretchBendType = ((ijBondType != 0 || (ijBondType == jkBondType)) ? 9 : 10);
                break;
            }

            case 8: {
                stretchBendType = 11;
                break;
            }
        }
        for (Float[] parameters : mmffParameters.stretchBendParameters) {
            if (parameters[2] == j && stretchBendType == parameters[0] && i == parameters[1] && k == parameters[3]) {
                if (swap) {
                    return new Float[]{parameters[0], parameters[1], parameters[2], parameters[3], parameters[5], parameters[4]};
                } else {
                    return parameters;
                }
            }
        }
        //no paramters available. we will use the empirical rule
        int iRow = Math.min(getPeriodicTableRow(iAtom), getPeriodicTableRow(kAtom));
        int jRow = getPeriodicTableRow(jAtom);
        int kRow = Math.max(getPeriodicTableRow(iAtom), getPeriodicTableRow(kAtom));
        swap = getPeriodicTableRow(iAtom) > getPeriodicTableRow(kAtom);
        for (Float[] stretchBendEmpiricalParamter : mmffParameters.stretchBendEmpiricalParamters) {
            if (iRow == stretchBendEmpiricalParamter[0] && jRow == stretchBendEmpiricalParamter[1] && kRow == stretchBendEmpiricalParamter[2]) {
                if (swap) {
                    return new Float[]{(float) stretchBendType, (float) i, (float) j, (float) k, stretchBendEmpiricalParamter[4], stretchBendEmpiricalParamter[3]};
                } else {
                    return new Float[]{(float) stretchBendType, (float) i, (float) j, (float) k, stretchBendEmpiricalParamter[3], stretchBendEmpiricalParamter[4]};
                }
            }
        }
        throw new RuntimeException(String.format("could not find or generate stretch-bend parameters for %s-%s-%s", i, j, k));
    }

    private int getPeriodicTableRow(IAtom atom) {
        if (atom.getAtomicNumber() <= 2) {
            return 0;
        } else if (atom.getAtomicNumber() <= 10) {
            return 1;
        } else if (atom.getAtomicNumber() <= 18) {
            return 2;
        } else if (atom.getAtomicNumber() <= 36) {
            return 3;
        } else if (atom.getAtomicNumber() <= 54) {
            return 4;
        } else if (atom.getAtomicNumber() <= 86) {
            return 5;
        } else {
            return 6;
        }
    }

    private int getHerschbachLauriePeriodicTableRow(IAtom atom) {
        int atomicNumber = atom.getAtomicNumber();
        if (atomicNumber == 2) {
            return 1;
        } else if (atomicNumber >= 3 && atomicNumber <= 10) {
            return 2;
        } else if (atomicNumber >= 11 && atomicNumber <= 18) {
            return 3;
        } else if (atomicNumber >= 19 && atomicNumber <= 36) {
            return (atomicNumber >= 21 && atomicNumber <= 30) ? 40 : 4;
        } else if (atomicNumber >= 37 && atomicNumber <= 54) {
            return (atomicNumber >= 39) && (atomicNumber <= 48) ? 50 : 5;
        }
        throw new UnsupportedOperationException(String.format("heavy atoms (%s) are not supported", atom.getAtomTypeName()));
    }

    private boolean isSeparatedByNOrMoreBonds(IAtom atom1, IAtom atom2, int n) {
        HashSet<IAtom> searched = new HashSet();
        HashSet<IAtom> lastDepth = new HashSet();
        lastDepth.add(atom1);
        for (int depth = 1; depth < n; depth++) {
            searched.addAll(lastDepth);
            HashSet<IAtom> lastDepthClone = (HashSet<IAtom>) lastDepth.clone();
            for (IAtom atom : lastDepthClone) {
                lastDepth.addAll(atom1.getContainer().getConnectedAtomsList(atom));
            }
            lastDepth.removeAll(searched);
        }
        searched.addAll(lastDepth);
        return !searched.contains(atom2);
    }

    private boolean isSeparatedByNBonds(IAtom atom1, IAtom atom2, int n) {
        HashSet<IAtom> searched = new HashSet();
        HashSet<IAtom> lastDepth = new HashSet();
        lastDepth.add(atom1);
        for (int depth = 1; depth <= n; depth++) {
            searched.addAll(lastDepth);
            HashSet<IAtom> lastDepthClone = (HashSet<IAtom>) lastDepth.clone();
            for (IAtom atom : lastDepthClone) {
                lastDepth.addAll(atom1.getContainer().getConnectedAtomsList(atom));
            }
            lastDepth.removeAll(searched);
        }
        return lastDepth.contains(atom2);
    }

    private static int binarySearch(List<Float[]> array, int keyIndex, float key) {
        int low = 0;
        int high = array.size() - 1;

        while (low <= high) {
            int mid = low + (high - low) / 2;
            float midValue = array.get(mid)[keyIndex];
            if (midValue > key) {
                high = mid - 1;
            } else if (midValue < key) {
                low = mid + 1;
            } else if (midValue == key) {
                //sihft back if the previous record has the key
                while (mid > 0 && array.get(mid - 1)[keyIndex] == key) {
                    mid--;
                }
                return mid;
            }
        }
        return -1;  // key not found.
    }

}
