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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Objects;
import java.util.Queue;
import org.jme.forcefield.EnergyComponent;
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
    static final String MMFF94_NON_BONDED_INTERACTION = "mmff.nonBondedInteraction";
    static final String MMFF94_VDW_INTERACTION = "mmff.vdwInteraction";
    static final String MMFF94_ANGLE_BENDING = "mmff.angleBending";
    static final String MMFF94_STRETCH_BEND = "mmff.stretchBend";
    static final String MMFF94_OUT_OF_PLANE = "mmff.outOfPlane";
    static final String MMFF94_TORSION = "mmff.torsion";
    private static final Logger LOGGER = Logger.getLogger(MMFF94.class.getName());

    private final Mmff mmff = new Mmff();
    private MMFF94Parameters mmffParameters;
    private boolean mmff94s;
    private MMFF94BondStretchingComponent bondStretchingComponent = new MMFF94BondStretchingComponent();
    private MMFF94StretchBendComponent stretchBendComponent = new MMFF94StretchBendComponent();
    private MMFF94AngleBendingComponent angleBendingComponent = new MMFF94AngleBendingComponent();
    private MMFF94OutOfPlaneComponent outOfPlaneComponent = new MMFF94OutOfPlaneComponent();
    private MMFF94TorsionComponent torsionComponent = new MMFF94TorsionComponent();
    private MMFF94VdwComponent vdwComponent = new MMFF94VdwComponent();
    private MMFF94ElectrostaticComponent electrostaticComponent = new MMFF94ElectrostaticComponent();
    private List<EnergyComponent> energyComponents;

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
        energyComponents = new ArrayList<>(Arrays.asList(bondStretchingComponent, stretchBendComponent, angleBendingComponent, outOfPlaneComponent, torsionComponent, vdwComponent, electrostaticComponent));
    }

    @Override
    public List<EnergyComponent> getEnergyComponents() {
        return energyComponents;
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
            if (Objects.equals(bond.getProperty("mmff.arom"), true)) {
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
        if (LOGGER.isLoggable(Level.FINE)) {
            for (IAtom atom : atomContainer.atoms()) {
                LOGGER.fine(atom.getAtomTypeName() + ":\t" + atom.getCharge());
            }
        }
        assignBondParameters(atomContainer);

        HashMap<List<IAtom>, Boolean> nonBondedInteraction = new HashMap<>();
        HashMap<List<IAtom>, Float[]> angleBending = new HashMap<>();
        HashMap<List<IAtom>, Float[]> stretchBend = new HashMap<>();
        HashMap<List<IAtom>, Float[]> outOfPlane = new HashMap<>();
        HashMap<List<IAtom>, Float[]> torsion = new HashMap<>();
        for (IAtom atom : atomContainer.atoms()) {
            int atomType = atom.getProperty(MMFF94_TYPE);
            atom.setProperty(MMFF94_VDW_INTERACTION, mmffParameters.vdwParameters.get(atomType - (atomType >= 87 ? 5 : 1)));
            //non bonded interaction
            for (IAtom atom2 : atomContainer.atoms()) {
                if (atom.getIndex() >= atom2.getIndex()) {
                    continue;
                }
                boolean[] separatedBy3OrMore = isSeparatedByNOrMoreBonds(atom, atom2, 3);
                if ((separatedBy3OrMore[0] || separatedBy3OrMore[1])) {
                    nonBondedInteraction.put(Arrays.asList(atom, atom2), separatedBy3OrMore[0]);
                }
            }
            //bonded interaction
            if (atom.getBondCount() >= 2) {
                for (IBond bond2 : atom.bonds()) {
                    IAtom atom2 = bond2.getOther(atom);
                    for (IBond bond3 : atom.bonds()) {
                        IAtom atom3 = bond3.getOther(atom);
                        if (atom3 == atom2) {
                            continue;
                        }
                        int j = atomType;
                        int i = atom2.getProperty(MMFF94_TYPE);
                        int k = atom3.getProperty(MMFF94_TYPE);
                        if (i < k || (i == k && atom2.getIndex() < atom3.getIndex())) {
                            angleBending.put(Arrays.asList(atom2, atom, atom3), findAngleBendingParameters(atom2, atom, atom3));
                            stretchBend.put(Arrays.asList(atom2, atom, atom3), findStretchBendParameters(atom2, atom, atom3));
                        }
                        if (atom.getBondCount() >= 3) {
                            for (IBond bond4 : atom.bonds()) {
                                IAtom atom4 = bond4.getOther(atom);
                                if (atom4 == atom2 || atom4 == atom3) {
                                    continue;
                                }
                                if (i < k || (k == i && atom2.getIndex() < atom3.getIndex())) {
                                    outOfPlane.put(Arrays.asList(atom2, atom, atom3, atom4), findOutOfPlaneParameters(atom2, atom, atom3, atom4));
                                }
                            }
                        }
                        //torsion
                        if (!atom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear && !atom3.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear && atom3.getBondCount() >= 2) {
                            for (IBond bond4 : atom3.bonds()) {
                                IAtom atom4 = bond4.getOther(atom3);
                                if (atom4 == atom || atom4 == atom2) {
                                    continue;
                                }
                                int l = atom4.getProperty(MMFF94_TYPE);
                                if (j < k || (j == k && i < l) || (j == k && i == l && atom.getIndex() < atom3.getIndex())) {
                                    torsion.put(Arrays.asList(atom2, atom, atom3, atom4), findTorsionParameters(atom2, atom, atom3, atom4));
                                }
                            }
                        }
                    }
                }
            }
        }
        atomContainer.setProperty(MMFF94_NON_BONDED_INTERACTION, nonBondedInteraction);
        atomContainer.setProperty(MMFF94_ANGLE_BENDING, angleBending);
        atomContainer.setProperty(MMFF94_STRETCH_BEND, stretchBend);
        atomContainer.setProperty(MMFF94_OUT_OF_PLANE, outOfPlane);
        atomContainer.setProperty(MMFF94_TORSION, torsion);
    }

    private void assignBondParameters(IAtomContainer atomContainer) {
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
            int index = Collections.binarySearch(mmffParameters.bondStretchParameters, new MMFF94Parameters.StretchParameters(bondType, i, j, 0, 0), (MMFF94Parameters.StretchParameters p1, MMFF94Parameters.StretchParameters p2) -> {
                if (p1.i != p2.i) {
                    return Integer.compare(p1.i, p2.i);
                } else if (p1.j != p2.j) {
                    return Integer.compare(p1.j, p2.j);
                } else {
                    return Integer.compare(p1.bondType, p2.bondType);
                }
            });
            if (index >= 0) {
                MMFF94Parameters.StretchParameters parameters = mmffParameters.bondStretchParameters.get(index);
                if (parameters.i == i && parameters.j == j && parameters.bondType == bondType) {
                    bond.setProperty(MMFF94_PARAMETER_STRETCH, parameters);
                    found = true;
                }
            }
            if (!found) {
                bond.setProperty(MMFF94_PARAMETER_STRETCH, calculateEmpiricalStretchParameters(iAtom, jAtom));
            }
        }
    }

    private Float[] findOutOfPlaneParameters(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
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
                            return parameters;
                        }
                    }
                }
            }
            if (equivalentIndex > 3) {
                //no parameters where found. We can omit oop
                return null;
            }
            ikl[0] = findAtomTypeEquivalence(iAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            ikl[1] = findAtomTypeEquivalence(kAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            ikl[2] = findAtomTypeEquivalence(lAtom.getProperty(MMFF94_TYPE), equivalentIndex);
            equivalentIndex++;
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
        int index = Collections.binarySearch(mmffParameters.bondStretchEmpiricalParameterses, new MMFF94Parameters.StretchEmpiricalParameters(atomicNum1, atomicNum2, 0, 0), (MMFF94Parameters.StretchEmpiricalParameters p1, MMFF94Parameters.StretchEmpiricalParameters p2) -> {
            if (p1.i != p2.i) {
                return Integer.compare(p1.i, p2.i);
            } else {
                return Integer.compare(p1.j, p2.j);
            }
        });
        if (index >= 0) {
            MMFF94Parameters.StretchEmpiricalParameters parameters = mmffParameters.bondStretchEmpiricalParameterses.get(index);
            if (parameters.i == atomicNum1 && parameters.j == atomicNum2) {
                mmffBndkParams = parameters;
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
        float[] r0_i = {radiiElectronegativityParameters[0][1], radiiElectronegativityParameters[1][1]};

        float r0 = (float) (r0_i[0] + r0_i[1] - c * Math.pow(Math.abs(radiiElectronegativityParameters[0][2] - radiiElectronegativityParameters[1][2]), n));
        float kb;
        if (mmffBndkParams != null) {
            kb = (float) (mmffBndkParams.kb * Math.pow(mmffBndkParams.r0 / r0, 6));
        } else {
            Float[] HerschbachLaurieParameters = findHerschbachLaurieParameters(iAtom, jAtom);
            kb = (float) Math.pow(10.0, -(r0 - HerschbachLaurieParameters[2]) / HerschbachLaurieParameters[3]);
        }

        return new MMFF94Parameters.StretchParameters(findBondType(iAtom, jAtom), iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), kb, r0);
    }

    private Float[] findHerschbachLaurieParameters(IAtom iAtom, IAtom jAtom) {
        int row1 = getHerschbachLauriePeriodicTableRow(iAtom);
        int row2 = getHerschbachLauriePeriodicTableRow(jAtom);
        if (row1 > row2) {
            int temp = row2;
            row2 = row1;
            row1 = temp;
        }
        int index = Collections.binarySearch(mmffParameters.herschbachLaurie, new Float[]{(float) row1, (float) row2}, (Float[] p1, Float[] p2) -> {
            if (!p1[0].equals(p2[0])) {
                return Integer.compare(p1[0].intValue(), p2[0].intValue());
            } else {
                return Integer.compare(p1[1].intValue(), p2[1].intValue());
            }
        });
        if (index >= 0) {
            Float[] parameters = mmffParameters.herschbachLaurie.get(index);
            if (row1 == parameters[0] && row2 == parameters[1]) {
                return parameters;
            }
        }
        throw new NullPointerException(String.format("could not find Herschbach-Laurie parameters for %s-%s bond", iAtom.getAtomTypeName(), jAtom.getAtomTypeName()));
    }

    static void checkParametersAssigned(IAtom iAtom) {
        if (iAtom.getProperty(MMFF94_TYPE) == null) {
            throw new RuntimeException("MMFF94 parameters need to be assigned first");
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
            checkParametersAssigned(atomContainer.getAtom(0));
        }
        double energy = 0;
        energy += bondStretchingComponent.calculateEnergy(atomContainer);
        energy += angleBendingComponent.calculateEnergy(atomContainer);
        energy += stretchBendComponent.calculateEnergy(atomContainer);
        energy += outOfPlaneComponent.calculateEnergy(atomContainer);
        energy += torsionComponent.calculateEnergy(atomContainer);
        energy += vdwComponent.calculateEnergy(atomContainer);
        energy += electrostaticComponent.calculateEnergy(atomContainer);
        return energy;
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
                    if (parameters[2] != j || parameters[3] > k) {
                        break;
                    }
                    if (parameters[0] == torsionType && parameters[3] == k) {
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
        throw new RuntimeException("no empirical torsion paramters were found for " + iAtom.getAtomTypeName());
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

    /**
     *
     * @param atom1
     * @param atom2
     * @param n
     * @return array with the first element describing if it was separated by N
     * bonds and the second one if it was more than N.
     */
    private static boolean[] isSeparatedByNOrMoreBonds(IAtom atom1, IAtom atom2, int n) {
        if (n < 0) {
            throw new IllegalArgumentException();
        } else if (n == 0) {
            return new boolean[]{atom1.equals(atom2), !atom1.equals(atom2)};
        }

        HashSet<IAtom> searched = new HashSet<>();
        Queue<IAtom> queue = new LinkedList<>();
        queue.add(atom1);
        searched.add(atom1);

        int depth = 0;

        while (!queue.isEmpty() && depth < n) {
            int levelSize = queue.size();
            depth++;

            for (int i = 0; i < levelSize; i++) {
                IAtom currentAtom = queue.poll();

                for (IBond bond : currentAtom.bonds()) {
                    IAtom neighbor = bond.getOther(currentAtom);
                    if (neighbor.equals(atom2)) {
                        return new boolean[]{depth == n, false};
                    }
                    if (searched.add(neighbor)) {
                        queue.add(neighbor);
                    }
                }
            }
        }

        return new boolean[]{false, true}; // If we exit the loop, atom2 was not found within n bonds.
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
