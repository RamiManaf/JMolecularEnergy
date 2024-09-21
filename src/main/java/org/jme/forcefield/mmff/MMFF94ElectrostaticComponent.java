/*
 * The MIT License
 *
 * Copyright 2024 Rami Manaf Abdullah.
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

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import static org.jme.forcefield.mmff.MMFF94.checkParametersAssigned;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.jme.forcefield.EnergyComponent;
import org.jme.forcefield.ForceField;
import static org.jme.forcefield.mmff.MMFF94.MMFF94_TYPE;

/**
 * The {@code MMFF94ElectrostaticComponent} class represents the electrostatic
 * interaction component in the MMFF94 force field. This component is
 * responsible for calculating the electrostatic energy between atoms with
 * partial and formal charges on atoms within a molecule.
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94ElectrostaticComponent extends EnergyComponent {

    private static final Logger LOGGER = Logger.getLogger(MMFF94ElectrostaticComponent.class.getName());

    /**
     * Calculates the electrostatic interaction between all the atoms in the
     * provided atom container if the atoms are separated by three or more bonds
     * or they are not in the same molecule. The dielectric constant is
     * considered 1.
     *
     * @param atomContainer
     * @return
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        Objects.requireNonNull(forceField);
        return calculateEnergy(atomContainer, forceField.getDielectricConstant());
    }

    /**
     * calculates the electrostatic interaction between all the atoms in the
     * provided atom container if the atoms are separated by three or more or
     * they are not in the same molecule.
     *
     * @param atomContainer
     * @param dielectricConstant
     * @return
     */
    public double calculateEnergy(IAtomContainer atomContainer, double dielectricConstant) {
        HashMap<List<IAtom>, Boolean> nonBondedInteraction = (HashMap<List<IAtom>, Boolean>) atomContainer.getProperty(MMFF94.MMFF94_NON_BONDED_INTERACTION);
        Objects.requireNonNull(nonBondedInteraction, "MMFF94 parameters need to be assigned first");
        double totalEnergy = 0;
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("    E L E C T R O S T A T I C     I N T E R A C T I O N      \n"
                    + "\n"
                    + " ------ATOMNAMES------   ATOM TYPES   I        J\n"
                    + "   I          J            I    J   CHARGE   CHARGE   DISTANCE    \n"
                    + " ----------------------------------------------------------------------------------------");
        }
        for (Map.Entry<List<IAtom>, Boolean> entry : nonBondedInteraction.entrySet()) {
            IAtom iAtom = entry.getKey().get(0);
            IAtom jAtom = entry.getKey().get(1);
            double energy = calculateEnergy(iAtom, jAtom, dielectricConstant, entry.getValue());
            totalEnergy += energy;
        }
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.fine(String.format("Total Electrostatic Interaction Energy =\t%.5f\n", totalEnergy));
        }
        return totalEnergy;
    }

    /**
     * Calculates electrostatic interaction energy between the two atoms.
     * Generally intermolecular interactions are calculated between non-bonded
     * atoms or atoms separated by three bonds (two atoms at least)
     *
     * @param iAtom
     * @param jAtom
     * @param dielectricConstant dielectric constant of the solvent. Default
     * value is 1
     * @return
     */
    public double calculateEnergy(IAtom iAtom, IAtom jAtom, double dielectricConstant) {
        checkParametersAssigned(iAtom);
        checkParametersAssigned(jAtom);
        if (iAtom.getContainer() != null && iAtom.getContainer().equals(jAtom.getContainer())) {
            HashMap<List<IAtom>, Boolean> nonBondedInteraction = (HashMap<List<IAtom>, Boolean>) iAtom.getContainer().getProperty(MMFF94.MMFF94_NON_BONDED_INTERACTION);
            List<IAtom> atoms;
            if (iAtom.getIndex() < jAtom.getIndex()) {
                atoms = Arrays.asList(iAtom, jAtom);
            } else {
                atoms = Arrays.asList(jAtom, iAtom);
            }
            if (LOGGER.isLoggable(Level.FINER)) {
                LOGGER.finer("    E L E C T R O S T A T I C     I N T E R A C T I O N      \n"
                        + "\n"
                        + " ------ATOMNAMES------   ATOM TYPES   I        J\n"
                        + "   I          J            I    J   CHARGE   CHARGE   DISTANCE    \n"
                        + " ----------------------------------------------------------------------------------------");
            }
            return calculateEnergy(iAtom, jAtom, dielectricConstant, nonBondedInteraction.get(atoms));
        } else {
            return calculateEnergy(iAtom, jAtom, dielectricConstant, false);
        }
    }

    private double calculateEnergy(IAtom iAtom, IAtom jAtom, double dielectricConstant, boolean separatedBy3Bonds) {
        Objects.requireNonNull(forceField);
        double qI = iAtom.getCharge();
        double qJ = jAtom.getCharge();
        double distance = iAtom.getPoint3d().distance(jAtom.getPoint3d());
        if (forceField.getCutoffDistance() > 0 && distance > forceField.getCutoffDistance()) {
            return 0;
        }
        double electrostaticEnergy = 332.0716 * qI * qJ / (dielectricConstant * (distance + 0.05));
        if (iAtom.getContainer() == jAtom.getContainer() && separatedBy3Bonds) {
            electrostaticEnergy *= 0.75;
        }
        electrostaticEnergy = ForceField.EnergyUnit.KCAL_PER_MOL.convertTo(electrostaticEnergy, forceField.getEnergyUnit());
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer(String.format("%s #%d\t%s #%d\t%d\t%d\t%.3f\t%.3f\t%.3f", iAtom.getSymbol(), iAtom.getIndex() + 1, jAtom.getSymbol(), jAtom.getIndex() + 1, iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), iAtom.getCharge(), jAtom.getCharge(), distance));
        }
        return electrostaticEnergy;
    }

}
