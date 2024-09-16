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

import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import static org.jme.forcefield.mmff.MMFF94.MMFF94_PARAMETER_STRETCH;
import static org.jme.forcefield.mmff.MMFF94.checkParametersAssigned;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.jme.forcefield.EnergyComponent;
import org.jme.forcefield.ForceField;
import static org.jme.forcefield.mmff.MMFF94.MMFF94_TYPE;

/**
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94BondStretchingComponent extends EnergyComponent {

    private static final Logger LOGGER = Logger.getLogger(MMFF94BondStretchingComponent.class.getName());

    /**
     * calculate the bond stretching energy for all the bonds in the atom
     * container
     *
     * @param atomContainer
     * @return
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        double totalEnergy = 0;
        if (!atomContainer.isEmpty()) {
            checkParametersAssigned(atomContainer.getAtom(0));
        }
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("        B O N D   S T R E T C H I N G      \n"
                    + "\n"
                    + " ------ATOMNAMES------   ATOM TYPES   FF     BOND     IDEAL             STRAIN     FORCE\n"
                    + "   I          J            I    J   CLASS   LENGTH   LENGTH    DIFF.    ENERGY   CONSTANT\n"
                    + " ----------------------------------------------------------------------------------------");
        }
        for (IBond bond : atomContainer.bonds()) {
            double energy = calculateEnergy(bond, bond.getProperty(MMFF94_PARAMETER_STRETCH));
            totalEnergy += energy;
        }
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.fine("Total Bond Stretching Energy =\t%.5f\n".formatted(totalEnergy));
        }
        return totalEnergy;
    }

    /**
     * calculate the bond stretching energy in the provided bond
     *
     * @param bond
     * @return
     */
    public double calculateEnergy(IBond bond) {
        checkParametersAssigned(bond.getBegin());
        MMFF94Parameters.StretchParameters parameters = bond.getProperty(MMFF94_PARAMETER_STRETCH);
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("        B O N D   S T R E T C H I N G      \n"
                    + "\n"
                    + " ------ATOMNAMES------   ATOM TYPES   FF     BOND     IDEAL             STRAIN     FORCE\n"
                    + "   I          J            I    J   CLASS   LENGTH   LENGTH    DIFF.    ENERGY   CONSTANT\n"
                    + " ----------------------------------------------------------------------------------------");
        }
        return calculateEnergy(bond, parameters);
    }

    private double calculateEnergy(IBond bond, MMFF94Parameters.StretchParameters parameters) {
        Objects.requireNonNull(forceField);
        IAtom iAtom = bond.getBegin();
        IAtom jAtom = bond.getEnd();
        double length = iAtom.getPoint3d().distance(jAtom.getPoint3d());
        double deltaR = length - parameters.r0;
        double energy = .5 * 143.9325 * parameters.kb * deltaR * deltaR * (1 - 2 * deltaR + (7d / 12d) * 4 * deltaR * deltaR);
        energy = ForceField.EnergyUnit.KCAL_PER_MOL.convertTo(energy, forceField.getEnergyUnit());
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer(String.format("%s #%d\t%s #%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", iAtom.getSymbol(), iAtom.getIndex(), jAtom.getSymbol(), jAtom.getIndex(), iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), parameters.bondType, length, parameters.r0, deltaR, energy, parameters.kb));
        }
        return energy;
    }

}
