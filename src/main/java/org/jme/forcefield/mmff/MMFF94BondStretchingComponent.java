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
        for (IBond bond : atomContainer.bonds()) {
            double energy = calculateEnergy(bond, bond.getProperty(MMFF94_PARAMETER_STRETCH));
            totalEnergy += energy;
            if (LOGGER.isLoggable(Level.FINE)) {
                LOGGER.fine(String.format("Stretch:\t%s-%s\t%f", bond.getBegin().getAtomTypeName(), bond.getEnd().getAtomTypeName(), energy));
            }
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
        return calculateEnergy(bond, parameters);
    }

    private double calculateEnergy(IBond bond, MMFF94Parameters.StretchParameters parameters) {
        Objects.requireNonNull(forceField);
        IAtom iAtom = bond.getBegin();
        IAtom jAtom = bond.getEnd();
        double deltaR = iAtom.getPoint3d().distance(jAtom.getPoint3d()) - parameters.r0;
        double energy = .5 * 143.9325 * parameters.kb * deltaR * deltaR * (1 - 2 * deltaR + (7d / 12d) * 4 * deltaR * deltaR);
        return ForceField.EnergyUnit.KCAL_PER_MOL.convertTo(energy, forceField.getEnergyUnit());
    }

}
