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
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.jme.forcefield.GeometryUtils;
import static org.jme.forcefield.mmff.MMFF94.MMFF94_PARAMETER_GEOMETRIC_PROPERTIES;
import static org.jme.forcefield.mmff.MMFF94.MMFF94_TYPE;
import static org.jme.forcefield.mmff.MMFF94.checkParametersAssigned;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.jme.forcefield.EnergyComponent;
import org.jme.forcefield.ForceField;

/**
 *
 * The {@code MMFF94AngleBendingComponent} class represents the angle bending
 * component in the MMFF94 force field. This component is responsible for
 * calculating the energy associated with deviations from the ideal bond angles
 * between three connected atoms in a molecule.
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94AngleBendingComponent extends EnergyComponent {

    private static final Logger LOGGER = Logger.getLogger(MMFF94AngleBendingComponent.class.getName());

    /**
     * calculates the angle bending energy for all angles formed by consecutive
     * three atoms in the atom container
     *
     * @param atomContainer
     * @return energy
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        Map<List<IAtom>, Float[]> angleBending = ((Map<List<IAtom>, Float[]>) atomContainer.getProperty(MMFF94.MMFF94_ANGLE_BENDING));
        Objects.requireNonNull(angleBending, "MMFF94 parameters need to be assigned first");
        double totalEnergy = 0;
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("      A N G L E   B E N D I N G         \n"
                    + "\n"
                    + " -------ATOMS-------   -ATOM TYPES-   FF     VALENCE      IDEAL                 STRAIN     FORCE\n"
                    + "  I       J       K     I    J    K  CLASS    ANGLE       ANGLE      DIFF.      ENERGY   CONSTANT\n"
                    + " -------------------------------------------------------------------------------------------------");
        }
        for (Map.Entry<List<IAtom>, Float[]> entry : angleBending.entrySet()) {
            List<IAtom> atoms = entry.getKey();
            double energy = calculateEnergy(atoms.get(0), atoms.get(1), atoms.get(2), entry.getValue());
            totalEnergy += energy;
        }
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.fine(String.format("Total Angle Bending Energy =\t%.5f\n", totalEnergy));
        }
        return totalEnergy;
    }

    /**
     * calculates the angle bending energy for the angle that is formed by three
     * consecutive atoms.
     *
     * @param iAtom the atom at one side
     * @param jAtom the atom in the middle
     * @param kAtom the atom at the other side
     * @return energy
     */
    public double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom) {
        if (iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom) == null) {
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k");
        }
        checkParametersAssigned(iAtom);
        Float[] parameters = ((Map<List<IAtom>, Float[]>) iAtom.getContainer().getProperty(MMFF94.MMFF94_ANGLE_BENDING)).get(Arrays.asList(iAtom, jAtom, kAtom));
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("      A N G L E   B E N D I N G         \n"
                    + "\n"
                    + " -------ATOMS-------   -ATOM TYPES-   FF     VALENCE      IDEAL                 STRAIN     FORCE\n"
                    + "  I       J       K     I    J    K  CLASS    ANGLE       ANGLE      DIFF.      ENERGY   CONSTANT\n"
                    + " -------------------------------------------------------------------------------------------------");
        }
        return calculateEnergy(iAtom, jAtom, kAtom, parameters);
    }

    private double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, Float[] parameters) {
        Objects.requireNonNull(forceField);
        double angle = GeometryUtils.calculateAngle(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d()) * 180 / Math.PI;
        double deltaAngle = angle - parameters[5];
        double energy;
        if (jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear) {
            energy = 143.9325 * parameters[4] * (1 + Math.cos(Math.toRadians(angle)));
        } else {
            energy = Math.toRadians(Math.toRadians(143.9325)) * parameters[4] * .5 * deltaAngle * deltaAngle * (1 - 0.006981317 * deltaAngle);
        }
        energy = ForceField.EnergyUnit.KCAL_PER_MOL.convertTo(energy, forceField.getEnergyUnit());
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer(String.format("%s\t%s #%d\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f", iAtom.getSymbol(), jAtom.getSymbol(), jAtom.getIndex() + 1, kAtom.getSymbol(), iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), kAtom.getProperty(MMFF94_TYPE), parameters[0].intValue(), angle, parameters[5], deltaAngle, energy, parameters[4]));
        }
        return energy;
    }

}
