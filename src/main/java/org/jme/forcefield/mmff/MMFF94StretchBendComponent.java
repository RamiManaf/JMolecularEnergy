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
import static org.jme.forcefield.mmff.MMFF94.MMFF94_PARAMETER_STRETCH;
import static org.jme.forcefield.mmff.MMFF94.MMFF94_TYPE;
import static org.jme.forcefield.mmff.MMFF94.checkParametersAssigned;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.jme.forcefield.EnergyComponent;
import org.jme.forcefield.ForceField;

/**
 * The {@code MMFF94StretchBendComponent} class represents the stretch-bend
 * interaction component in the MMFF94 force field. This component is
 * responsible for calculating the energy associated with coupled bond
 * stretching and angle bending deformations in a molecular structure.
 * <p>
 * Stretch-bend interactions arise when the length of a bond and the angle
 * between two bonds involving a central atom change simultaneously.
 * </p>
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94StretchBendComponent extends EnergyComponent {

    private static final Logger LOGGER = Logger.getLogger(MMFF94StretchBendComponent.class.getName());

    /**
     * Calculates the energy term that couple the stretching and bending
     * energies of a three consecutively bonded atoms.
     *
     * @param atomContainer
     * @return
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        Map<List<IAtom>, Float[]> angleBending = ((Map<List<IAtom>, Float[]>) atomContainer.getProperty(MMFF94.MMFF94_ANGLE_BENDING));
        Map<List<IAtom>, Float[]> stretchBend = ((Map<List<IAtom>, Float[]>) atomContainer.getProperty(MMFF94.MMFF94_STRETCH_BEND));
        Objects.requireNonNull(stretchBend, "MMFF94 parameters need to be assigned first");
        double totalEnergy = 0;
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("      S T R E T C H   B E N D I N G         \n"
                    + "\n"
                    + " -------ATOMS-------  -ATOM TYPES-    FF     VALENCE      DELTA      DELTA     STRAIN      F CON \n"
                    + "  I       J       K     I    J    K  CLASS    ANGLE       ANGLE      R(I,J)    ENERGY       I-J\n"
                    + " -------------------------------------------------------------------------------------------------");
        }
        for (Map.Entry<List<IAtom>, Float[]> entry : stretchBend.entrySet()) {
            List<IAtom> atoms = entry.getKey();
            double energy = calculateEnergy(atoms.get(0), atoms.get(1), atoms.get(2), angleBending.get(atoms), entry.getValue());
            totalEnergy += energy;
        };
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.fine(String.format("Total Stretch-Bend Strain Energy =\t%.5f\n", totalEnergy));
        }
        return totalEnergy;
    }

    /**
     * Calculates stretch bend energy that couples i-j and j-k bonds stretching
     * to the angle formed by i-j-k. This energy favors bonds stretching when
     * the angle decrease.
     *
     * @param iAtom first atom
     * @param jAtom middle atom
     * @param kAtom end atom
     * @return
     */
    public double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom) {
        if (iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom) == null) {
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k");
        }
        checkParametersAssigned(iAtom);
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("      S T R E T C H   B E N D I N G         \n"
                    + "\n"
                    + " -------ATOMS-------  -ATOM TYPES-    FF     VALENCE      DELTA      DELTA     STRAIN      F CON \n"
                    + "  I       J       K     I    J    K  CLASS    ANGLE       ANGLE      R(I,J)    ENERGY       I-J\n"
                    + " -------------------------------------------------------------------------------------------------");
        }
        Float[] angleBendingParameters = ((Map<List<IAtom>, Float[]>) iAtom.getContainer().getProperty(MMFF94.MMFF94_ANGLE_BENDING)).get(Arrays.asList(iAtom, jAtom, kAtom));
        Float[] stretchBendParameters = ((Map<List<IAtom>, Float[]>) iAtom.getContainer().getProperty(MMFF94.MMFF94_STRETCH_BEND)).get(Arrays.asList(iAtom, jAtom, kAtom));
        return calculateEnergy(iAtom, jAtom, kAtom, angleBendingParameters, stretchBendParameters);
    }

    public double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, Float[] angleBendingParameters, Float[] stretchBendParameters) {
        Objects.requireNonNull(forceField);
        if (jAtom.<MMFF94Parameters.GeometricParameters>getProperty(MMFF94_PARAMETER_GEOMETRIC_PROPERTIES).ideallyLinear) {
            return 0;
        }
        double angle = (GeometryUtils.calculateAngle(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d()) * 180 / Math.PI);
        double deltaAngle = angle - angleBendingParameters[5];
        MMFF94Parameters.StretchParameters ijStretch = iAtom.getContainer().getBond(iAtom, jAtom).<MMFF94Parameters.StretchParameters>getProperty(MMFF94_PARAMETER_STRETCH);
        MMFF94Parameters.StretchParameters jkStretch = iAtom.getContainer().getBond(jAtom, kAtom).<MMFF94Parameters.StretchParameters>getProperty(MMFF94_PARAMETER_STRETCH);
        double deltaRij = iAtom.getPoint3d().distance(jAtom.getPoint3d()) - ijStretch.r0;
        double deltaRjk = jAtom.getPoint3d().distance(kAtom.getPoint3d()) - jkStretch.r0;
        double energy = 2.51210 * (stretchBendParameters[4] * deltaRij + stretchBendParameters[5] * deltaRjk) * deltaAngle;
        energy = ForceField.EnergyUnit.KCAL_PER_MOL.convertTo(energy, forceField.getEnergyUnit());
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer(String.format("%s\t%s #%d\t%s\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t", iAtom.getSymbol(), jAtom.getSymbol(), jAtom.getIndex() + 1, kAtom.getSymbol(), iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), kAtom.getProperty(MMFF94_TYPE), stretchBendParameters[0].intValue(), angle, deltaAngle, deltaRij, energy, stretchBendParameters[4]));
        }
        return energy;
    }

}
