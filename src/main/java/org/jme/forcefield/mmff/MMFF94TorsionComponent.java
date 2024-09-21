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
import static org.jme.forcefield.mmff.MMFF94.MMFF94_TYPE;
import static org.jme.forcefield.mmff.MMFF94.checkParametersAssigned;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.jme.forcefield.EnergyComponent;
import org.jme.forcefield.ForceField;

/**
 * The {@code MMFF94TorsionComponent} class represents the torsional interaction
 * component in the MMFF94 force field. This class is responsible for
 * calculating the energy associated with the rotation around a bond between two
 * atoms, which influences the molecular conformation.
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94TorsionComponent extends EnergyComponent {

    private static final Logger LOGGER = Logger.getLogger(MMFF94TorsionComponent.class.getName());

    /**
     * Calculates the torsional energy associated with the rotation around a
     * bond, as defined by the dihedral angle formed by four consecutively
     * bonded atoms in a Newman projection.
     *
     * @param atomContainer
     * @return
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        Map<List<IAtom>, Float[]> torsion = ((Map<List<IAtom>, Float[]>) atomContainer.getProperty(MMFF94.MMFF94_TORSION));
        Objects.requireNonNull(torsion, "MMFF94 parameters need to be assigned first");
        double totalEnergy = 0;
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("    T O R S I O N A L    \n"
                    + "\n"
                    + " --------------ATOMS--------------  ---ATOM TYPES---   FF     TORSION    STERIC    --FORCE CONSTANTS--\n"
                    + "  I       J          K       L        I   J   K   L   CLASS    ANGLE     ENERGY     V1      V2      V3\n"
                    + " ------------------------------------------------------------------------------------------------------");
        }
        for (Map.Entry<List<IAtom>, Float[]> entry : torsion.entrySet()) {
            List<IAtom> atoms = entry.getKey();
            double energy = calculateEnergy(atoms.get(0), atoms.get(1), atoms.get(2), atoms.get(3), entry.getValue());
            totalEnergy += energy;
        };
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.fine(String.format("Total Torsion Strain Energy: =\t%.5f\n", totalEnergy));
        }
        return totalEnergy;

    }

    /**
     * Calculates the torsion energy for the torsion angle formed by i-j-k-l
     * atoms. The calculation is base on the angle formed by the rotation around
     * the j-k bond.
     *
     * @param iAtom
     * @param jAtom
     * @param kAtom
     * @param lAtom
     * @return
     */
    public double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
        if (iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom) == null && kAtom.getBond(lAtom) == null) {
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k-l");
        }
        checkParametersAssigned(iAtom);
        Float[] torsionParameters = ((Map<List<IAtom>, Float[]>) iAtom.getContainer().getProperty(MMFF94.MMFF94_TORSION)).get(Arrays.asList(iAtom, jAtom, kAtom, lAtom));
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("    T O R S I O N A L    \n"
                    + "\n"
                    + " --------------ATOMS--------------  ---ATOM TYPES---   FF     TORSION    STERIC    --FORCE CONSTANTS--\n"
                    + "  I       J          K       L        I   J   K   L   CLASS    ANGLE     ENERGY     V1      V2      V3\n"
                    + " ------------------------------------------------------------------------------------------------------");
        }
        return calculateEnergy(iAtom, jAtom, kAtom, lAtom, torsionParameters);
    }

    private double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom, Float[] torsionParameters) {
        Objects.requireNonNull(forceField);
        int i = iAtom.getProperty(MMFF94_TYPE);
        int j = jAtom.getProperty(MMFF94_TYPE);
        int k = kAtom.getProperty(MMFF94_TYPE);
        int l = lAtom.getProperty(MMFF94_TYPE);
        if (k < j || (j==k && i > l) || (j==k && i==l && jAtom.getIndex() > kAtom.getIndex())) {
            IAtom temp = jAtom;
            jAtom = kAtom;
            kAtom = temp;
            temp = iAtom;
            iAtom = lAtom;
            lAtom = temp;
        }
        double torsionAngle = GeometryUtils.calculateTorsionAngle(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d(), lAtom.getPoint3d());//in radians
        double energy = 0.5 * (torsionParameters[5] * (1 + Math.cos(torsionAngle)) + torsionParameters[6] * (1 - Math.cos(2 * torsionAngle)) + torsionParameters[7] * (1 + Math.cos(3 * torsionAngle)));
        energy = ForceField.EnergyUnit.KCAL_PER_MOL.convertTo(energy, forceField.getEnergyUnit());
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer(String.format("%s\t%s #%d\t%s #%d\t%s\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f", iAtom.getSymbol(), jAtom.getSymbol(), jAtom.getIndex() + 1, kAtom.getSymbol(), kAtom.getIndex() + 1, lAtom.getSymbol(), iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), kAtom.getProperty(MMFF94_TYPE), lAtom.getProperty(MMFF94_TYPE), torsionParameters[0].intValue(), Math.toDegrees(torsionAngle), energy, torsionParameters[5], torsionParameters[6], torsionParameters[7]));
        }
        return energy;
    }

}
