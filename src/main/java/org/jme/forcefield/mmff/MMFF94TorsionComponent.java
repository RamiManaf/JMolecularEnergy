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

/**
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94TorsionComponent implements EnergyComponent {

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
        for (Map.Entry<List<IAtom>, Float[]> entry : torsion.entrySet()) {
            List<IAtom> atoms = entry.getKey();
            double energy = calculateEnergy(atoms.get(0), atoms.get(1), atoms.get(2), atoms.get(3), entry.getValue());
            if (LOGGER.isLoggable(Level.FINE)) {
                LOGGER.fine(String.format("Torsion:\t%d-%d-%d-%d\t%f", atoms.get(0).getProperty(MMFF94_TYPE), atoms.get(1).getProperty(MMFF94_TYPE), atoms.get(2).getProperty(MMFF94_TYPE), atoms.get(3).getProperty(MMFF94_TYPE), energy));
            }
            totalEnergy += energy;
        };
        LOGGER.fine("total torsion:\t" + totalEnergy);
        return totalEnergy;

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
    public double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
        if (iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom) == null && kAtom.getBond(lAtom) == null) {
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k-l");
        }
        checkParametersAssigned(iAtom);
        Float[] torsionParameters = ((Map<List<IAtom>, Float[]>) iAtom.getContainer().getProperty(MMFF94.MMFF94_TORSION)).get(Arrays.asList(iAtom, jAtom, kAtom, lAtom));
        return calculateEnergy(iAtom, jAtom, kAtom, lAtom, torsionParameters);
    }

    public double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom, Float[] torsionParameters) {
        if ((Integer) kAtom.getProperty(MMFF94_TYPE) < (Integer) jAtom.getProperty(MMFF94_TYPE)
                || (jAtom.getProperty(MMFF94_TYPE).equals((Integer) kAtom.getProperty(MMFF94_TYPE)) && (Integer) iAtom.getProperty(MMFF94_TYPE) > (Integer) lAtom.getProperty(MMFF94_TYPE))) {
            IAtom temp = jAtom;
            jAtom = kAtom;
            kAtom = temp;
            temp = iAtom;
            iAtom = lAtom;
            lAtom = temp;
        }
        double torsionAngle = GeometryUtils.calculateTorsionAngle(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d(), lAtom.getPoint3d());//in radians
        return 0.5 * (torsionParameters[5] * (1 + Math.cos(torsionAngle)) + torsionParameters[6] * (1 - Math.cos(2 * torsionAngle)) + torsionParameters[7] * (1 + Math.cos(3 * torsionAngle)));
    }

}
