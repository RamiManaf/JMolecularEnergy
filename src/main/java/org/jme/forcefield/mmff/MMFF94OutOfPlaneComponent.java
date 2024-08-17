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
public class MMFF94OutOfPlaneComponent implements EnergyComponent {

    private static final Logger LOGGER = Logger.getLogger(MMFF94OutOfPlaneComponent.class.getName());

    /**
     * calculates the out of plane energy component in all the atoms that are
     * connected by four atoms.
     *
     * @param atomContainer
     * @return
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        Map<List<IAtom>, Float[]> outOfPlane = ((Map<List<IAtom>, Float[]>) atomContainer.getProperty(MMFF94.MMFF94_OUT_OF_PLANE));
        Objects.requireNonNull(outOfPlane, "MMFF94 parameters need to be assigned first");
        double totalEnergy = outOfPlane.entrySet().stream().mapToDouble((entry) -> {
            List<IAtom> atoms = entry.getKey();
            double energy = calculateEnergy(atoms.get(0), atoms.get(1), atoms.get(2), atoms.get(3), entry.getValue());
            if (LOGGER.isLoggable(Level.FINE)) {
                LOGGER.fine(String.format("OOP:\t%d-%d-%d-%d\t%.3f", atoms.get(0).getProperty(MMFF94_TYPE), atoms.get(1).getProperty(MMFF94_TYPE), atoms.get(2).getProperty(MMFF94_TYPE), atoms.get(3).getProperty(MMFF94_TYPE), energy));
            }
            return energy;
        }).sum();
        LOGGER.fine("total OOP:\t" + totalEnergy);
        return totalEnergy;
    }

    /**
     * calculates out of plane energy for the atoms
     *
     * @param iAtom first atom forming the plane
     * @param jAtom the middle atom forming the plane
     * @param kAtom the last atom forming the plane
     * @param lAtom out of the plane atom
     * @return
     */
    public double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom) {
        if (jAtom.getBond(iAtom) == null || jAtom.getBond(kAtom) == null && jAtom.getBond(lAtom) == null) {
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k/l where l is bonded to j");
        }
        checkParametersAssigned(iAtom);
        Float[] parameters = ((HashMap<List<IAtom>, Float[]>) iAtom.getContainer().getProperty(MMFF94.MMFF94_OUT_OF_PLANE)).get(Arrays.asList(iAtom, jAtom, kAtom, lAtom));
        return calculateEnergy(iAtom, jAtom, kAtom, lAtom, parameters);
    }

    private double calculateEnergy(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom, Float[] parameters) {
        if (parameters == null) {
            return 0;
        }
        double angle = GeometryUtils.calculateOutOfPlaneAngle(iAtom.getPoint3d(), jAtom.getPoint3d(), kAtom.getPoint3d(), lAtom.getPoint3d()) * 180 / Math.PI;
        return Math.toRadians(Math.toRadians(143.9325)) * .5 * parameters[4] * angle * angle;
    }
}