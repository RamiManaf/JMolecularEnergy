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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.logging.Level;
import java.util.logging.Logger;
import static org.jme.forcefield.mmff.MMFF94.MMFF94_TYPE;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.jme.forcefield.EnergyComponent;
import org.jme.forcefield.ForceField;

/**
 * The {@code MMFF94VdwComponent} class represents the van der Waals interaction
 * component in the MMFF94 force field. This class is responsible for
 * calculating the non-bonded interactions between atoms due to van der Waals
 * forces, which arise from transient dipole-induced dipole interactions.
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94VdwComponent extends EnergyComponent {

    private static final Logger LOGGER = Logger.getLogger(MMFF94VdwComponent.class.getName());

    /**
     * Calculates the Van Der Waals interaction energy between all the atoms in
     * the provided atom container if the atoms are separated by three or more
     * bonds or they are not from the same molecule.
     *
     * @param atomContainer
     * @return
     */
    @Override
    public double calculateEnergy(IAtomContainer atomContainer) {
        HashMap<List<IAtom>, Boolean> nonBondedInteraction = (HashMap<List<IAtom>, Boolean>) atomContainer.getProperty(MMFF94.MMFF94_NON_BONDED_INTERACTION);
        Objects.requireNonNull(nonBondedInteraction, "MMFF94 parameters need to be assigned first");
        double totalEnergy = 0;
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("    V D W   E N E R G Y      \n"
                    + "\n"
                    + " ------ATOMNAMES------   ATOM TYPES   \n"
                    + "   I          J            I    J     DISTANCE   R*   EPS\n"
                    + " ----------------------------------------------------------------------------------------");
        }
        for (Map.Entry<List<IAtom>, Boolean> entry : nonBondedInteraction.entrySet()) {
            IAtom iAtom = entry.getKey().get(0);
            IAtom jAtom = entry.getKey().get(1);
            double energy = calculateVdwEnergy(iAtom, jAtom);
            totalEnergy += energy;
        };
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.fine(String.format("Total Vdw Energy =\t%.5f\n", totalEnergy));
        }
        return totalEnergy;
    }

    /**
     * Calculate Van Der Waals energy between two atoms. Generally
     * intermolecular interactions are calculated between non-bonded atoms or
     * atoms separated by three bonds (two atoms at least)
     *
     * @param iAtom
     * @param jAtom
     * @return
     */
    public double calculateEnergy(IAtom iAtom, IAtom jAtom) {
        MMFF94.checkParametersAssigned(iAtom);
        MMFF94.checkParametersAssigned(jAtom);
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer("    V D W   E N E R G Y      \n"
                    + "\n"
                    + " ------ATOMNAMES------   ATOM TYPES   \n"
                    + "   I          J            I    J     DISTANCE   R*   EPS\n"
                    + " ----------------------------------------------------------------------------------------");
        }
        return calculateVdwEnergy(iAtom, jAtom);
    }

    private double calculateVdwEnergy(IAtom iAtom, IAtom jAtom) {
        Objects.requireNonNull(forceField);
        int i = iAtom.getProperty(MMFF94_TYPE);
        int j = jAtom.getProperty(MMFF94_TYPE);
        //there is a 4 number gap in mmff atom type number (83-86)
        MMFF94Parameters.VdwParameters iParameters = iAtom.getProperty(MMFF94.MMFF94_VDW_INTERACTION);
        MMFF94Parameters.VdwParameters jParameters = jAtom.getProperty(MMFF94.MMFF94_VDW_INTERACTION);
        double RStarIJ;
        if (iAtom.getAtomTypeName().equals(jAtom.getAtomTypeName())) {
            //like pairs
            RStarIJ = iParameters.A * Math.pow(iParameters.alpha, 0.25d);
        } else {
            //unlike pairs
            double RStarII = iParameters.A * Math.pow(iParameters.alpha, 0.25d);
            double RStarJJ = jParameters.A * Math.pow(jParameters.alpha, 0.25d);

            double B = (iParameters.DA == MMFF94Parameters.VdwParameters.DONER || jParameters.DA == MMFF94Parameters.VdwParameters.DONER) ? 0 : 0.2;
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
        if (forceField.getCutoffDistance() > 0 && distance > forceField.getCutoffDistance()) {
            return 0;
        }
        double vdwEnergy = epsilonIJ * Math.pow((1.07 * RStarIJ) / (distance + 0.07 * RStarIJ), 7);
        vdwEnergy = vdwEnergy * (((1.12 * Math.pow(RStarIJ, 7)) / (Math.pow(distance, 7) + 0.12 * Math.pow(RStarIJ, 7))) - 2);
        vdwEnergy = ForceField.EnergyUnit.KCAL_PER_MOL.convertTo(vdwEnergy, forceField.getEnergyUnit());
        if (LOGGER.isLoggable(Level.FINER)) {
            LOGGER.finer(String.format("%s #%d\t%s #%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t", iAtom.getSymbol(), iAtom.getIndex() + 1, jAtom.getSymbol(), jAtom.getIndex() + 1, iAtom.getProperty(MMFF94_TYPE), jAtom.getProperty(MMFF94_TYPE), distance, RStarIJ, epsilonIJ));
        }
        return vdwEnergy;
    }

}
