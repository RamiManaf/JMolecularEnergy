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
package org.jme.forcefield;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * The {@code ForceField} class serves as an abstract base class for defining
 * different types of force fields. A force field represents a mathematical
 * model that describes the potential energy of a molecular system as a function
 * of atomic positions.
 *
 * <p>
 * <b>Key responsibilities:</b></p>
 * <ul>
 * <li>Defining the potential energy function for molecular systems.</li>
 * <li>Atom typing and providing parameters for various energy terms, such as
 * bond stretching, angle bending, and nonbonded interactions.</li>
 * <li>Managing the energy components of the force field.</li>
 * </ul>
 *
 * @author Rami Manaf Abdullah
 */
public abstract class ForceField {

    private double dielectricConstant = 1;
    private double cutoffDistance = 0;
    protected EnergyUnit energyUnit = EnergyUnit.KCAL_PER_MOL;
    protected List<EnergyComponent> components = new ArrayList<>();

    /**
     * Assigns the force field parameters to the given molecule.
     * <p>
     * This method must be called before any energy calculations can be
     * performed. It sets up the necessary parameters specific to the molecule.
     * If the molecule undergoes structural changes, such as the formation or
     * breakup of bonds, this method needs to be called again to update the
     * parameters accordingly.
     * </p>
     *
     * @param atomContainer The container holding the atoms of the molecule.
     */
    public abstract void assignParameters(IAtomContainer atomContainer);

    /**
     * Calculates the potential energy for the given molecule.
     *
     * @param atomContainer The container holding the atoms of the molecule.
     * @return The potential energy.
     */
    public abstract double calculateEnergy(IAtomContainer atomContainer);

    /**
     * Retrieves the list of individual energy components associated with this
     * force field.
     *
     * @return An unmodifiable {@code List} of {@code EnergyComponent} objects
     * associated with this force field.
     */
    public List<EnergyComponent> getEnergyComponents() {
        return Collections.unmodifiableList(components);
    }

    /**
     * Adds a new energy component to the force field.
     *
     * @param component The {@code EnergyComponent} to be added.
     */
    public void addEnergyComponent(EnergyComponent component) {
        Objects.requireNonNull(component);
        component.setForceField(this);
        components.add(component);
    }

    /**
     * Removes an existing energy component from the force field.
     *
     * @param component The {@code EnergyComponent} to be removed.
     */
    public void removeEnergyComponent(EnergyComponent component) {
        Objects.requireNonNull(component);
        component.setForceField(null);
        components.remove(component);
    }

    /**
     * Sets the unit of energy to be used by this force field. The default
     * energy unit is kcal/mol.
     *
     * @param energyUnit The {@code EnergyUnit} to set.
     */
    public void setEnergyUnit(EnergyUnit energyUnit) {
        this.energyUnit = Objects.requireNonNull(energyUnit);
    }

    /**
     * Retrieves the current unit of energy used by this force field. The
     * default energy unit is kcal/mol.
     *
     * @return The {@code EnergyUnit} currently set for this force field.
     */
    public EnergyUnit getEnergyUnit() {
        return energyUnit;
    }

    /**
     * Enum representing common energy units used in force field calculations.
     */
    public enum EnergyUnit {
        KCAL_PER_MOL(4184.0),
        KJ_PER_MOL(1000.0),
        JOULE_PER_MOL(1.0);

        private final double toJouleFactor;

        private EnergyUnit(double toJouleFactor) {
            this.toJouleFactor = toJouleFactor;
        }

        /**
         * Converts a value from the current energy unit to the target energy
         * unit.
         *
         * @param value The energy value in the current unit.
         * @param toUnit The target energy unit to convert to.
         * @return The converted value in the target unit.
         */
        public double convertTo(double value, EnergyUnit toUnit) {
            if (this == toUnit) {
                return value;
            }
            return (value * this.toJouleFactor) / toUnit.toJouleFactor;
        }

    }

    /**
     * Sets the non-bonded cutoff distance for the force field.
     *
     * <p>
     * The non-bonded cutoff distance defines the maximum distance at which
     * non-bonded interactions (e.g., van der Waals and electrostatic forces)
     * are considered. Interactions beyond this distance will not be calculated,
     * improving computational efficiency. A value of 0 indicates an unlimited
     * cutoff distance, meaning that all non-bonded interactions will be
     * considered without restriction.</p>
     *
     * @param cutoffDistance in angstroms
     */
    public void setCutoffDistance(double cutoffDistance) {
        this.cutoffDistance = cutoffDistance;
    }

    /**
     * Retrieves the current non-bonded cutoff distance in angstroms.
     *
     * @return
     */
    public double getCutoffDistance() {
        return cutoffDistance;
    }

    /**
     * Sets the dielectric constant for the force field.
     *
     * <p>
     * The dielectric constant is used to model the effect of the environment on
     * electrostatic interactions within the force field. A higher dielectric
     * constant typically represents a more polar solvent, while a lower value
     * approximates vacuum conditions.</p>
     *
     * @param dielectricConstant
     */
    public void setDielectricConstant(double dielectricConstant) {
        this.dielectricConstant = dielectricConstant;
    }

    /**
     * Retrieves the current dielectric constant.
     *
     * @return
     */
    public double getDielectricConstant() {
        return dielectricConstant;
    }

}
