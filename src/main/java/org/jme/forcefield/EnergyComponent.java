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
package org.jme.forcefield;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Represents a single energy component within a force field.
 * <p>
 * Each {@code EnergyComponent} contributes a distinct portion of the total
 * potential energy for a molecule. Extending this class encapsulate the
 * calculations relevant to their specific type of energy, such as bond
 * stretching, angle bending, or torsional forces. After instantiation, the
 * energy component must be added to a force field to be functional.
 * </p>
 *
 * @author Rami Manaf Abdullah
 */
public abstract class EnergyComponent {

    protected ForceField forceField;

    void setForceField(ForceField forceField) {
        this.forceField = forceField;
    }

    /**
     * Calculates the potential energy with respect to atomic positions for this
     * component.
     *
     * @param atomContainer The container holding the molecule's atoms.
     * @return The potential energy as a {@code double}.
     */
    public abstract double calculateEnergy(IAtomContainer atomContainer);

}
