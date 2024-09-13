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

import java.util.List;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * A general interface for defining a force field
 *
 * @author Rami Manaf Abdullah
 */
public interface ForceField {

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
    public void assignParameters(IAtomContainer atomContainer);

    /**
     * Calculates the potential energy for the given molecule.
     *
     * @param atomContainer The container holding the atoms of the molecule.
     * @return The potential energy.
     */
    public double calculateEnergy(IAtomContainer atomContainer);

    /**
     * Retrieves the individual energy components of the force field.
     *
     * @return A {@code List} of {@code EnergyComponent} objects.
     */
    public List<EnergyComponent> getEnergyComponents();

}
