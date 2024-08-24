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
package org.jme.constraint;

import org.jme.forcefield.ForceField;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Restriction for molecules atom movement based on a specific criteria. This
 * class provide ready to use implementation of commonly used constrains such as
 * distance constrain between atoms.
 *
 * @author Rami Manaf Abdullah
 * @param <F>
 */
public interface Constraint<F extends ForceField> {

    /**
     * checks if the container meet the constraint. If the container violates
     * the constraint the method will return false. Otherwise the container
     * meets the constraint and the method will return true.
     *
     * @param forcefield forcefield
     * @param container
     * @return
     */
    public boolean check(F forcefield, IAtomContainer container);

    /**
     * combine the two constraints using AND logical operator in a new
     * constraint.
     *
     * @param constraint
     * @return
     */
    public default Constraint and(Constraint<F> constraint) {
        return (Constraint<F>) (F forcefield, IAtomContainer container) -> Constraint.this.check(forcefield, container) && constraint.check(forcefield, container);
    }

    /**
     * combine the two constraints using OR logical operator in a new
     * constraint.
     *
     * @param constraint
     * @return
     */
    public default Constraint or(Constraint<F> constraint) {
        return (Constraint<F>) (F forcefield, IAtomContainer container) -> Constraint.this.check(forcefield, container) || constraint.check(forcefield, container);
    }

    /**
     * invert this constraint using NOT logical operator and return a new
     * constraint.
     *
     * @return
     */
    public default Constraint not() {
        return (Constraint<F>) (F forcefield, IAtomContainer container) -> !Constraint.this.check(forcefield, container);
    }

}
