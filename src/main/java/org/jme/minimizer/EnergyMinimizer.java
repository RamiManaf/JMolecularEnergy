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
package org.jme.minimizer;

import java.util.function.BiConsumer;
import org.jme.constraint.Constraint;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * The {@code EnergyMinimizer} abstract class serves as a base class for
 * implementing various energy minimization algorithms aimed at optimizing the
 * geometry of molecular systems. This class provides a framework for minimizing
 * the potential energy surface to achieve stable molecular conformations.
 *
 * @author Rami Manaf Abdullah
 */
public abstract class EnergyMinimizer {

    protected int step;
    protected BiConsumer<Integer, Double> onStepConsumer;
    protected Constraint constraint;

    /**
     * A consumer that is called after each minimization step. It receives two
     * parameters, the step number and the energy after that step.
     *
     * @param onStepConsumer
     */
    public void setOnStepConsumer(BiConsumer<Integer, Double> onStepConsumer) {
        this.onStepConsumer = onStepConsumer;
    }

    /**
     * set a constraint for atoms movement while energy minimization. The
     * minimizer will prevent violation of constraints. However, in the case of
     * any further movement would case a violation, the minimization will stop.
     *
     * @param constraint
     */
    public void setConstraint(Constraint constraint) {
        this.constraint = constraint;
    }

    public Constraint getConstraint() {
        return constraint;
    }

    /**
     * Minimizes the potential energy of the specified molecular container
     * within the current thread.
     *
     * <p>
     * This method iteratively adjusts the molecular geometry to reduce the
     * potential energy, employing the specified step size and stopping criteria
     * based on the energy difference.</p>
     *
     * @param container the molecule to minimize
     * @param maximumSteps maximum number of steps
     * @param stepSize the step size
     * @param threshold the energy difference which the minimizer will stop if
     * last step minimization reduced energy less than it
     */
    public abstract void minimize(IAtomContainer container, int maximumSteps, double stepSize, double threshold);

    /**
     * Returns the step minimizer is currently at it
     *
     * @return
     */
    public int getStep() {
        return step;
    }

}
