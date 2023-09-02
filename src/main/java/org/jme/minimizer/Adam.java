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

import org.jme.forcefield.ForceField;
import javax.vecmath.Point3d;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * energy minimization based on Adam optimization algorithm. This is commonly
 * used for avoiding gradient descent explosion
 *
 * @author Rami Manaf Abdullah
 */
public class Adam extends EnergyMinimizer {

    private final ForceField forcefield;
    private final double beta1;
    private final double beta2;
    private final double epsilon;

    /**
     * creates a new instance of Adam optimizer
     *
     * @param forceField the force field used in optimization
     * @param beta1
     * @param beta2
     * @param epsilon
     */
    public Adam(ForceField forceField, double beta1, double beta2, double epsilon) {
        this.forcefield = forceField;
        this.beta1 = beta1;
        this.beta2 = beta2;
        this.epsilon = epsilon;
    }

    @Override
    public void minimize(IAtomContainer container, int maximumSteps, double stepSize, double threshold) {
        step = 0;
        double[][] atomContainerGradients = new double[container.getAtomCount()][3];
        AdamOptimizer[][] adams = new AdamOptimizer[container.getAtomCount()][3];
        for (int i = 0; i < atomContainerGradients.length; i++) {
            atomContainerGradients[i] = new double[]{stepSize, stepSize, stepSize};
            for (int j = 0; j < 3; j++) {
                adams[i][j] = new AdamOptimizer(stepSize, beta1, beta2, epsilon);
            }
        }
        boolean initializedGradients = false;
        double lastCalculatedEnergy = forcefield.calculateEnergy(container);
        double energyPerStep = 0;
        if (onStepConsumer != null) {
            onStepConsumer.accept(step, lastCalculatedEnergy);
        }
        while (step < maximumSteps) {
            if (step != 0 && Math.abs(energyPerStep - lastCalculatedEnergy) <= threshold) {
                break;
            }
            energyPerStep = lastCalculatedEnergy;
            for (IAtom atom : container.atoms()) {
                Point3d point = atom.getPoint3d();
                double[] atomGradients = atomContainerGradients[atom.getIndex()];
                for (int i = 0; i < atomGradients.length; i++) {
                    double oneAxisGradient = atomGradients[i];
                    if (oneAxisGradient == 0) {
                        continue;
                    }
                    double initialEnergy = forcefield.calculateEnergy(container);
                    double distance = adams[atom.getIndex()][i].optimize((i == 0 ? point.x : (i == 1 ? point.y : point.z)), oneAxisGradient) - (i == 0 ? point.x : (i == 1 ? point.y : point.z));
                    movePoint(point, distance, i);
                    if (constraint == null || constraint.check(forcefield, container) || !initializedGradients) {
                        double difference = forcefield.calculateEnergy(container) - initialEnergy;
                        atomGradients[i] = -difference / Math.abs(distance);
                        if (!initializedGradients) {
                            movePoint(point, -distance, i);
                        } else {
                            if (difference > 0) {
                                movePoint(point, -distance, i);
                            } else {
                                lastCalculatedEnergy += difference;
                            }
                        }
                    } else {
                        movePoint(point, -distance, i);
                    }
                }
            }
            if (!initializedGradients) {
                initializedGradients = true;
            } else {
                step++;
                if (onStepConsumer != null) {
                    onStepConsumer.accept(step, lastCalculatedEnergy);
                }
            }
        }
    }

    private void movePoint(Point3d point, double distance, int i) {
        if (i == 0) {
            point.x += distance;
        } else if (i == 1) {
            point.y += distance;
        } else {
            point.z += distance;
        }
    }

    private static class AdamOptimizer {

        private final double learningRate;
        private final double beta1;
        private final double beta2;
        private final double epsilon;

        private double m;
        private double v;
        private int t;

        public AdamOptimizer(double learningRate, double beta1, double beta2, double epsilon) {
            this.learningRate = learningRate;
            this.beta1 = beta1;
            this.beta2 = beta2;
            this.epsilon = epsilon;
            this.m = 0.0;
            this.v = 0.0;
            this.t = 1;
        }

        public double optimize(double x, double gradient) {
            t++;
            m = beta1 * m + (1 - beta1) * gradient;
            v = beta2 * v + (1 - beta2) * gradient * gradient;
            double m_hat = m / (1 - Math.pow(beta1, t));
            double v_hat = v / (1 - Math.pow(beta2, t));
            x = x - learningRate * m_hat / (Math.sqrt(v_hat) + epsilon);
            return x;
        }
    }

}
