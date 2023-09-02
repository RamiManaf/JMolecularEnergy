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
 * Gradient descent energy minimizer
 *
 * @author Rami Manaf Abdullah
 */
public class GradientDescent extends EnergyMinimizer {

    private ForceField forcefield;
    private int maximumSteps;
    private float stepSize;

    public GradientDescent(ForceField forcefield) {
        this.forcefield = forcefield;
    }

    @Override
    public void minimize(IAtomContainer container, int maximumSteps, double stepSize, double threshold) {
        step = 0;
        double[][] atomContainerGradients = new double[container.getAtomCount()][3];
        for (int i = 0; i < atomContainerGradients.length; i++) {
            atomContainerGradients[i] = new double[]{stepSize, stepSize, stepSize};
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
                    double distance;
                    distance = Math.min(Math.abs(stepSize * oneAxisGradient), .3) * Math.signum(oneAxisGradient);
                    movePoint(point, distance, i);
                    if (constraint == null || constraint.check(forcefield, container) || !initializedGradients) {
                        double difference = forcefield.calculateEnergy(container) - lastCalculatedEnergy;
                        atomGradients[i] = -difference / (oneAxisGradient);
                        if (!initializedGradients) {
                            movePoint(point, -distance, i);
                        } else {
                            if (difference > 0) {
                                movePoint(point, -distance, i);
                            } else {
                                lastCalculatedEnergy += difference;
                            }
                        }
                    }else{
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
}
