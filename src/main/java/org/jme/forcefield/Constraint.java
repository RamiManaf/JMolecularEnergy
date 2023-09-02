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

import javax.vecmath.Point3d;
import org.openscience.cdk.interfaces.IAtom;
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

    /**
     * creates a distance constraint between the two atoms with minimum and
     * maximum distance allowed.
     *
     * @param iAtom
     * @param jAtom
     * @param minimumDistance minimum distance in angstrom
     * @param maximumDistance maximum distance in angstrom
     * @return
     */
    public static Constraint constrainDistance(IAtom iAtom, IAtom jAtom, double minimumDistance, double maximumDistance) {
        return new Constraint() {
            private Point3d iPointCache = new Point3d(), jPointCache = new Point3d();
            private boolean checkCache;

            @Override
            public boolean check(ForceField forcefield, IAtomContainer container) {
                if (iAtom.getPoint3d().equals(iPointCache) && jAtom.getPoint3d().equals(jPointCache)) {
                    return checkCache;
                }
                iPointCache.set(iAtom.getPoint3d());
                jPointCache.set(jAtom.getPoint3d());
                double distance = iAtom.getPoint3d().distance(jAtom.getPoint3d());
                checkCache = distance >= minimumDistance && distance <= maximumDistance;
                return checkCache;
            }
        };
    }

    /**
     * creates an angle constraint for the angle between the three bonded atoms i-j-k with minimum and
     * maximum angle allowed.
     * 
     * @param iAtom
     * @param jAtom
     * @param kAtom
     * @param minimumAngle minimum angle in radians
     * @param maximumAngle maximum angle in radians
     * @return 
     */
    public static Constraint constrainAngle(IAtom iAtom, IAtom jAtom, IAtom kAtom, double minimumAngle, double maximumAngle) {
        if (iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom) == null) {
            throw new RuntimeException("the atoms must be bonded in the order i-j-k");
        }
        return new Constraint() {
            private Point3d iPointCache = new Point3d(), jPointCache = new Point3d(), kPointCache = new Point3d();
            private boolean checkCache;

            @Override
            public boolean check(ForceField forcefield, IAtomContainer container) {
                Point3d iPoint = iAtom.getPoint3d();
                Point3d jPoint = jAtom.getPoint3d();
                Point3d kPoint = kAtom.getPoint3d();
                if (iPoint.equals(iPointCache) && jPoint.equals(jPointCache) && kPoint.equals(kPointCache)) {
                    return checkCache;
                }
                iPointCache.set(iPoint);
                jPointCache.set(jPoint);
                kPointCache.set(kPoint);
                double angle = GeometryUtils.calculateAngle(iPoint, jPoint, kPoint);
                checkCache = angle >= minimumAngle && angle <= maximumAngle;
                return checkCache;
            }
        };
    }

    /**
     * creates a torsion angle constraint for the four bonded atoms i-j-k-l with minimum and
     * maximum angle allowed.
     * 
     * @param iAtom
     * @param jAtom
     * @param kAtom
     * @param lAtom
     * @param minimumAngle minimum angle in radians
     * @param maximumAngle maximum angle in radians
     * @return 
     */
    public static Constraint constrainTorsionAngle(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom, double minimumAngle, double maximumAngle) {
        if (iAtom.getBond(jAtom) == null || jAtom.getBond(kAtom) == null && kAtom.getBond(lAtom) == null) {
            throw new IllegalArgumentException("the atoms must be bonded in the order i-j-k-l");
        }
        return new Constraint() {
            private Point3d iPointCache = new Point3d(), jPointCache = new Point3d(), kPointCache = new Point3d(), lPointCache = new Point3d();
            private boolean checkCache;

            @Override
            public boolean check(ForceField forcefield, IAtomContainer container) {
                Point3d iPoint = iAtom.getPoint3d();
                Point3d jPoint = jAtom.getPoint3d();
                Point3d kPoint = kAtom.getPoint3d();
                Point3d lPoint = lAtom.getPoint3d();
                if (iPoint.equals(iPointCache) && jPoint.equals(jPointCache) && kPoint.equals(kPointCache)) {
                    return checkCache;
                }
                iPointCache.set(iPoint);
                jPointCache.set(jPoint);
                kPointCache.set(kPoint);
                lPointCache.set(lPoint);
                double angle = GeometryUtils.calculateTorsionAngle(iPoint, jPoint, kPoint, lPoint);
                checkCache = angle >= minimumAngle && angle <= maximumAngle;
                return checkCache;
            }
        };
    }

    /**
     * creates a torsion angle constraint for the four bonded atoms i-j-k/l with minimum and
     * maximum angle allowed. The angle is between the i-j-k plane and the l atom.
     * 
     * @param iAtom
     * @param jAtom
     * @param kAtom
     * @param lAtom
     * @param minimumAngle minimum angle in radians
     * @param maximumAngle maximum angle in radians
     * @return 
     */
    public static Constraint constrainOutOfPlaneAngle(IAtom iAtom, IAtom jAtom, IAtom kAtom, IAtom lAtom, double minimumAngle, double maximumAngle) {
        if (jAtom.getBond(iAtom) == null || jAtom.getBond(kAtom) == null && jAtom.getBond(lAtom) == null) {
            throw new RuntimeException("the atoms must be bonded in the order i-j-k/l where l is bonded to j");
        }
        return new Constraint() {
            private Point3d iPointCache = new Point3d(), jPointCache = new Point3d(), kPointCache = new Point3d(), lPointCache = new Point3d();
            private boolean checkCache;

            @Override
            public boolean check(ForceField forcefield, IAtomContainer container) {
                Point3d iPoint = iAtom.getPoint3d();
                Point3d jPoint = jAtom.getPoint3d();
                Point3d kPoint = kAtom.getPoint3d();
                Point3d lPoint = lAtom.getPoint3d();
                if (iPoint.equals(iPointCache) && jPoint.equals(jPointCache) && kPoint.equals(kPointCache) && lPoint.equals(lPointCache)) {
                    return checkCache;
                }
                iPointCache.set(iPoint);
                jPointCache.set(jPoint);
                kPointCache.set(kPoint);
                lPointCache.set(lPoint);
                double angle = GeometryUtils.calculateOutOfPlaneAngle(iPoint, jPoint, kPoint, lPoint);
                checkCache = angle >= minimumAngle && angle <= maximumAngle;
                return checkCache;
            }
        };
    }

}
