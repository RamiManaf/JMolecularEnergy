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
import javax.vecmath.Vector3d;

/**
 * A utils class for geometrical calculations on the atoms
 * 
 * @author Rami Manaf Abdullah
 */
public class GeometryUtils {

    /**
     * Calculate the angle between the i-j-k connected points in radians
     * @param i
     * @param j
     * @param k
     * @return 
     */
    public static double calculateAngle(Point3d i, Point3d j, Point3d k) {
        double ux = i.x - j.x;
        double uy = i.y - j.y;
        double uz = i.z - j.z;
        double vx = k.x - j.x;
        double vy = k.y - j.y;
        double vz = k.z - j.z;
        return Math.acos(clampToOne((ux * vx + uy * vy + uz * vz) / (i.distance(j) * k.distance(j))));
    }
    
    /**
     * Calculate the torsion angle between the i-j-k-l connected points in radians
     * @param i
     * @param j
     * @param k
     * @param l
     * @return 
     */
    public static double calculateTorsionAngle(Point3d i, Point3d j, Point3d k, Point3d l) {
        Vector3d vectorIJ = new Vector3d(i.x - j.x, i.y - j.y, i.z - j.z);
        Vector3d vectorKJ = new Vector3d(k.x - j.x, k.y - j.y, k.z - j.z);
        Vector3d vectorJK = new Vector3d(vectorKJ);
        vectorJK.negate();
        Vector3d vectorLK = new Vector3d(l.x - k.x, l.y - k.y, l.z - k.z);
        Vector3d n1 = new Vector3d(), n2 = new Vector3d();
        n1.cross(vectorIJ, vectorKJ);
        n2.cross(vectorJK, vectorLK);
        return Math.acos(clampToOne(n1.dot(n2) / (n1.length() * n2.length())));
    }

    private static double clampToOne(double n) {
        return n > 1 ? 1 : (n < -1 ? -1 : n);
    }

    /**
     * Calculate the angle between the i-j-k/l connected points in radians where l is the out of plane vector
     * @param i
     * @param j
     * @param k
     * @param l
     * @return 
     */
    public static double calculateOutOfPlaneAngle(Point3d i, Point3d j, Point3d k, Point3d l) {
        Vector3d ijVector = new Vector3d(), kjVector = new Vector3d(), ljVector = new Vector3d(), perpendicular = new Vector3d();
        ijVector.sub(i, j);
        ijVector.normalize();
        kjVector.sub(k, j);
        kjVector.normalize();
        perpendicular.cross(ijVector, kjVector);
        ljVector.sub(l, j);
        ljVector.normalize();
        return Math.asin(perpendicular.dot(ljVector) / Math.sin(GeometryUtils.calculateAngle(i, j, k)));
    }
}
