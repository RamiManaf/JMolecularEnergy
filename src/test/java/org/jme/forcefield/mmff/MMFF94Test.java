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
package org.jme.forcefield.mmff;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.config.atomtypes.OWLAtomTypeMappingReader;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

/**
 * based on validation suites:
 * https://ccl.net/cca/data/MMFF94/index.shtml
 * https://server.ccl.net/cca/data/MMFF94s/index.shtml
 * 
 *
 * @author Rami Manaf Abdullah
 */
public class MMFF94Test {

    private static MMFF94 mmff94 = new MMFF94(false);
    private static MMFF94 mmff94s = new MMFF94(true);
    private static Map<String, String> sybylToCdk = new HashMap<>();

    @BeforeClass
    public static void init() {
        System.setProperty("cdk.logging.level", "error");
        Map<String, String> map = new OWLAtomTypeMappingReader(new InputStreamReader(MMFF94Test.class.getResourceAsStream("/org/openscience/cdk/dict/data/cdk-sybyl-mappings.owl"))).readAtomTypeMappings();
        for (Map.Entry<String, String> entry : map.entrySet()) {
            sybylToCdk.put(entry.getValue(), entry.getKey());
        }
    }

    @Test
    public void testMMFF94() throws IOException, CDKException {
        ArrayList<Float> energyValidation = parseEnergyValidationFile("MMFF94.energies");
        List<IAtomContainer> containers = ChemFileManipulator.getAllAtomContainers(new Mol2Reader(getClass().getResourceAsStream("MMFF94_dative.mol2")).read(new ChemFile()));
        fixMMFF94Molecules(containers);
        double error = 0;
        int n = 0;
        int total = 0;
        for (int i = 0; i < containers.size(); i++) {
            IAtomContainer container = containers.get(i);
            double energy = calculateMMFF94Energy(container, mmff94);
            total++;
            if (Math.abs(energyValidation.get(i) - energy) > .01) {
                System.out.println(i + "-" + container.getTitle() + ":\tExpected:\t" + energyValidation.get(i) + "\tFound:\t" + energy);
                error += Math.abs(energyValidation.get(i) - energy);
                n++;
            }
        }
        System.out.println("Error:\t" + (error) + "\tN:" + n);
        System.out.println("Percentage:\t" + (n / (float) total) * 100 + "\t" + n + "/" + total);
        Assert.assertEquals(0, n);
    }

    private void fixMMFF94Molecules(List<IAtomContainer> containers) {
        for (int i = 0; i < containers.size(); i++) {
            IAtomContainer container = containers.get(i);
            if (i == 726) {
                //fix bond order
                container.getBond(6).setOrder(IBond.Order.DOUBLE);
            }
            if (i == 255 || i == 324 || i == 326) {
                //fix amine double bond in sulfonamide
                for (IAtom atom : container.atoms()) {
                    if (atom.getAtomicNumber() == 16) {
                        for (IAtom atom2 : container.getConnectedAtomsList(atom)) {
                            if (atom2.getAtomicNumber() == 7 && atom2.getBondCount() == 2) {
                                List<IAtom> connectedAtoms = container.getConnectedAtomsList(atom2);
                                IAtom other = connectedAtoms.get(0) != atom ? connectedAtoms.get(0) : connectedAtoms.get(1);
                                if (atom2.getBond(atom).getOrder().equals(IBond.Order.SINGLE)
                                        && atom2.getBond(other).getOrder().equals(IBond.Order.SINGLE)) {
                                    container.getBond(atom, atom2).setOrder(IBond.Order.DOUBLE);
                                }
                            }
                        }
                    }
                }
            } else if (i == 445) {
                //assign formal charges
                container.getAtom(0).setFormalCharge(-1);
            } else if (i == 737) {
                container.getAtom(2).setFormalCharge(-1);
            } else if (i == 741) {
                container.getAtom(3).setFormalCharge(1);
            } else if (i == 742) {
                container.getAtom(9).setFormalCharge(2);
            } else if (i == 743) {
                container.getAtom(9).setFormalCharge(2);
            } else if (i == 744) {
                container.getAtom(9).setFormalCharge(3);
            } else if (i == 752) {
                container.getAtom(9).setFormalCharge(2);
            }
        }
    }

    @Test
    public void testMMFF94s() throws IOException, CDKException {
        ArrayList<Float> energyValidation = parseEnergyValidationFile("MMFF94s.energies");
        List<IAtomContainer> containers = ChemFileManipulator.getAllAtomContainers(new Mol2Reader(getClass().getResourceAsStream("MMFF94s_dative.mol2")).read(new ChemFile()));
        double error = 0;
        int n = 0;
        int total = 0;
        for (int i = 0; i < energyValidation.size(); i++) {
            IAtomContainer container = containers.get(i);
            double energy = calculateMMFF94Energy(container, mmff94s);
            total++;
            if (Math.abs(energyValidation.get(i) - energy) > .01) {
                System.out.println(i + "-" + container.getTitle() + ":\tExpected:\t" + energyValidation.get(i) + "\tFound:\t" + energy);
                error += Math.abs(energyValidation.get(i) - energy);
                n++;
            }
        }
        System.out.println("Error:\t" + (error) + "\tN:" + n);
        System.out.println("Percentage:\t" + (n / (float) total) * 100 + "\t" + n + "/" + total);
        Assert.assertEquals(0, n);
    }

    private double calculateMMFF94Energy(IAtomContainer container, MMFF94 mmff) throws CDKException {
        container.atoms().forEach((IAtom t) -> {
            t.setImplicitHydrogenCount(0);
            for (IAtom atom : container.atoms()) {
                //assign proper formal charge
                if (atom.getAtomicNumber() == 7) {
                    int sum = 0;
                    for (IBond bond : atom.bonds()) {
                        sum += bond.getOrder().numeric();
                    }
                    if (sum == 4) {
                        atom.setFormalCharge(1);
                    }
                } else if (atom.getAtomicNumber() == 16) {
                    //fix bond order
                    for (IAtom atom2 : container.getConnectedAtomsList(atom)) {
                        if ((((atom2.getAtomicNumber() == 8 || atom2.getAtomicNumber() == 16) && atom2.getBondCount() == 1))
                                && atom2.getBond(atom).getOrder().equals(IBond.Order.SINGLE)) {
                            container.getBond(atom, atom2).setOrder(IBond.Order.DOUBLE);
                        }
                    }
                }
            }
        });
        mmff.assignParameters(container);
        return mmff.calculateEnergy(container);
    }

    private ArrayList<Float> parseEnergyValidationFile(String file) throws IOException {
        ArrayList<Float> energyValidationFile = new ArrayList<>();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(getClass().getResourceAsStream(file), StandardCharsets.UTF_8))) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (line.isEmpty() || line.charAt(0) == '*') {
                    continue;
                }
                energyValidationFile.add(Float.parseFloat(line.split("\t")[1]));
            }
        }
        return energyValidationFile;
    }

}
