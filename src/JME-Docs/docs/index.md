# Quick Start Tutorial
In this tutorial, we'll show you how to use Java Molecular Energy to calculate the potential energy of a molecule using the MMFF94 force field.
## Install Java Molecular Energy
Before we can use Java Molecular Energy, we need to install it. Just add it as a dependency in you maven or gradle build system.
```xml
<dependency>
    <groupId>io.github.ramimanaf</groupId>
    <artifactId>jme</artifactId>
    <version>1.0.0</version>
</dependency>
```
## Load the Molecule
Next, we'll load a molecule into our program. You can use CDK readers to import your molecule as `AtomContainer`:
```java
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

ChemFileManipulator.getAllAtomContainers(new Mol2Reader(getClass().getResourceAsStream("molecule.mol2")).read(new ChemFile()));
IAtomContainer container = containers.get(0);
```
Make sure that your molecule doesn't have implicit hydrogens and all it's atoms have 3D coordinates. If you still have them you can see the following page.
## Calculate The Potential Energy in MMFF94
Create a MMFF94 object and assign MMFF94 parameters to it. After that you can calculate the potential energy for your molecule:
```java
MMFF94 mmff94 = new MMFF94(false);
mmff.assignParameters(container);
System.out.println("Energy: "+mmff.calculateEnergy(container)+" kcal/mol");
```
And that's it!
