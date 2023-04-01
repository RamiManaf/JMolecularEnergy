# Calculate Potential Energy without 3D Coordinates
As you know, force fields depends on the 3D structure of the molecule for calculating potential energy. However, in some cases we don't have initial coordinates as we are importing from chemical formats that doesn't support encoding 3D coordinates such as SMILES or we are have converted implicit hydrogens to explicit ones. In that case you have to follow some steps to get an initial 3D coordinates for your molecule. We will demonstrate loading molecule from SMILES as the worst case.
## Load the molecule
Load the molecule from SMILES using CDK:
```java
IAtomContainer container = new SmilesParser(SilentChemObjectBuilder.getInstance()).parseSmiles("CCCC");
```
## Assign Atom Types
If your atoms are not propertly configured in CDK, then JME could give wrong results. This is true especially for assigning formal charges for nitrogens in your molecule.
```java
AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
```
## Add Explicit Hydrogen
After assigning atom types you can add implicit hydrogens based on atom types assigned. Then make sure to convert the implicit hydrogens to explicit ones.
```java
CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance()).addImplicitHydrogens(container);
AtomContainerManipulator.convertImplicitToExplicitHydrogens(container);
```
## Assign 3D Coordinates For The Molecule
You can assign initial 3D Coordinates manually or use the CDK `ModelBuilder3D` to assign preferred atom positions based on your molecule. 
```java
ModelBuilder3D.getInstance(TemplateHandler3D.getInstance(), "mmff94", SilentChemObjectBuilder.getInstance()).generate3DCoordinates(container, false);
```
Finally you can calculate potential energy as normal.