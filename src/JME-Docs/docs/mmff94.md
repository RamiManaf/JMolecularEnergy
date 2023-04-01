# MMFF94 Force Field
MMFF94  is an empirical force field that is optimized for calculating potential energy of organic molecules. The force fields have two subdivisions MMFF94 & MMFF94s (stands for static). MMFF94 is the first one that is intended for simulation while MMFF94s is optimized for energy minimization. You can choose which one through the constructor boolean parameter of the force field. Use `true` for MMFF94s and `false` for MMFF94.
```java
MMFF94 mmff94 = new MMFF94(false);
```
## MMFF94 Energy Terms
The force field use several terms to calculate potential energy of the molecule.
![Force Field](img/forcefield.png)
The forces are divided to bonded and non-bonded interactions:
* Bonded Interactions:
  * Bond Stretching: for calculating energy based on the distance between two bonded atoms.
  * Angle Bending: the angle between three bonded atoms.
  * Stretch Bend Interaction: A force that increase the bonds when the angle between three bonded atoms decrease. This is used to correct the structure of the molecule.
  * Torsion: This term forces angle formed between four bonded atoms to be at certain angle. This angle is also called dihydral angle.
  * Out Of Plane: This term is used to describe the effect of the angle formed when three atoms at least connect to a central atom. Two end atoms and the central atom can form a plane which the third end atom is out of it by certain angle. This term is used usually to correct the Sp2 hybridized central atom planarity.
* Non-bonded Interactions:
  * Van Der Waals: Calculate the van der waals force between two non-bonded atoms.
  * Electrostatic Interaction: This term calculates electrostatic interaction energy between atoms. Note that not only atoms with formal charges are included as MMFF94 assign partial charges for unionized atoms.
 ## Implementation
 MMFF94 instance when first created loads all the force fields parameters needed. Before calculating energy to any molecule you have to assign parameters to this molecule first. This process prevent MMFF94 from searching for proper paramters for the molecule each time you want to calculate the energy.
 
  The partial charges are assigned to the molecule and you can retrive them using the following code:
 ```java
 atomContainer.getAtom(0).getCharge();
 ```
 For more informations about MMFF94 implementation you can see the original papers:

* [Merck molecular force field. I. Basis, form, scope, parameterization, and performance of MMFF94](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<490::AID-JCC1>3.0.CO;2-P)
* [Merck molecular force field. II. MMFF94 van der Waals and electrostatic parameters for intermolecular interactions](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6%3C520::AID-JCC2%3E3.0.CO;2-W)
* [Merck molecular force field. III. Molecular geometries and vibrational frequencies for MMFF94.](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<553::AID-JCC3>3.0.CO;2-T)
* [Merck molecular force field. IV. conformational energies and geometries for MMFF94](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<587::AID-JCC4>3.0.CO;2-Q)
* [Merck molecular force field. V. Extension of MMFF94 using experimental data, additional computational data, and empirical rules](https://doi.org/10.1002/(SICI)1096-987X(199604)17:5/6<616::AID-JCC5>3.0.CO;2-X)
* [MMFF VI. MMFF94s option for energy minimization studies](https://doi.org/10.1002/(SICI)1096-987X(199905)20:7<730::AID-JCC8>3.0.CO;2-T)
* [MMFF VII. Characterization of MMFF94, MMFF94s, and other widely available force fields for conformational energies and for intermolecular-interaction energies and geometries](https://doi.org/10.1002/(SICI)1096-987X(199905)20:7<720::AID-JCC7>3.0.CO;2-X)
