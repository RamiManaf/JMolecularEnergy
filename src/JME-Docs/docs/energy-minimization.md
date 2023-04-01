# Energy Minimization
Previously we saw how to calculate the potential energy for a molecule using Java Molecular Energy. This showed how the molecule conformation can affect the potential energy. An important application for calculating potential energy for the molecule is minimizing it to obtain the most stable and favorable structure for the molecule. Here we will show you how to energy minimize a molecule using Java Molecular Energy.
## Create a Minimizer For The Molecule
The energy minimizer use the force field to calculate the energy and try to change the molecule structure to minimize it. We have different choices for choosing energy minimization algorithm. We will use `GradientDescent` as the most popular method.
```java
MMFF94 mmff = new MMFF94(true);
GradientDescent minimizer = new GradientDescent(mmff);
```
## Assign Force Field Parameters To The Molecule
This step is required as the minimizer will keep calculating the energy of the molecule.
```java
mmff.assignParameters(container);
```
## Minimize The Energy For The Molecule
When calling `minimize` method the minimizer will block the thread untlin the process finishes. The process is an iterative cycle with each cycle called a step. The minimizer will recieve the maximum steps at wich it will stop after them. Each step the minimizer will try to change the conformation of the molecule depending on the step size. The larger the step the faster the molecule will move to the energy minimum. However, small steps help the minimizer reach the smallest possible energy but with more steps.

In this example we called the minimizer to do minimization for 1000 step with 0.1 step size. The final parameter is the threshold. This parameter determins the energy difference that when reached between two steps the minimizer will stop. E.g.: here if the last step minimized the energy by 10kcal, the minimizer will stop.
```java
minimizer.minimize(container, 1000, .1f, 10);
```