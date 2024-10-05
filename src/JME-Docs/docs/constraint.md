# Constraints

Constraints are locks that can prevent the molecule from violating certain conditions.
 These locks are useful in certain situations such as in energy minimization and simulaiton.
 The locks can be based on molecule confomration such as fixing two atoms at certain
 distance from each other or allowing a range of angle that a dihydral angle is permitted to be
 in.

```java
import org.jme.constraint.*;
...
ConformationalConstraintFactory.constrainDistance(atom1, atom2, 1.4, 1.6);//locks the bond distance between 1.4 to 1.6 Å
```

`Constraints` can be combined with boolean operators such as not, and, or which are provided
 as methods in the `Constraints` interface.

```java
Constraint constraint = bondConstraint.and(angleConstraint);
```
