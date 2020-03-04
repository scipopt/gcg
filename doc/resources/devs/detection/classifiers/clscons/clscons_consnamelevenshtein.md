# Name Classifier (Levenshtein) {#clscons_consnamelevenshtein}

This classifier adds constraints to classes according to the **differences in their names** with a certain threshold corresponding to the **Levenshtein distance**.

### Classification
#### Dantzig-Wolfe and Benders Decomposition
In both modes, the handling is as follows:
 * For each constraint \f$i\f$, the constraint name is read and the distance to each other constraint name is calculated (complexity \f$\mathcal{O}(n^2)\f$)
  * A breadth first search is conducted until every constraint is assigned to a class
   * Iff the Levenshtein distance bewtween two constraints is higher than `connectivity` (by default `1`), they aren't grouped into the same class.
 * If the number of classes exceeds `nmaxconss` (by default `5000`), the classifier will abort.

### Example
