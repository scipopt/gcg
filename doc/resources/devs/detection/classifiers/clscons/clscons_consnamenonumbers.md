# Name Classifier (no numbers) {#clscons_consnamenonumbers}

This classifier adds constraints to classes according to their **names** through removing the numbers and grouping **exact matches**.

### Classification
#### Dantzig-Wolfe and Benders Decomposition
In both modes, the handling is as follows:
 * For each constraint, remove all digits.
 * For each constraint \f$i\f$, the constraint name is checked for whether there is already a class for it.
  * If it exists, add constraint \f$i\f$ to it.
  * If it doesn't exist, create one and add the constraint \f$i\f$ to it.

### Example
