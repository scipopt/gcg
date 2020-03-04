# Coefficient Classifier (N-nonzeros) {#clscons_nnonzeros}

This classifier adds constraints to classes according to the **number of different variable coefficients**.

### Classification
#### Dantzig-Wolfe and Benders Decomposition
In both modes, the handling is as follows:
 * For each constraint \f$i\f$, the respective number of different variable coefficients \f$\Delta_i\f$ is determined.
  * If the representative class for \f$\Delta_i\f$ exists, constraint \f$i\f$ is added to this class.
  * If it doesn't exist, the class for \f$\Delta_i\f$ is created and constraint \f$x_i\f$ is added to this class.

### Example
