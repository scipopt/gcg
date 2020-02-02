# Objective Function Coefficient Classifier (Values) {#clsvar_objvalues}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

This classifier adds variables \f$x_i\f$ to classes according to the **value of their respective
coefficient** \f$c_i\f$ in the objective function.

### Classification
#### Dantzig-Wolfe and Benders Decomposition

In both modes, the handling is as follows:

 * For each \f$x_i\f$, the respective \f$c_i\f$ is determined.
  * If the representative class for \f$c_i\f$ exists, \f$x_i\f$ is added to this class.
  * If it doesn't exist, the class for \f$c_i\f$ is created and \f$x_i\f$ is added to this class.

### Example
