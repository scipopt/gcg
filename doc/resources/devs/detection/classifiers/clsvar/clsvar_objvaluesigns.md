# Objvaluesigns Classifier {#clsvar_objvaluesigns}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

This classifier adds variables \f$x_i\f$ to classes according to the **sign of their respective
coefficient** \f$c_i\f$ in the objective function.

### Classification
#### Dantzig-Wolfe and Benders Decomposition

In both modes, the handling is as follows:

 * Three classes are created:
  * variables \f$x_i\f$ with \f$c_i = 0\f$
  * variables \f$x_i\f$ with \f$c_i > 0\f$
  * variables \f$x_i\f$ with \f$c_i < 0\f$
 * For each \f$x_i\f$, the respective \f$c_i\f$ is determined and added to the corresponding class (see above).

### Background


### Example
