# Type Classifier (SCIP Variable Type) {#clsvar_scipvartypes}

This classifier adds variables to classes according to their **domain given by SCIP**.
There are four different variable types supported:

 * Binary (`bin`)
 * Integer (`int`)
 * Implicit Integer (`impl`)
 * Continuous (`cont`)

### Classification
#### Dantzig-Wolfe Decomposition

The handling for Dantzig-Wolfe decomposition mode is as follows:

 * A new class is generated for each of the aforementioned types.
 * By default, the class with name `newVartype` is generated

#### Benders' Decomposition
The handling for Benders' decomposition is as follows:

 * If `detection/benders/onlycontsubpr` is set:
  * Integer variables will be added to the linking variables
  * Continuous variables will be added to the block variables

 * If `detection/benders/onlybinmaster` is set:
  * Binary variables will be added to the linking variables
  * Integer Variables, Implicit Integer Variables and Continuous Variables will be added to the block variables

### Example
Different domains of variables usually hint to a different purpose of those.
Thus, we find an example for a problem that one usually models using different variable types.<br>
We look at a "text-book" mixed-integer program for the Lot-Sizing problem.

**Given:** single product, \f$T\f$ periods, holding cost \f$l\f$ per unit, setup cost \f$r\f$ per lot, demand \f$b_t\f$ in period \f$t\f$ <br>
**Problem:** how much to produce in each period in order to satisfy demands at minimum setup and holding cost  <br>
**Variables:** \f$x_t\f$ how much to store in time \f$t\f$, \f$y_t\f$ how much to produce in time \f$t\f$, \f$z_t\f$ produce in \f$t\f$

\f{align}{
  & \text{min}
  & & \sum_{t=0}^T(l y_t+r z_t) \\
  & \text{s.t.} & & y_{t-1}+x_t = b_t+y_t && t=0,\cdots,T-1 \\
  & & & x_t\leq z_t\cdot M && t=0,\cdots,T \\
  & & & \mathbf{x_t, y_t \in \mathbb{Z}_{+}} && t=0,\cdots,T \\
  & & & \mathbf{z_t\in\{0,1\}} && t=0,\cdots,T
  \f}

We see that there are two different types of variables, which is detected by GCG:

    Varclassifier "vartypes" yields a classification with 2 different variable classes
