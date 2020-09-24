# Subset Row Separator {#sepa_subsetrow}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**
### Overview
This separator is our first approach to add strengthening cuts directly to the master problem. At first we just tried to add [Subsetrow cuts]( https://www.researchgate.net/publication/220243923_Subset-Row_Inequalities_Applied_to_the_Vehicle-Routing_Problem_with_Time_Windows). Now it received a small update and uses the heuristic of [Koster et al.]( https://www.researchgate.net/publication/220770325_Algorithms_to_Separate_0_12-Chvatal-Gomory_Cuts) to generate arbitrary zero-half cuts in the solving process. 
> Please note that the subsetrowcut separator is disabled by default! \n
> Look at the 'Algorithmic Details' and 'Parameters' sections to see why it is enabled and how to enable it.

### Theoretical background
Consider a general linear program with two sets of constraints:
\f{align}{
    \min\quad           & c^T x &\\
    \text{s.t.} \quad   &  Ax &\ge b \\
                        & Dx &\ge d \\
                        & x &\ge 0 \\
\f}
Assume that GCG chose the constraints \f$ Dx \ge d \f$ to derive the decomposition from. Set \f$ X = \{x \colon Dx \ge d, x \ge 0\} \f$ and let \f$P\f$ be the set of extreme points of \f$X\f$. Then the restricted master problem has the form 

\f{align}{ 
    \min\quad           & \sum \limits_{p \in P'} c_p \lambda_p & \\
    \text{s.t.}\quad    & \sum \limits_{p \in P'}  a_p\lambda_p &\ge b  & \qquad [\pi] \\
                        & \sum \limits_{p \in P'} \lambda_p&= 1         & \qquad [\pi_0] \\
                        & \lambda &\ge 0 \\
\f}

with \f$ P' \subseteq P , a_p := A p, c_p = c^T p \f$. The corresponding pricing problem look like this: 

\f{align}{ 
    \min\quad           & (c^T -\pi^TA) x - \pi_0\\
    \text{s.t.}\quad    & Dx &\ge d \\
                        & x  &\ge 0 \\
\f}

In our setting we add a Gomory cut (or more accurately a zero-half cut) to the master:

\f{align}{
    \min\quad           & \sum \limits_{p \in P'} c_p  \lambda_p   &\\
    \text{s.t.}\quad    & \sum \limits_{p \in P'}  a_p\lambda_p    &\ge b   & \qquad [\pi] \\
                        & \sum \limits_{p \in P'} \lambda_p        &= 1     & \qquad [\pi_0] \\
                        & \sum \limits_{p \in P'}  g_p\lambda_p    & \le h  & \qquad [\beta \le 0] \\
                        & \lambda &\ge 0 \\
\f}

Note that \f$ g_p = \lfloor A x \rfloor =: g(A)\f$. As we have a restricted master problem, we still need to price variables and while doing so pay attention to the new constrain we added. Therefore, we need to adjust our pricing problem:
\f{align}{
    \min\quad        & (c^T -\pi^TA) x - \pi_0 - \beta g(A)\\
    \text{s.t.}\quad & Dx &\ge d \\
                     & x  &\ge 0 \\
\f}
If we do it exactly this way, it results in a nonlinear cost function. But we can apply a simple trick to save it. To this end, we introduce a new integer variable \f$y\f$ to our pricing problem and reformulate it the following way:
\f{align}{
    \min\quad        & (c^T -\pi^TA) x - \pi_0 - \beta y&\\
    \text{s.t.}\quad & Dx &\ge d \\
                     & y &\le Ax +1 -\varepsilon \\
                     & x &\ge 0 \\
                     & y &\in \mathbb{Z} \\
\f}
Because \f$- \beta \ge 0\f$, y will be set as small as possible. If we choose \f$ \varepsilon > 0\f$ small enough, we will obtain \f$ y = \lfloor A x \rfloor \f$ in any optimal solution.   
Observe that all above considerations can be extended to general linear integer programs and to problems with multiple pricing problems.  
### Algorithmic Details 
From the algorithmic point of view the above theory is most of all problematic because of the \f$ \varepsilon \f$ and the resulting numerical instability. If we set \f$ \varepsilon < 10^{-5} \f$, SCIP tends to ignore it completely. Hence, we need to set \f$ \varepsilon \f$ quite large, sometimes resulting in wrong results. To prevent it from happening, the separator is disabled by default and considered an experimental feature. It can provide amazing good results for some cases, so if you are stuck at some point it is worth trying out if the separator works in your case. In order to do this please set the ```SEPA_FREQ``` variable of the separator to ```1```, see the [SCIP documentation]( https://www.scipopt.org/doc-7.0.1/html/SEPA.php) for more details. 

[Scientists, Developers]
Please keep in mind that the separator is only called in the root node, as the heuristic we use to find a zero-half cut is currently only tested for this case. If you want to change this, please comment out the relevant if clause in the ```SCIP_DECL_SEPAEXECLP```. Further, the separator always adds at most one zero-half cut in a single iteration, no matter how many the corresponding heuristic found. This can be changed in the method ```selectConstraintsForSubsetrowCut_ZEROHALF_Kosteretall```, a detailed tutorial for this will lead to far at this point. In addition to the heuristic of Koster et al. we provide to further methods: 
*  ```selectConstraintsForSubsetrowCut_RANDOM``` selects three random constrains for a Subsetrowcut.  
* ```selectConstraintsForSubsetrowCut``` lets you specify three indices for rows from which the subsetrowcut will be created.   
We do not recommend using them, but they provide an example how to implement an own method for selecting rows for the zero-half cuts. If you want to create general Gomory cuts and not just zero-half cuts, you just need to modify the ```calculateMultipliersAndRHS``` method to meat your needs. 

At last, we would like to direct your attention to the fact that the separator destroys any kind of structure the pricing problem. Because of this, **we absolutely don't recommend to use this separator if you want to use specialized solvers for the pricing problems!**

### Parameters
Here we list all the different parameters you can tune for this separator which are not already listed in the [SCIP documentation]( https://www.scipopt.org/doc-7.0.1/html/SEPA.php):
* ```SUBSETROW_EPSILON``` represents the value of the \f$ \varepsilon\f$. Its default value is ```10^{-4}```. **Do not set it smaller!** If possible, we recommend to even raise its value to ```0.1```.  
* ```MAX_NUMBER_SUBSETROWCUTS``` is the maximal number of zero-half cuts we will create. Its default value is ```1```. If your root relaxation does not take up too much time, we do not recommend increasing this number. Otherwise it can make sense to set it to ```100``` or even larger, since every cut we generate will definitely cut of some fractional solutions.
* ```STARTMAXCUTS``` is the initial size of the array containing the generated cuts (more exactly pointers to them) from our separators. Because it is quite bothersome to reallocate it, the default value is ```100``` in anticipation that the user will want to generate more that one cut. As long as ```MAX_NUMBER_SUBSETROWCUTS``` is not to large, try to keep \f$ \texttt{STARTMAXCUTS} \ge \texttt{MAX_NUMBER_SUBSETROWCUTS}\f$.

