# Original Problem Separator {#sepa_original}

### Overview
This separator finds cuts in the original problem and then transfers them to the master.

### Theoretical Background
Consider a general linear program with two sets of constraints:
\f{align}{
\min \quad& c^T x &\\
s.t. \quad &  Ax &\ge b \\
 & Dx &\ge d \\
 & x &\ge 0 \\
\f}
Assume that GCG chose the constraints \f$ Dx \ge d \f$ to derive the decomposition from. Set \f$ X = \{x \colon Dx \ge d, x \ge 0\} \f$ and let \f$P\f$ be the set of extreme points of \f$ X \f$. Then the restricted master problem has the form
 \f{align}{
\min & \sum \limits_{p \in P’} c_p  \lambda_p &\\
s.t. & \sum \limits_{p \in P’}  a_p\lambda_p &\ge b & [\pi] \\
 & \sum \limits_{p \in P’} \lambda_p&= 1& [\pi_0] \\
 & \lambda &\ge 0 \\
\f}
with \f$ P‘ \subseteq P , a_p = A p, c_p = c^T p \f$. The corresponding pricing problem look like this:
\f{align}{
\min \quad & (c^T -\pi^TA) x - \pi_0&\\
s.t. \quad & Dx &\ge d \\
 & x &\ge 0 \\
\f}
In our setting we find additional constraints \f$ Fx\ge g\f$ for the original problem which now has the form:
\f{align}{
\min \quad& c^T x &\\
s.t. \quad&  Ax &\ge b \\
 & Dx &\ge d \\
 & Fx & \ge g \\
 & x &\ge 0 \\
\f}
If we transfer them to our restricted master problem, we obtain
 \f{align}{
\min & \sum \limits_{p \in P’} c_p  \lambda_p &\\
s.t. & \sum \limits_{p \in P’}  a_p\lambda_p &\ge b & [\pi] \\
 & \sum \limits_{p \in P’} \lambda_p&= 1& [\pi_0] \\
 & \sum \limits_{p \in P’}  a_p\lambda_p &\ge b & [\pi] \\
 & \sum \limits_{p \in P’}  f_p\lambda_p &\ge g & [\alpha] \\

 & \lambda &\ge 0 \\
\f}
where \f$ f_p = F_p \f$. At last, we need to correct the pricing problem
\f{align}{
\min & (c^T -\pi^TA - \alpha^T F) x - \pi_0&\\
s.t. & Dx &\ge d \\
 & x &\ge 0 \\
\f}

### Implementation
The separator 'original' does implement exactly the above described procedure. In case we do not have an aggregated pricing problem, the separator calls the intern SCIP separators for the original problem and then transfers the found cuts to the master problem.
This separator can be tuned like any normal SCIP separator, see the [SCIP documentation](https://www.scipopt.org/doc-7.0.1/html/SEPA.php). **Please disable this separator if you implemented an own solver for the pricing problems**, because it can unintentionally destroy the special structure you want to use.
