# How to create generic master cuts {#generic-mastercuts}
> Sometimes, we developers would like to add constraints to the reformulation that have no
> (known) counterpart in the original problem. To facilitate this, \GCG provides an interface to add
> such constraints. This page defines the structure of such constraints, which we call _generic
> master cuts_, and provides an example of how to add them to the reformulation.

## Theoretical background
Generic master cuts are constraints that are added to the master problem of the Dantzig-Wolfe
reformulation. Ordinarily, we require such constraints when no counterpart in the original
problem exists. In general, a generic master cut consists of some constraint to be added to the
master:

\f{align}{
\sum_{p \in P} f(p) \lambda_p + \sum_{r \in R} f(r) \lambda_p \leq h \quad \left[ \gamma \right]
\f}

where \(f\) is a function used to determine the coefficient of each column in the constraint.

To ensure this constraint is valid and correct throughout the column generation process, the
relevant pricing problems are now responsible for generating the correct column coefficient of
the master constraint for each new column. In doing so, we can also consider the dual value of the
constraint in the pricing problem, ensuring only columns with negative reduced cost are generated.
For this we create a new variable \(y\), the _coefficient variable_ in the pricing problem, which
we will force to take the coefficient value of the new column in the master constraint. Expressing
this constraint \(y = f(x)\) might require additional constraints and auxiliary variables.

\f{align}{
& \text{min}
& & \left( c^\top - \pi^\top A \right) x - \gamma y - \pi_0 & \\
& \text{s.t.} & & D x &\geq d \\
& & & y &= f(x) \\
& & & x &\in X \\
& & & y &\in Y
\f}

## Interface
In \GCG, a generic mastercut (`GCG_MASTERCUTDATA`) is a wrapper around either a `SCIP_CONS` or a
`SCIP_ROW` to be added to the master problem, in addition to one set of pricing modifications
(`GCG_PRICINGMODIFICATION`) for each relevant pricing problem. A pricing problem consists of a
coefficient variable (`SCIP_VAR`), as well as any additional constraints (`SCIP_CONS`) and auxiliary
variables (`SCIP_VAR`) required to express the constraint \(y = f(x)\). It is required that these
variables have the type `GCG_VARTYPE_INFERREDPRICING`. We advise using the following interfaces:
 - create inferred pricing variables with `GCGcreateInferredPricingVar()`
 - create modifications with `GCGpricingmodificationCreate()`
 - create a generic mastercut around a `SCIP_CONS` with `GCGmastercutCreateFromCons()`
 - create a generic mastercut around a `SCIP_ROW` with `GCGmastercutCreateFromRow()`


## Usage
```C
int blocknr = ...; // which pricing problem to add the pricing modifications to
SCIP_VAR* x1 = ...; // some existing variable in the pricing problem
SCIP_VAR* x2 = ...; // some existing variable in the pricing problem

/* create the master constraint (alternatively SCIP_ROW) */
SCIP_CONS* mastercons = NULL;
// create the constraint here

/* create the coefficient variable of type GCG_VARTYPE_INFERREDPRICING */
SCIP_VAR* coefvar = NULL;
char coefvarName[SCIP_MAXSTRLEN];
(void) SCIPsnprintf(coefvarName, SCIP_MAXSTRLEN, "y");
// here, we create a binary variable
SCIP_CALL( GCGcreateInferredPricingVar(pricingscip, coefvar, coefvarName, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, blocknr) );

/* create additional constraints (and auxiliary variables) */
SCIP_CONS** additionalcons = NULL;
int nadditionalcons = 2;
SCIP_CALL( SCIPallocBlockMemoryArray(pricingscip, &additionalcons, nadditionalcons) );
char consname[SCIP_MAXSTRLEN];
// y + x1 <= 1
(void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cons1");
SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[0], consname, coefvar, x1, 1, -SCIPinfinity(pricingscip), 1,
            TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
// y + x2 <= 1
(void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cons2");
SCIP_CALL( SCIPcreateConsVarbound(pricingscip, &additionalcons[1], consname, coefvar, x2, 1, -SCIPinfinity(pricingscip), 1,
            TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

/* create the pricing modification for blocknr */
GCG_PRICINGMODIFICATION** pricingmods = NULL;
int npricingmods = 1;
SCIP_CALL( SCIPallocBlockMemoryArray(masterscip, &pricingmods, npricingmods) );
SCIP_CALL( GCGpricingmodificationCreate(
   masterscip,
   &pricingmods[0],
   blocknr,
   coefvar,
   NULL, // we have no auxiliary variables
   0,
   additionalcons,
   nadditionalcons
) );

/* create the master cut */
GCG_MASTERCUTDATA* mastercutdata = NULL;
SCIP_CALL( GCGmastercutCreateFromCons(
   masterscip,
   mastercutdata,
   mastercons,
   pricingmods,
   npricingmods,
   branchdata,
   mastercutGetCoeffCompBnd
) );
```

\GCG needs to be made aware of the new master cut. This can be easily achieved by implementing callbacks of the other interfaces, e.g. for branching, implement:
```C
static
GCG_DECL_BRANCHGETMASTERCUT(branchGetMastercutCompBnd)
{
   // grab the mastercut, e.g. from the branchdata
   *mastercutdata = branchdata->mastercutdata;
   return SCIP_OKAY;
}
```

Don't forget to free the memory at some point:
```C
SCIP_CALL( GCGmastercutFree(masterscip, &mastercutdata) );
```