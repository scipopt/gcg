/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   struct_solver.h
 * @brief  datastructures for solvers
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SOLVER_H__
#define __SCIP_STRUCT_SOLVER_H__

#include "type_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/** branching rule */
struct GCG_Solver
{
   char*                 name;
   char*                 description;
   int                   priority;
   GCG_SOLVERDATA*       solverdata;     
   GCG_DECL_SOLVERSOLVE((*solversolve));
   GCG_DECL_SOLVERSOLVEHEUR((*solversolveheur));
   GCG_DECL_SOLVERFREE   ((*solverfree));
   GCG_DECL_SOLVERINIT   ((*solverinit));
   GCG_DECL_SOLVEREXIT   ((*solverexit));
   GCG_DECL_SOLVERINITSOL((*solverinitsol));
   GCG_DECL_SOLVEREXITSOL((*solverexitsol));
   SCIP_CLOCK*           optfarkasclock;
   SCIP_CLOCK*           optredcostclock;
   SCIP_CLOCK*           heurfarkasclock;
   SCIP_CLOCK*           heurredcostclock;
   int                   optfarkascalls;
   int                   optredcostcalls;
   int                   heurfarkascalls;
   int                   heurredcostcalls;
};


#ifdef __cplusplus
}
#endif

#endif
