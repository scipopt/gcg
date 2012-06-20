/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_masterbranch.h
 * @brief  constraint handler for storing the branching decisions at each node of the tree
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

#ifndef CONSMASTERBRANCH_H
#define CONSMASTERBRANCH_H

#include "scip/scip.h"
#include "type_branchgcg.h"

/** returns the masterbranch constraint of the current node */
extern
SCIP_CONS* GCGconsMasterbranchGetActiveCons(
   SCIP* scip                 /**< SCIP data structure */
   );

/** creates the handler for masterbranch constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrMasterbranch(
   SCIP* scip                 /**< SCIP data structure */
   );

/** creates and captures a masterbranch constraint */
extern
SCIP_RETCODE GCGcreateConsMasterbranch(
   SCIP*       scip,          /**< SCIP data structure */
   SCIP_CONS** cons,          /**< pointer to hold the created constraint */
   SCIP_NODE*  node,          /**< node at which the constraint should be created */
   SCIP_CONS*  parentcons     /**< parent constraint */
   );


/** returns the stack and the number of elements on it */
extern
void GCGconsMasterbranchGetStack(
   SCIP*        scip,            /**< SCIP data structure */
   SCIP_CONS*** stack,           /**< return value: pointer to the stack */
   int*         nstackelements   /**< return value: pointer to int, for number of elements on the stack */
   );

/** returns the number of elements on the stack */
extern
int GCGconsMasterbranchGetNStackelements(
   SCIP* scip                 /**< SCIP data structure */
   );

/** returns the branching data for a given masterbranch constraint */
extern
GCG_BRANCHDATA* GCGconsMasterbranchGetBranchdata(
   SCIP_CONS* cons            /**< constraint pointer */
   );

/** returns the node in the B&B tree at which the given masterbranch constraint is sticking */
extern
SCIP_NODE* GCGconsMasterbranchGetNode(
   SCIP_CONS* cons            /**< constraint pointer */
   );

/** returns the masterbranch constraint of the B&B father of the node at which the
    given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetParentcons(
   SCIP_CONS* cons            /**< constraint pointer */
   );

/** returns the masterbranch constraint of the first child of the node at which the
    given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetChild1cons(
   SCIP_CONS* cons            /**< constraint pointer */
   );

/** returns the masterbranch constraint of the second child of the node at which the
    given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetChild2cons(
   SCIP_CONS* cons            /**< constraint pointer */
   );

/** returns the origbranch constraint of the node in the original program corresponding to the node
    at which the given masterbranch constraint is sticking */
extern
SCIP_CONS* GCGconsMasterbranchGetOrigcons(
   SCIP_CONS* cons            /**< constraint pointer */
   );

/** sets the origbranch constraint of the node in the master program corresponding to the node
    at which the given masterbranchbranch constraint is sticking */
extern
void GCGconsMasterbranchSetOrigcons(
   SCIP_CONS* cons,           /**< constraint pointer */
   SCIP_CONS* origcons        /**< original origbranch constraint */
   );

/** checks the consistency of the masterbranch constraints in the problem */
extern
void GCGconsMasterbranchCheckConsistency(
   SCIP* scip                 /**< SCIP data structure */
   );

/** adds initial constraint to root node */
extern
SCIP_RETCODE SCIPconsMasterbranchAddRootCons(
   SCIP*                 scip                /**< SCIP data structure */
   );


#endif
