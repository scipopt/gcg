/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_decomp.h
 * @brief  structure information for decomposition information in GCG projects
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef STRUCT_DECOMP_H_
#define STRUCT_DECOMP_H_
#include "scip/scip.h"
#include "type_decomp.h"

#ifdef __cplusplus
extern "C" {
#endif


/** decomposition structure information */
struct DecDecomp
{
   int            nblocks;       /**< number of blocks in this decomposition */
   SCIP_VAR***    subscipvars;   /**< two dimensional array of variables in each block */
   int*           nsubscipvars;  /**< array of number of variables in each block */
   SCIP_CONS***   subscipconss;  /**< two dimensional array of constraints in each block */
   int*           nsubscipconss; /**< array of number of constraints in each block */
   SCIP_CONS**    linkingconss;  /**< array of constraints linking the blocks */
   int            nlinkingconss; /**< number of linking constraints */
   SCIP_VAR**     linkingvars;   /**< array of variables linking the blocks */
   int            nlinkingvars;  /**< number of linking variables */
   SCIP_HASHMAP*  vartoblock;    /**< hashmap mapping variables to their blocks (from 1 to nblocks) */
   SCIP_HASHMAP*  constoblock;   /**< hashmap mapping constraints to their blocks (from 1 to nblocks) */
   SCIP_HASHMAP*  varindex;      /**< hashmap mapping variables to indeces for a visual ordering */
   SCIP_HASHMAP*  consindex;     /**< hashmap mapping constraints to indices for visual ordering */
   DEC_DECTYPE    type;          /**< type of the decomposition */
};

#ifdef __cplusplus
}
#endif

#endif /* STRUCT_DECOMP_H_ */
