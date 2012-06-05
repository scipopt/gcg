/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file stat.h
 * @brief  Prints information about the best decomposition
 * @author Alexander Gross
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#ifndef STAT_H_
#define STAT_H_

#include "scip/type_scip.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** prints information about the best decomposition*/
SCIP_RETCODE writeDecompositionData(
   SCIP* scip           /**< SCIP data structure */
   );

/** prints information about the creation of the Vars*/
SCIP_RETCODE writeVarCreationDetails(
   SCIP* scip           /**< SCIP data structure */
);

int checkNodes(
   int* nodes,
   int node
   );


#ifdef __cplusplus
}
#endif
#endif /* STAT_H_ */
