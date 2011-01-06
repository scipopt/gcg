/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**@file    interface.h
 * @brief   scip interface for gcg (header file)
 * @author  Martin Bergner
 */
#ifndef INTERFACE_H_
#define INTERFACE_H_

#include "scip/scip.h"
extern SCIP_RETCODE GCGsetBlocksForProblem(
		SCIP *scip,
		int* blocksPerVar,
		int nblocks,
		int* masterConstraints,
		int nMasterConstraints
		);

#endif /* INTERFACE_H_ */
