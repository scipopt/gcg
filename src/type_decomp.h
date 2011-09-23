/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_decomp.h 198 2011-01-06 16:58:56Z ggamrath $"

/**@file   type_decomp.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for branching rules in gcg projects
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_DECDECOMP_H__
#define __SCIP_TYPE_DECDECOMP_H__

typedef struct DecDecomp DECDECOMP;

#ifdef __cplusplus
extern "C" {
#endif


/**
 * initialize data for a detector
 */
#define DEC_DECL_INITDETECTOR(x) SCIP_RETCODE x (SCIP* scip)

/**
 * free data from a detector
 */
#define DEC_DECL_EXITDETECTOR(x) SCIP_RETCODE x (SCIP* scip)

/**
 * detects the structure of a the problem
 */
#define DEC_DECL_DETECTSTRUCTURE(x) SCIP_RETCODE x (SCIP* scip, DEC_DETECTORDATA* decdata, SCIP_RESULT* result)



#ifdef __cplusplus
}
#endif

#endif
