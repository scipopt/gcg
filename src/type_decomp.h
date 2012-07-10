/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_decomp.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for decomposition information in GCG projects
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_DEC_DECOMP_H__
#define __SCIP_TYPE_DEC_DECOMP_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DecDecomp DEC_DECOMP; /**< decomposition structure */

/** type of the decomposition */
enum Dectype
{
   DEC_DECTYPE_UNKNOWN,    /**< unknown structure (used for initialization) */
   DEC_DECTYPE_ARROWHEAD,  /**< arrowhead structure (linking variables and constraints) */
   DEC_DECTYPE_STAIRCASE,  /**< staircase structure (linking variables between consecutive blocks) */
   DEC_DECTYPE_DIAGONAL,   /**< block diagonal structure (no linking variables and constraints) */
   DEC_DECTYPE_BORDERED    /**< bordered block diagonal structure (linking constraints only) */
};

typedef enum Dectype DEC_DECTYPE; /**< decomposition type */

#ifdef __cplusplus
}
#endif

#endif
