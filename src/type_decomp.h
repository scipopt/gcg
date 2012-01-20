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
 * @brief  type definitions for branching rules in gcg projects
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_DECDECOMP_H__
#define __SCIP_TYPE_DECDECOMP_H__

typedef struct DecDecomp DECDECOMP;
enum Dectype
{
   DEC_DECTYPE_ARROWHEAD, DEC_DECTYPE_STAIRCASE, DEC_DECTYPE_DIAGONAL, DEC_DECTYPE_BORDERED, DEC_DECTYPE_UNKNOWN
};

typedef enum Dectype DEC_DECTYPE;

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif

#endif
