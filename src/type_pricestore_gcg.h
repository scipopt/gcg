/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_pricestore.h
 * @brief  type definitions for storing priced cols
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_TYPE_PRICESTORE_H__
#define __GCG_TYPE_PRICESTORE_H__

#ifdef __cplusplus
extern "C" {
#endif

/** possible settings for specifying the solution for which cuts are selected */
enum GCG_Efficiacychoice
{
   GCG_EFFICIACYCHOICE_DANTZIG = 0,          /**< use Dantzig's rule (reduced cost) to base efficacy on */
   GCG_EFFICIACYCHOICE_STEEPESTEDGE = 1,     /**< use steepest edge rule s( to base efficacy on */
   GCG_EFFICIACYCHOICE_LAMBDA = 2            /**< use lambda pricing to base efficacy on */
};
typedef enum GCG_Efficiacychoice GCG_EFFICIACYCHOICE;

typedef struct GCG_PriceStore GCG_PRICESTORE;     /**< storage for priced variables */

#ifdef __cplusplus
}
#endif

#endif
