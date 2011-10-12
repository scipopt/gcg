/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_master.h
 * @brief  master separator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_MASTER_H__
#define __SCIP_SEPA_MASTER_H__


#include "scip/scip.h"


/** creates the master separator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeSepaMaster(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* returns the array of original cuts saved in the separator data */
extern
SCIP_ROW** GCGsepaGetOrigcuts(
   SCIP*                 scip
   );

/* returns the number of original cuts saved in the separator data */
extern
int GCGsepaGetNOrigcuts(
   SCIP*                 scip
   );

/* returns the array of master cuts saved in the separator data */
extern
SCIP_ROW** GCGsepaGetMastercuts(
   SCIP*                 scip
   );

/* returns the number of master cuts saved in the separator data */
extern
int GCGsepaGetNMastercuts(
   SCIP*                 scip
   );

#endif
