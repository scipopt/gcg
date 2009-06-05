/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   sepa_lowerbound.c
 * @ingroup SEPARATORS
 * @brief  lowerbound separator
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepa_lowerbound.h"
#include "probdata_gcg.h"

#define SEPA_NAME              "lowerbound"
#define SEPA_DESC              "separator for cutting off nodes due to the lower bound"
#define SEPA_PRIORITY           1000000
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */




/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
};





/*
 * Callback methods of separator
 */

/* define not used callbacks as NULL */
#define sepaFreeLowerbound NULL
#define sepaInitLowerbound NULL
#define sepaExitLowerbound NULL
#define sepaInitsolLowerbound NULL
#define sepaExitsolLowerbound NULL


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpLowerbound)
{  
   //printf("sepaexeclp\n");

   *result = SCIP_DIDNOTFIND;

   if ( SCIPisObjIntegral(GCGprobGetOrigprob(scip)) && SCIPisFeasLT(scip, SCIPgetPrimalbound(scip), SCIPgetLocalDualbound(scip)+0.999) )
      *result = SCIP_CUTOFF;

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolLowerbound)
{  
   //printf("sepaexecsol\n");

   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}



/*
 * separator specific interface methods
 */

/** creates the lowerbound separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create lowerbound separator data */
   sepadata = NULL;

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeLowerbound, sepaInitLowerbound, sepaExitLowerbound, 
         sepaInitsolLowerbound, sepaExitsolLowerbound,
         sepaExeclpLowerbound, sepaExecsolLowerbound,
         sepadata) );

   return SCIP_OKAY;
}
