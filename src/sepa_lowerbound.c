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
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepa_lowerbound.h"
#include "probdata_gcg.h"

#define SEPA_NAME              "lowerbound"
#define SEPA_DESC              "separator template"
#define SEPA_PRIORITY           1000000
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */




/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** separator data */
struct SCIP_SepaData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of separator
 */

/* TODO: Implement all necessary separator methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_SEPAFREE(sepaFreeLowerbound)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lowerbound separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeLowerbound NULL
#endif


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitLowerbound)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lowerbound separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitLowerbound NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitLowerbound)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lowerbound separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitLowerbound NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolLowerbound)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lowerbound separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolLowerbound NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolLowerbound)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of lowerbound separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolLowerbound NULL
#endif


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
   /* TODO: (optional) create separator specific data here */

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, SEPA_DELAY,
         sepaFreeLowerbound, sepaInitLowerbound, sepaExitLowerbound, 
         sepaInitsolLowerbound, sepaExitsolLowerbound,
         sepaExeclpLowerbound, sepaExecsolLowerbound,
         sepadata) );

   /* add lowerbound separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
