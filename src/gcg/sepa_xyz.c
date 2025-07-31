/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    sepa_xyz.c
 *
 * @brief   xyz separator for master problem (put your description here)
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sepa_xyz.h"
#include "gcg/gcg.h"


#define SEPA_NAME              "xyz"
#define SEPA_DESC              "separator template"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */


/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** SCIP separator data */
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


/** copy method for separator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_SEPACOPY(sepaCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaCopyXyz NULL
#endif


/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_SEPAFREE(sepaFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeXyz NULL
#endif


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitXyz NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitXyz NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolXyz NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolXyz NULL
#endif


/** LP solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECLP(sepaExeclpXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExeclpXyz NULL
#endif


/** arbitrary primal solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExecsolXyz NULL
#endif


/*
 * Callback methods of MASTER separator
 */

/* TODO: Implement separator methods. */

/** method for adding new master variable to cut */
static
GCG_DECL_SEPAGETCOLCOEFFICIENT(sepagetcolcoefficientXyz)
{  /*lint --e{715}*/
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** method for adding new master variable to cut */
static
GCG_DECL_SEPAGETVARCOEFFICIENT(sepagetvarcoefficientXyz)
{  /*lint --e{715}*/
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** method for adapting pricing objectives to consider cut */
static
GCG_DECL_SEPASETOBJECTIVE(sepasetobjectiveXyz)
{  /*lint --e{715}*/
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** method for modifying the outdated values of a gcg column */
static
GCG_DECL_SEPAADJUSTCOL(sepaadjustcolXyz)
{  /*lint --e{715}*/
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** callback to delete the sepamastercutdata */
static
GCG_DECL_SEPAMASTERCUTDELETE(sepamastercutdeleteXyz)
{  /*lint --e{715}*/
   assert(strcmp(SEPA_NAME, SCIPsepaGetName(GCGsepaGetScipSeparator(sepa))) == 0);
   SCIPerrorMessage("method of xyz separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/*
 * separator specific interface methods
 */

/** creates the xyz separator and includes it in SCIP */
SCIP_RETCODE GCGincludeSepaXyz(
   GCG*                 gcg                /**< GCG data structure */
)
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;
   GCG_SEPA* gcgsepa;
   SCIP* scip;

   scip = GCGgetMasterprob(gcg);

   /* create xyz separator data */
   sepadata = NULL;
   /* TODO: (optional) create separator specific data here */

   sepa = NULL;

   SCIP_CALL( GCGrelaxIncludeSepa(gcg, &sepa, &gcgsepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
      SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpXyz, sepaExecsolXyz, sepadata, sepaadjustcolXyz, sepagetcolcoefficientXyz,
      sepagetvarcoefficientXyz, sepasetobjectiveXyz, sepamastercutdeleteXyz) );

   assert(sepa != NULL);
   assert(gcgsepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyXyz) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeXyz) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitXyz) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitXyz) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolXyz) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolXyz) );

   /* add xyz separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}