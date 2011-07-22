/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_decomp.c,v 1.54.2.1 2011/01/02 11:19:45 bzfheinz Exp $"

/**@file   cons_decomp.c
 * @ingroup CONSHDLRS 
 * @brief  constraint handler for decomp constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>

#include "cons_decomp.h"

//#include "dec_cutpacking.h"
//#include "dec_stairexact.h"
#include "dec_stairheur.h"
//#include "reader_dec.h"
#include "reader_gp.h"
#include "reader_ref.h"
#include "relax_gcg.h"
#include "struct_detector.h"
#include "struct_decomp.h"
#include "string.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "decomp"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS         0 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA         TRUE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP         TRUE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL       TRUE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_DETECTION             1 /**< Which detection scheme should be used as default 0 = arrowhead, 1 = bordered, 2 = staircase */


/*
 * Data structures
 */

/* TODO: fill in the necessary constraint data */

/** constraint data for decomp constraints */
struct SCIP_ConsData
{
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   DECDECOMP* decdecomp;

   SCIP_STAIRHEURDATA *stairheurdata;
   DEC_DETECTOR** detectors;
   int *priorities;
   int ndetectors;
   int usedetection;
};




/*
 * Local methods
 */

/** converts the structure to the gcg format by setting the appropriate blocks and master constraints */
static
SCIP_RETCODE DECOMPconvertStructToGCG(
      SCIP*         scip,     /**< SCIP data structure          */
      DECDECOMP*    decdecomp /**< decdecom data structure      */
   )
{
   int i;
   int j;
   int nvars;
   SCIP_VAR** origvars;
   SCIP_HASHMAP* transvar2origvar;

   assert(decdecomp != NULL);
   assert(decdecomp->linkingconss != NULL);
   assert(decdecomp->nsubscipvars != NULL);
   assert(decdecomp->subscipvars != NULL);

   origvars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);

   SCIP_CALL(SCIPhashmapCreate(&transvar2origvar, SCIPblkmem(scip), nvars));
   GCGrelaxSetNPricingprobs(scip, decdecomp->nblocks);
   SCIP_CALL( GCGrelaxCreateOrigVarsData(scip) );

   /* set master constraints */
   for(i = 0; i < decdecomp->nlinkingconss; ++i)
   {
      assert(decdecomp->linkingconss[i] != NULL);
      SCIP_CALL(GCGrelaxMarkConsMaster(scip, decdecomp->linkingconss[i]));
   }

   /* prepare the map from transformed to original variables */
   for(i = 0; i < nvars; ++i)
   {
      SCIP_VAR* transvar;
      SCIP_CALL(SCIPgetTransformedVar(scip, origvars[i], &transvar));
      assert(transvar != NULL);
      SCIPhashmapInsert(transvar2origvar, transvar, origvars[i]);
   }

   for( i = 0; i < decdecomp->nblocks; ++i)
   {
      assert(decdecomp->subscipvars[i] != NULL);
      for( j = 0; j < decdecomp->nsubscipvars[i]; ++j)
      {
         assert(decdecomp->subscipvars[i][j] != NULL);
         if(SCIPvarIsActive(decdecomp->subscipvars[i][j]) && !SCIPvarIsDeleted(decdecomp->subscipvars[i][j]))
         {
            if(SCIPhashmapGetImage(transvar2origvar, decdecomp->subscipvars[i][j]) != NULL)
               SCIP_CALL(GCGrelaxSetOriginalVarBlockNr(scip, SCIPhashmapGetImage(transvar2origvar, decdecomp->subscipvars[i][j]), i));
         }
      }
   }
   SCIPhashmapFree(&transvar2origvar);
   return SCIP_OKAY;
}


/* put your local methods here, and declare them static */
static
SCIP_RETCODE decdecompCreate(
   SCIP* scip,           /**< Pointer to the SCIP instance */
   DECDECOMP** decomp /**< Pointer to the decdecomp instance */
   )
{
   assert(scip != NULL);
   assert(decomp != NULL);
   SCIP_CALL( SCIPallocMemory(scip, decomp) );

   (*decomp)->type = DEC_UNKNOWN;
   (*decomp)->constoblock = NULL;
   (*decomp)->vartoblock = NULL;
   (*decomp)->subscipvars = NULL;
   (*decomp)->subscipconss = NULL;
   (*decomp)->nsubscipconss = NULL;
   (*decomp)->nsubscipvars = NULL;
   (*decomp)->linkingconss = NULL;
   (*decomp)->nlinkingconss = 0;
   (*decomp)->linkingvars = NULL;
   (*decomp)->nlinkingvars = 0;
   (*decomp)->nblocks = 0;
   (*decomp)->consindex = NULL;
   (*decomp)->varindex = NULL;

   return SCIP_OKAY;
}

static
void decdecompFree(
   SCIP* scip,           /**< Pointer to the SCIP instance */
   DECDECOMP** decdecomp /**< Pointer to the decdecomp instance */
   )
{
   DECDECOMP* decomp;
   int i;

   assert( scip!= NULL );
   assert( decdecomp != NULL);
   decomp = *decdecomp;

   assert(decomp != NULL);

   for( i = 0; i < decomp->nblocks; ++i)
   {
      SCIPfreeMemoryArray(scip, &decomp->subscipvars[i]);
      SCIPfreeMemoryArray(scip, &decomp->subscipconss[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipvars);
   SCIPfreeMemoryArrayNull(scip, &decomp->subscipconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->nsubscipconss);

   /* free hashmaps if they are not NULL */
   if(decomp->constoblock != NULL)
      SCIPhashmapFree(&decomp->constoblock);
   if(decomp->vartoblock != NULL)
      SCIPhashmapFree(&decomp->vartoblock);
   if(decomp->varindex != NULL)
      SCIPhashmapFree(&decomp->varindex);
   if(decomp->consindex != NULL)
      SCIPhashmapFree(&decomp->consindex);

   SCIPfreeMemoryArrayNull(scip, &decomp->linkingconss);
   SCIPfreeMemoryArrayNull(scip, &decomp->linkingvars);

   SCIPfreeMemory(scip, decdecomp);
}

/*
 * Callback methods of constraint handler
 */

#define conshdlrCopyDecomp NULL
#define consFreeDecomp NULL
#define consInitDecomp NULL
#define consExitDecomp NULL
#define consInitpreDecomp NULL
#define consExitpreDecomp NULL
#define consDeleteDecomp NULL
#define consTransDecomp NULL
#define consInitlpDecomp NULL
#define consSepalpDecomp NULL
#define consSepasolDecomp NULL
#define consPropDecomp NULL
#define consPresolDecomp NULL
#define consRespropDecomp NULL
#define consActiveDecomp NULL
#define consDeactiveDecomp NULL
#define consEnableDecomp NULL
#define consDisableDecomp NULL
#define consPrintDecomp NULL
#define consCopyDecomp NULL
#define consParseDecomp NULL

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolDecomp)
{
   SCIP_READER* reader;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT result;
   int i;
   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   SCIP_CALL(decdecompCreate(scip, &(conshdlrdata->decdecomp)));

   reader = SCIPfindReader(scip, "refreader");
   if(reader != NULL)
   {
      SCIP_CALL(SCIPReaderREFSetDecomp(scip, reader, conshdlrdata->decdecomp));
   }
//   SCIP_CALL(SCIPReaderDecSetDecomp(scip, conshdlrdata->decdecomp));
   SCIP_CALL(SCIPReaderGpSetDecomp(scip, conshdlrdata->decdecomp));
//   SCIP_CALL(SCIPArrowHeurSetDecomp(scip, conshdlrdata->arrowheurdata, conshdlrdata->decdecomp));
//   SCIP_CALL(SCIPBorderheurSetDecomp(scip, conshdlrdata->borderheurdata, conshdlrdata->decdecomp));
   //   SCIP_CALL(SCIPCutpackingSetDecomp(scip, conshdlrdata->cutpackingdata, conshdlrdata->decdecomp));

   if( GCGrelaxGetNPricingprobs(scip) <= 0 )
   {
      SCIPdebugMessage("Trying %d detectors.\n", conshdlrdata->ndetectors);
      for(i = 0; i < conshdlrdata->ndetectors; ++i)
      {
         DEC_DETECTOR *detector;
         detector = conshdlrdata->detectors[i];
         assert(detector != NULL);

         if(detector->initDetection != NULL)
         {
            SCIPdebugMessage("Calling initDetection of %s\n", detector->name);
            SCIP_CALL((*detector->initDetection)(scip));
         }

         (*detector->setStructDecomp)(scip, conshdlrdata->decdecomp);
         SCIPdebugMessage("Calling detectStructure of %s: ", detector->name);
         SCIP_CALL((*detector->detectStructure)(scip, detector->decdata, &result));
         if(result == SCIP_SUCCESS)
         {
            SCIPdebugPrintf("Success!\n");
            break;
         }
         else
         {
            SCIPdebugPrintf("Failure!\n");
         }
      }

      SCIP_CALL(DECOMPconvertStructToGCG(scip, conshdlrdata->decdecomp));
   }

//   SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s_%d_dec.lp", SCIPgetProbName(scip), conshdlrdata->decdecomp->nblocks);
//   SCIP_CALL(SCIPwriteOrigProblem(scip, "prob_dec.lp", "lp", FALSE));

   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */

static
SCIP_DECL_CONSEXITSOL(consExitsolDecomp)
{  /*lint --e{715}*/

   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   for(i = 0; i < conshdlrdata->ndetectors; ++i)
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if(detector->exitDetection != NULL)
      {
         SCIPdebugMessage("Calling exitDetection of %s\n", detector->name);
         SCIP_CALL((*detector->exitDetection)(scip));
      }
   }
   decdecompFree(scip, &conshdlrdata->decdecomp);
   SCIPfreeMemoryArray(scip, &conshdlrdata->priorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->detectors);
   SCIPfreeBlockMemory(scip, &conshdlrdata);
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDecomp)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDecomp)
{
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckDecomp)
{

   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}



/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDecomp)
{
   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for decomp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create decomp constraint handler data */
   SCIP_CALL(SCIPallocBlockMemory(scip, &conshdlrdata));
   assert(conshdlrdata != NULL);

   conshdlrdata->decdecomp = NULL;
   conshdlrdata->ndetectors = 0;
   conshdlrdata->priorities = NULL;
   conshdlrdata->detectors = NULL;

   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         conshdlrCopyDecomp,
         consFreeDecomp, consInitDecomp, consExitDecomp,
         consInitpreDecomp, consExitpreDecomp, consInitsolDecomp, consExitsolDecomp,
         consDeleteDecomp, consTransDecomp, consInitlpDecomp,
         consSepalpDecomp, consSepasolDecomp, consEnfolpDecomp, consEnfopsDecomp, consCheckDecomp,
         consPropDecomp, consPresolDecomp, consRespropDecomp, consLockDecomp,
         consActiveDecomp, consDeactiveDecomp,
         consEnableDecomp, consDisableDecomp,
         consPrintDecomp, consCopyDecomp, consParseDecomp,
         conshdlrdata) );

   /* add decomp constraint handler parameters */
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */
   SCIP_CALL(SCIPaddIntParam(scip, "cons/decomp/usedetection", "Which detection scheme should be used 0 = arrowhead (default), 1 = bordered, 2 = staircase.\n", &conshdlrdata->usedetection, FALSE, DEFAULT_DETECTION, 0, 2, NULL, NULL ));
   return SCIP_OKAY;
}

/** creates and captures a decomp constraint */
SCIP_RETCODE SCIPcreateConsDecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name               /**< name of constraint */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsDecomp() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   SCIPerrorMessage("method of decomp constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527} --e{715}*/

   /* find the decomp constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("decomp constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   /* TODO: create and store constraint specific data here */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         FALSE, FALSE, FALSE, TRUE, FALSE) );

   return SCIP_OKAY;
}

/** returns the decomposition structure **/
extern
DECDECOMP* SCIPconshdlrDecompGetDecdecomp(
      SCIP *scip                            /**< SCIP data structure */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL)
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   return conshdlrdata->decdecomp;
}

extern
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*  detector
   )
{
   assert(detector != NULL);
   return detector->decdata;

}

extern
const char* DECdetectorGetName(
   DEC_DETECTOR*  detector
   )
{
   assert(detector != NULL);
   return detector->name;
}


extern
DEC_DETECTOR* DECfindDetector(
   SCIP *      scip,
   const char* name
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL)
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for(i = 0; i < conshdlrdata->ndetectors; ++i)
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if(strcmp(detector->name, name) == 0)
      {
         return detector;
      }
   }

   return NULL;
}

extern
SCIP_RETCODE DECincludeDetector(
   SCIP* scip,
   const char *name,
   int priority,
   DEC_DETECTORDATA *detectordata,
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)),
   DEC_DECL_SETSTRUCTDECOMP((*setStructDecomp)),
   DEC_DECL_INITDETECTOR((*initDetector)),
   DEC_DECL_EXITDETECTOR((*exitDetector))
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DETECTOR *detector;
   assert(scip != NULL);
   assert(name != NULL);
   assert(detectStructure != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL)
      return SCIP_ERROR;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL(SCIPallocBlockMemory(scip, &detector));
   assert(detector != NULL);

   SCIPdebugMessage("Adding detector %i: %s\n", conshdlrdata->ndetectors+1, name);


#ifndef NDEBUG
   assert(DECfindDetector(scip, name) == NULL);
#endif

   detector->decdata = detectordata;
   detector->name = name;
   detector->priority = priority;
   detector->detectStructure = detectStructure;
   detector->initDetection = initDetector;
   detector->setStructDecomp = setStructDecomp;
   detector->exitDetection = exitDetector;

   SCIP_CALL(SCIPreallocMemoryArray(scip, &conshdlrdata->detectors, conshdlrdata->ndetectors+1));
   SCIP_CALL(SCIPreallocMemoryArray(scip, &conshdlrdata->priorities, conshdlrdata->ndetectors+1));

   conshdlrdata->detectors[conshdlrdata->ndetectors] = detector;
   conshdlrdata->priorities[conshdlrdata->ndetectors] = priority;
   conshdlrdata->ndetectors = conshdlrdata->ndetectors+1;

   SCIPdebugMessage("Sorting %i detectors\n", conshdlrdata->ndetectors);
   SCIPsortIntPtr(conshdlrdata->priorities, (void**)conshdlrdata->detectors, conshdlrdata->ndetectors);

   return SCIP_OKAY;

}
