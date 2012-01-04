/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_decomp.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 *
 * This constraint handler will run all registered structure detectors in
 * increasing priority until the first detector finds a suitable structure.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>

#include "cons_decomp.h"

#include "reader_gp.h"
#include "reader_ref.h"
#include "reader_dec.h"
#include "cons_connected.h"
#include "relax_gcg.h"
#include "struct_detector.h"
#include "struct_decomp.h"
#include "string.h"
#include "scip_misc.h"
#include "scip/clock.h"
#include "pub_gcgvar.h"
#include "pub_decomp.h"
/* constraint handler properties */
#define CONSHDLR_NAME          "decomp"
#define CONSHDLR_DESC          "constraint handler template"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
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

/** constraint data for decomp constraints */
struct SCIP_ConsData
{
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   DECDECOMP* decdecomp;
   DECDECOMP** decdecomps;
   DEC_DETECTOR** detectors;
   int *priorities;
   int ndetectors;
   int usedetection;
   SCIP_CLOCK* detectorclock;
   SCIP_Bool hasrun;
   int ndecomps;
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of constraint handler
 */

#define conshdlrCopyDecomp NULL
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
#define consDelvarsDecomp NULL
#define consPrintDecomp NULL
#define consCopyDecomp NULL
#define consParseDecomp NULL

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeDecomp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL(SCIPfreeClock(scip, &conshdlrdata->detectorclock));
   DECdecdecompFree(scip, &conshdlrdata->decdecomp);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if( detector->exitDetection != NULL)
      {
         SCIPdebugMessage("Calling exitDetection of %s\n", detector->name);
         SCIP_CALL((*detector->exitDetection)(scip, detector));
      }
   }

   SCIPfreeMemoryArray(scip, &conshdlrdata->priorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->detectors);
   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolDecomp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( !conshdlrdata->hasrun )
   {
      SCIP_CALL( DECdetectStructure(scip) );
      assert( conshdlrdata->hasrun );
   }
   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolDecomp)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->decdecomps != NULL )
   {
      int i;
      for( i = 0; i < conshdlrdata->ndecomps; ++i )
      {
         DECdecdecompFree(scip, &conshdlrdata->decdecomps[i]);
      }
      SCIPfreeMemoryArray(scip, &conshdlrdata->decdecomps);
   }
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
   SCIP_READER* reader;

   /* create decomp constraint handler data */
   SCIP_CALL(SCIPallocMemory(scip, &conshdlrdata));
   assert(conshdlrdata != NULL);

   conshdlrdata->decdecomp = NULL;
   conshdlrdata->decdecomps = NULL;
   conshdlrdata->ndecomps = 0;
   conshdlrdata->ndetectors = 0;
   conshdlrdata->priorities = NULL;
   conshdlrdata->detectors = NULL;
   conshdlrdata->hasrun = FALSE;

   SCIP_CALL(SCIPcreateWallClock(scip, &conshdlrdata->detectorclock));

   /* TODO: (optional) create constraint handler specific data here */

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         SCIP_PROPTIMING_AFTERLPNODE, conshdlrCopyDecomp,
         consFreeDecomp, consInitDecomp, consExitDecomp,
         consInitpreDecomp, consExitpreDecomp, consInitsolDecomp, consExitsolDecomp,
         consDeleteDecomp, consTransDecomp, consInitlpDecomp,
         consSepalpDecomp, consSepasolDecomp, consEnfolpDecomp, consEnfopsDecomp, consCheckDecomp,
         consPropDecomp, consPresolDecomp, consRespropDecomp, consLockDecomp,
         consActiveDecomp, consDeactiveDecomp,
         consEnableDecomp, consDisableDecomp,
         consDelvarsDecomp, consPrintDecomp, consCopyDecomp, consParseDecomp,
         conshdlrdata) );

   SCIP_CALL(DECdecdecompCreate(scip, &(conshdlrdata->decdecomp)));

   reader = SCIPfindReader(scip, "refreader");
   if( reader != NULL)
   {
      SCIP_CALL(SCIPReaderREFSetDecomp(scip, reader, conshdlrdata->decdecomp));
   }
   reader = SCIPfindReader(scip,"decreader");
   if(reader!=NULL)
   {
      SCIP_CALL(SCIPReaderDecSetDecomp(scip, conshdlrdata->decdecomp));
   }
   SCIP_CALL(SCIPReaderGpSetDecomp(scip, conshdlrdata->decdecomp));

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
      SCIP *scip                             /**< SCIP data structure */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->decdecomp != NULL);

   return conshdlrdata->decdecomp;
}

/** returns the decomposition structure **/
extern
DECDECOMP** SCIPconshdlrDecompGetDecdecomps(
      SCIP *scip                             /**< SCIP data structure */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->decdecomps;
}

/** returns the decomposition structure **/
extern
int SCIPconshdlrDecompGetNDecdecomps(
      SCIP *scip                             /**< SCIP data structure */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->ndecomps;
}

/** returns the data of the provided detector */
extern
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*  detector                   /**< detector data structure */
   )
{
   assert(detector != NULL);
   return detector->decdata;

}

/** returns the name of the provided detector */
extern
const char* DECdetectorGetName(
   DEC_DETECTOR*  detector
   )
{
   assert(detector != NULL);
   return detector->name;
}

/** Searches for the detector and returns it or returns NULL if detector is not found*/
extern
DEC_DETECTOR* DECfindDetector(
   SCIP *scip,                               /**< SCIP data structure */
   const char *name                          /**< Name of the detector to return */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for(i = 0; i < conshdlrdata->ndetectors; ++i)
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if( strcmp(detector->name, name) == 0 )
      {
         return detector;
      }
   }

   return NULL;
}

/** includes the detector */
extern
SCIP_RETCODE DECincludeDetector(
   SCIP* scip,                                     /**< SCIP data structure */
   const char *name,                               /**< name of the detector */
   DEC_DETECTORDATA *detectordata,                 /**< the associated detector data (or NULL) */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)),   /**< the method that will detect the structure (must not be NULL)*/
   DEC_DECL_SETSTRUCTDECOMP((*setStructDecomp)),   /**< interface method to tell detector where to store structure information (must not be NULL) */
   DEC_DECL_INITDETECTOR((*initDetector)),         /**< initialization method of detector (or NULL) */
   DEC_DECL_EXITDETECTOR((*exitDetector)),         /**< deinitialization method of detector (or NULL) */
   DEC_DECL_GETPRIORITY((*getPriority))            /**< interface method to get priority of detector (must not be NULL) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DETECTOR *detector;
   assert(scip != NULL);
   assert(name != NULL);
   assert(detectStructure != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   assert(detectStructure != NULL);
   assert(setStructDecomp != NULL);
   assert(getPriority != NULL);

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

   detector->detectStructure = detectStructure;

   detector->initDetection = initDetector;
   detector->setStructDecomp = setStructDecomp;
   detector->exitDetection = exitDetector;
   detector->getPriority = getPriority;
   detector->i = conshdlrdata->ndetectors;
   SCIP_CALL(SCIPreallocMemoryArray(scip, &conshdlrdata->detectors, conshdlrdata->ndetectors+1));
   SCIP_CALL(SCIPreallocMemoryArray(scip, &conshdlrdata->priorities, conshdlrdata->ndetectors+1));

   conshdlrdata->detectors[conshdlrdata->ndetectors] = detector;
   conshdlrdata->ndetectors = conshdlrdata->ndetectors+1;

   return SCIP_OKAY;

}

/** returns the remaning time of scip that the decomposition may use */
extern
SCIP_Real DECgetRemainingTime(
   SCIP* scip                    /**< SCIP data structure */
   )
{
   SCIP_Real timelimit;
   assert(scip != NULL);
   SCIP_CALL_ABORT(SCIPgetRealParam(scip, "limits/time", &timelimit));
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   return timelimit;
}


/** converts the structure to the gcg format by setting the appropriate blocks and master constraints */
extern
SCIP_RETCODE DECOMPconvertStructToGCG(
      SCIP*         scip,     /**< SCIP data structure          */
      DECDECOMP*    decdecomp /**< decdecom data structure      */
   )
{
   int i;
   int j;
   int k;
   int v;
   int nvars;
   SCIP_VAR** origvars;
   SCIP_HASHMAP* transvar2origvar;

   assert(decdecomp != NULL);
   assert(decdecomp->linkingconss != NULL || decdecomp->nlinkingconss == 0);
   assert(decdecomp->nsubscipvars != 0);
   assert(decdecomp->subscipvars != NULL);

   origvars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);

   SCIP_CALL(SCIPhashmapCreate(&transvar2origvar, SCIPblkmem(scip), nvars));
   GCGrelaxSetNPricingprobs(scip, decdecomp->nblocks);
   SCIP_CALL( GCGcreateOrigVarsData(scip) );

   /* set master constraints */
   for( i = 0; i < decdecomp->nlinkingconss; ++i )
   {
      assert(decdecomp->linkingconss[i] != NULL);
      SCIP_CALL(GCGrelaxMarkConsMaster(scip, decdecomp->linkingconss[i]));
   }

   /* prepare the map from transformed to original variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* transvar;
      SCIP_CALL(SCIPgetTransformedVar(scip, origvars[i], &transvar));
      assert(transvar != NULL);
      SCIPhashmapInsert(transvar2origvar, transvar, origvars[i]);
   }

   for( i = 0; i < decdecomp->nblocks; ++i )
   {
      assert(decdecomp->subscipvars[i] != NULL);
      for( j = 0; j < decdecomp->nsubscipvars[i]; ++j )
      {
         assert(decdecomp->subscipvars[i][j] != NULL);
         assert(isVarRelevant(decdecomp->subscipvars[i][j]));
         if(SCIPhashmapGetImage(transvar2origvar, decdecomp->subscipvars[i][j]) != NULL)
         {
            SCIP_CALL(GCGrelaxSetOriginalVarBlockNr(scip, SCIPhashmapGetImage(transvar2origvar, decdecomp->subscipvars[i][j]), i));
         }
         else
         {
            SCIP_CALL(GCGorigVarCreateData(scip, getRelevantVariable(decdecomp->subscipvars[i][j])));
            SCIP_CALL(GCGrelaxSetOriginalVarBlockNr(scip, getRelevantVariable(decdecomp->subscipvars[i][j]), i));
         }
      }
   }
   for( i = 0; i < decdecomp->nlinkingvars; ++i )
   {
      if( SCIPvarGetData(decdecomp->linkingvars[i]) == NULL)
      {
         int found;
         SCIP_CALL(GCGorigVarCreateData(scip, getRelevantVariable(decdecomp->linkingvars[i])));
         /* HACK; TODO: find out constraint blocks */
         for( j = 0; j < decdecomp->nblocks; ++j )
         {
            found = FALSE;
            for( k = 0; k < decdecomp->nsubscipconss[j]; ++k )
            {
               SCIP_VAR** curvars;
               int        ncurvars;
               curvars = SCIPgetVarsXXX(scip, decdecomp->subscipconss[j][k]);
               ncurvars = SCIPgetNVarsXXX(scip, decdecomp->subscipconss[j][k]);

               for( v = 0; v < ncurvars; ++v )
               {
                  if( SCIPvarGetProbvar(curvars[v]) == decdecomp->linkingvars[i] )
                  {
                     SCIPdebugMessage("%s is in %d\n", SCIPvarGetName(SCIPvarGetProbvar(curvars[v])), j);
                     SCIP_CALL(GCGrelaxSetOriginalVarBlockNr(scip, decdecomp->linkingvars[i], j));
                     found = TRUE;
                     break;
                  }
               }
               SCIPfreeMemoryArray(scip, &curvars);
               if( found )
                  break;
            }
         }
      }
   }

   SCIPhashmapFree(&transvar2origvar);
   return SCIP_OKAY;
}


/** interface method to detect the structure */
SCIP_RETCODE DECdetectStructure(
   SCIP *scip
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT result;
   int i;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists, cannot detect structure!\n");

      if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
         conshdlrdata->hasrun = TRUE;
      return SCIP_OKAY;
   }

   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      SCIP_CALL( SCIPtransformProb(scip) );

   SCIP_CALL(SCIPresetClock(scip, conshdlrdata->detectorclock));
   SCIP_CALL(SCIPstartClock(scip, conshdlrdata->detectorclock));

   if( GCGrelaxGetNPricingprobs(scip) <= 0 )
   {
      for( i = 0; i < conshdlrdata->ndetectors; ++i )
      {
         DEC_DETECTOR *detector;
         detector = conshdlrdata->detectors[i];
         assert(detector != NULL);
         conshdlrdata->priorities[i] = detector->getPriority(scip);
      }

      SCIPdebugMessage("Sorting %i detectors\n", conshdlrdata->ndetectors);
      SCIPsortIntPtr(conshdlrdata->priorities, (void**)conshdlrdata->detectors, conshdlrdata->ndetectors);

      SCIPdebugMessage("Trying %d detectors.\n", conshdlrdata->ndetectors);
      for( i = 0; i < conshdlrdata->ndetectors; ++i )
      {
         DEC_DETECTOR* detector;
         DECDECOMP** decdecomps;
         int ndecdecomps;

         ndecdecomps = -1;
         detector = conshdlrdata->detectors[i];
         assert(detector != NULL);

         if(detector->initDetection != NULL)
         {
            SCIPdebugMessage("Calling initDetection of %s\n", detector->name);
            SCIP_CALL((*detector->initDetection)(scip, detector));
         }

         (*detector->setStructDecomp)(scip, conshdlrdata->decdecomp);
         SCIPdebugMessage("Calling detectStructure of %s: ", detector->name);
         SCIP_CALL((*detector->detectStructure)(scip, detector->decdata, &decdecomps, &ndecdecomps,  &result));
         if( result == SCIP_SUCCESS )
         {
            assert(ndecdecomps >= 0);
            assert(decdecomps != NULL || ndecdecomps == 0);
            SCIPdebugPrintf("we have %d decompositions!\n", ndecdecomps);

            SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->decdecomps, conshdlrdata->ndecomps+ndecdecomps) );
            BMScopyMemoryArray(&conshdlrdata->decdecomps[conshdlrdata->ndecomps], decdecomps, ndecdecomps);
            SCIPfreeMemoryArray(scip, &decdecomps);
            conshdlrdata->ndecomps += ndecdecomps;

            break;
         }
         else
         {
            SCIPdebugPrintf("Failure!\n");
         }
      }

      SCIP_CALL(DECOMPconvertStructToGCG(scip, conshdlrdata->decdecomps[0]));
/*    for( i = 0; i < conshdlrdata->ndecomps; ++i)
      {
         DECdecdecompFree(scip, &(conshdlrdata->decdecomps[i]));
      }
*/
   }
   SCIP_CALL(SCIPstopClock(scip, conshdlrdata->detectorclock));
   SCIPdebugMessage("Detection took %fs\n", SCIPclockGetTime(conshdlrdata->detectorclock));

   conshdlrdata->hasrun = TRUE;
   return SCIP_OKAY;
}

/** write out all known decompositions **/
SCIP_RETCODE DECwriteAllDecomps(
   SCIP* scip,
   char* extension
   )
{
   int i;
   char name[SCIP_MAXSTRLEN];
   const char *pname;

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   assert(extension != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   pname = strrchr(SCIPgetProbName(scip), '/');
   if( pname == NULL )
      pname = SCIPgetProbName(scip);
 
   for ( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      FILE* file;
      SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d.%s", pname, DECdecdecompGetNBlocks(conshdlrdata->decdecomps[i]), extension);
      file = SCIPfopen(name, "w");
      assert(file != NULL);
      SCIP_CALL( SCIPwriteTransProblem(scip, name, extension, FALSE) );
   }

   return SCIP_OKAY;
}
