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

#include <assert.h>

#include "cons_decomp.h"
#include "dec_arrowheur.h"
#include "dec_borderheur.h"

//#include "dec_cutpacking.h"
//#include "dec_stairexact.h"
#include "dec_stairheur.h"
//#include "reader_dec.h"
#include "reader_gp.h"
#include "relax_gcg.h"

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
   SCIP_ARROWHEURDATA* arrowheurdata;
//   SCIP_CUTPACKINGDATA* cutpackingdata;
   SCIP_BORDERHEURDATA* borderheurdata;

   SCIP_STAIRHEURDATA *stairheurdata;
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
   SCIP_VAR** transvar2origvar;

   assert(decdecomp != NULL);
   assert(decdecomp->linkingconss != NULL);
   assert(decdecomp->nsubscipvars != NULL);
   assert(decdecomp->subscipvars != NULL);

   origvars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);

   SCIP_CALL(SCIPallocBufferArray(scip, &transvar2origvar, nvars));
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
      int index;
      transvar = SCIPvarGetTransVar(origvars[i]);
      assert(transvar != NULL);
      index = SCIPvarGetProbindex(transvar);
      if(index >= 0)
         transvar2origvar[index] = origvars[i];
   }

   for( i = 0; i < decdecomp->nblocks; ++i)
   {
      assert(decdecomp->subscipvars[i] != NULL);
      for( j = 0; j < decdecomp->nsubscipvars[i]; ++j)
      {
         int index;
         assert(decdecomp->subscipvars[i][j] != NULL);
         index = SCIPvarGetProbindex(decdecomp->subscipvars[i][j]);
         assert(index >= 0);
         assert(index < nvars);
         //SCIP_CALL(GCGrelaxSetOriginalVarBlockNr(scip, decdecomp->subscipvars[i][j], i+1));
         SCIP_CALL(GCGrelaxSetOriginalVarBlockNr(scip, transvar2origvar[index], i));
      }
   }
   SCIPfreeBufferArray(scip, &transvar2origvar);
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
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT result;
   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   SCIP_CALL(decdecompCreate(scip, &(conshdlrdata->decdecomp)));

//   SCIP_CALL(SCIPReaderDecSetDecomp(scip, conshdlrdata->decdecomp));
   SCIP_CALL(SCIPReaderGpSetDecomp(scip, conshdlrdata->decdecomp));
   SCIP_CALL(SCIPArrowHeurSetDecomp(scip, conshdlrdata->arrowheurdata, conshdlrdata->decdecomp));
   SCIP_CALL(SCIPBorderheurSetDecomp(scip, conshdlrdata->borderheurdata, conshdlrdata->decdecomp));
   //   SCIP_CALL(SCIPCutpackingSetDecomp(scip, conshdlrdata->cutpackingdata, conshdlrdata->decdecomp));

   if( GCGrelaxGetNPricingprobs(scip) <= 0 )
   {
//      SCIP_CALL(detectAndBuildBordered(scip, conshdlrdata->borderheurdata, &result));
      SCIP_CALL(detectAndBuildArrowHead(scip, conshdlrdata->arrowheurdata, &result));
//      SCIP_CALL(detectStructureCutpacking(scip, conshdlrdata->cutpackingdata, &result));
      //SCIP_CALL(detectStructureStairheur(scip, conshdlrdata->stairheurdata, &result));
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

   assert(conshdlr != NULL);
   assert(scip != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   freeArrowheurData(scip, &conshdlrdata->arrowheurdata);
   freeBorderheurData(scip, &conshdlrdata->borderheurdata);
   freeStairheurData(scip, &conshdlrdata->stairheurdata);
//   freeCutpackingData(scip, &conshdlrdata->cutpackingdata);
   decdecompFree(scip, &conshdlrdata->decdecomp);
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

   SCIP_CALL(createArrowheurData(scip, &conshdlrdata->arrowheurdata));
   SCIP_CALL(createBorderheurData(scip, &conshdlrdata->borderheurdata));
   SCIP_CALL(createStairheurData(scip, &conshdlrdata->stairheurdata));
//   SCIP_CALL(createCutpackingData(scip, &conshdlrdata->cutpackingdata));
   SCIP_CALL(SCIPincludeDetectionArrowheur(scip, conshdlrdata->arrowheurdata));
   SCIP_CALL(SCIPincludeDetectionStairheur(scip, conshdlrdata->stairheurdata));
   SCIP_CALL(SCIPincludeDetectionBorderheur(scip, conshdlrdata->borderheurdata));

//   SCIP_CALL(SCIPincludeDetectionCutpacking(scip, conshdlrdata->cutpackingdata));

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
