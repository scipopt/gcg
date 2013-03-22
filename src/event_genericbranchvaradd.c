/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_genericbranchvaradd.c
 * @ingroup EVENTS
 * @brief  eventhdlr for Genericbranchvaradd event
 * @author Marcel Schmickerath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG

#include <assert.h>
#include <string.h>

#include "event_genericbranchvaradd.h"
#include "branch_generic.h"
#include "relax_gcg.h"
#include "cons_masterbranch.h"
#include "pricer_gcg.h"
#include "scip/cons_linear.h"
#include "type_branchgcg.h"
#include "pub_gcgvar.h"

#define EVENTHDLR_NAME         "genericbranchvaradd"
#define EVENTHDLR_DESC         "event handler for adding a new generated mastervar into the right branching constraints by using Vanderbecks generic branching scheme"


/** event handler data */
struct SCIP_EventhdlrData
{
};

/*
 * Local methods
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
#define eventCopyGenericbranchvaradd NULL

/** destructor of event handler to free user data (called when SCIP is exiting) */
#define eventFreeGenericbranchvaradd NULL

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#define eventInitsolGenericbranchvaradd NULL

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#define eventExitsolGenericbranchvaradd NULL

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitGenericbranchvaradd)
{  /*lint --e{715}*/
	assert(scip != NULL);
	assert(eventhdlr != NULL);
	assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

	/* notify SCIP that your event handler wants to react on the event type */
	SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

	return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitGenericbranchvaradd)
{  /*lint --e{715}*/
	assert(scip != NULL);
	assert(eventhdlr != NULL);
	assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

	/* notify SCIP that your event handler wants to drop the event type */
	SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );

	return SCIP_OKAY;
}

/** frees specific event data */
#define eventDeleteGenericbranchvaradd NULL

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecGenericbranchvaradd)
{  /*lint --e{715}*/
   SCIP* origscip;
   SCIP_CONS* masterbranchcons;
   SCIP_CONS* parentcons;
   SCIP_Bool varinS;
   SCIP_VAR* mastervar;
   SCIP_VAR** allorigvars;
   SCIP_VAR** mastervars;
   GCG_BRANCHDATA* branchdata;
   int p;
   int allnorigvars;
   int nmastervars;

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   varinS = TRUE;
   p = 0;
   mastervar = SCIPeventGetVar(event);

   origscip = GCGpricerGetOrigprob(scip);
   assert(origscip != NULL);

   //SCIPdebugMessage("exec method of event_genericbranchvaradd\n");

   masterbranchcons = GCGconsMasterbranchGetActiveCons(scip);
   assert(masterbranchcons != NULL);
   SCIP_CALL( SCIPgetVarsData(origscip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetVarsData(scip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

   parentcons = masterbranchcons;
   branchdata = GCGconsMasterbranchGetBranchdata(parentcons);

   if( GCGvarIsMaster(mastervar) )
   {
      while( parentcons != NULL && branchdata != NULL && ( strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetbranchrule(parentcons)), "generic") == 0 || strcmp(SCIPbranchruleGetName(GCGconsMasterbranchGetOrigbranchrule(parentcons)), "generic") == 0 )
            && GCGbranchGenericBranchdataGetConsS(branchdata) != NULL && GCGbranchGenericBranchdataGetConsSsize(branchdata) > 0 )
      {
         SCIP_Bool blockfound;
         SCIP_VAR** pricingvars;
         int k;

         assert(branchdata != NULL);

         varinS = TRUE;

         if( (GCGbranchGenericBranchdataGetConsblocknr(branchdata) != GCGvarGetBlock(mastervar) && GCGvarGetBlock(mastervar) != -1)
               || (GCGvarGetBlock(mastervar) == -1 && !GCGvarIsLinking(mastervar)) )
         {
            parentcons = GCGconsMasterbranchGetParentcons(parentcons);

            if(parentcons != NULL)
               branchdata = GCGconsMasterbranchGetBranchdata(parentcons);

            continue;
         }

         blockfound = TRUE;

         if( GCGvarGetBlock(mastervar) == -1 )
         {
            assert( GCGvarIsLinking(mastervar) );
            blockfound = FALSE;

            pricingvars = GCGlinkingVarGetPricingVars(mastervar);
            assert(pricingvars != NULL );

            for( k=0; k<GCGlinkingVarGetNBlocks(mastervar); ++k )
            {
               if( pricingvars[k] != NULL )
               {
                  if( GCGvarGetBlock(pricingvars[k]) == GCGbranchGenericBranchdataGetConsblocknr(branchdata) )
                  {
                     blockfound = TRUE;
                     break;
                  }
               }
            }
         }
         if( !blockfound )
         {
            parentcons = GCGconsMasterbranchGetParentcons(parentcons);

            if(parentcons != NULL)
               branchdata = GCGconsMasterbranchGetBranchdata(parentcons);

            continue;
         }


        // SCIPdebugMessage("consSsize = %d\n", GCGbranchGenericBranchdataGetConsSsize(branchdata));

         for( p = 0; p < GCGbranchGenericBranchdataGetConsSsize(branchdata); ++p )
         {
            SCIP_Real generatorentry;

            generatorentry = getGeneratorEntry(scip, mastervar, GCGbranchGenericBranchdataGetConsS(branchdata)[p].component);

            if( GCGbranchGenericBranchdataGetConsS(branchdata)[p].sense == GCG_COMPSENSE_GE )
            {
               if( SCIPisLT(scip, generatorentry, GCGbranchGenericBranchdataGetConsS(branchdata)[p].bound) )
               {
                  varinS = FALSE;
                  break;
               }
            }
            else
            {
               if( SCIPisGE(scip, generatorentry, GCGbranchGenericBranchdataGetConsS(branchdata)[p].bound) )
               {
                  varinS = FALSE;
                  break;
               }
            }
         }
         if( varinS )
         {
            SCIPdebugMessage("mastervar is added\n");
            SCIP_CALL( SCIPaddCoefLinear(scip, GCGbranchGenericBranchdataGetMastercons(branchdata), mastervar, 1.0) );
         }

         parentcons = GCGconsMasterbranchGetParentcons(parentcons);
         branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
      }
   }

   return SCIP_OKAY;
}

/** includes event handler for best solution found */
SCIP_RETCODE SCIPincludeEventHdlrGenericbranchvaradd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   eventhdlrdata = NULL;

   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventCopyGenericbranchvaradd, eventFreeGenericbranchvaradd, eventInitGenericbranchvaradd, eventExitGenericbranchvaradd,
         eventInitsolGenericbranchvaradd, eventExitsolGenericbranchvaradd, eventDeleteGenericbranchvaradd, eventExecGenericbranchvaradd,
         eventhdlrdata) );

   return SCIP_OKAY;
}
