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

#include <assert.h>
#include <string.h>

#include "event_genericbranchvaradd.h"
#include "branch_generic.h"
//#include "branch_generic.c"
#include "relax_gcg.h"
#include "cons_masterbranch.h"
#include "pricer_gcg.h"
#include "scip/cons_linear.h"
#include "type_branchgcg.h"
#include "pub_gcgvar.h"

#define EVENTHDLR_NAME         "genericbranchvaradd"
#define EVENTHDLR_DESC         "event handler for adding a new generated mastervar into the right branching constraints by using Vanderbecks generic branching scheme"


typedef int ComponentBoundSequence[3];

/*
 * Data structures
 */
struct GCG_BranchData
{
	ComponentBoundSequence**   C;             /**< S[k] bound sequence for block k */ //!!! sort of each C[i]=S[i] is important !!!
	int*               sequencsizes;                 /**< number of bounds in S[k] */
	int                Csize;
	ComponentBoundSequence*   S;             /**< component bound sequence which induce the current branching constraint */
	int                Ssize;
	ComponentBoundSequence*   childS;       /**< component bound sequence which induce the child nodes, need for prune by dominance */
	int                childSsize;
	int                blocknr;             /**< number of block branching was performed */
	int                childnr;
	int                lhs;
	int                nchildNodes;
	int*               childlhs;
	//SCIP_Real*         gerenicPseudocostsOnOrigvars;  /**< giving the branchingpriorities */
	SCIP_CONS*         mastercons;          /**< constraint enforcing the branching restriction in the master problem */
};
/* TODO: fill in the necessary event handler data */

/** event handler data */
struct SCIP_EventhdlrData
{
};

/*
 * Local methods
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyGenericbranchvaradd) 
{  /*lint --e{715}*/
	assert(scip != NULL);
	assert(eventhdlr != NULL);
	assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

	/* call inclusion method of event handler */
	SCIP_CALL( SCIPincludeEventHdlrGenericbranchvaradd(scip) );

	return SCIP_OKAY;
}

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
	SCIP* masterscip;
	SCIP_CONS* masterbranchcons;
	SCIP_CONS* parentcons;
	SCIP_Bool varinS;
	SCIP_VAR* mastervar;
	SCIP_VAR** allorigvars;
	GCG_BRANCHDATA* branchdata; 
	int p;
	int j;
	int allnorigvars;

	assert(eventhdlr != NULL);
	assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
	assert(event != NULL);
	assert(scip != NULL);
	assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

	SCIPdebugMessage("exec method of event handler Genericbranchvaradd for adding a pricingvar\n");

	varinS = TRUE;
	p = 0;
	mastervar = SCIPeventGetVar(event); 
	masterscip = GCGrelaxGetMasterprob(scip);
	masterbranchcons = GCGconsMasterbranchGetActiveCons(scip);
	SCIP_CALL( SCIPgetVarsData(scip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );

	if(masterbranchcons != NULL && GCGvarIsMaster(mastervar))
	{
		parentcons = masterbranchcons;
		while( parentcons != NULL )
		{
			branchdata = GCGconsMasterbranchGetBranchdata(parentcons);

			if( branchdata->blocknr != GCGvarGetBlock(mastervar) )
			{
				parentcons = GCGconsMasterbranchGetParentcons(masterbranchcons);
				continue;
			}

			for( p=0; p<branchdata->Ssize; ++p)
			{
				SCIP_VAR** presentmastervars;
				SCIP_VAR** origvars;
				int npresentmastervars;
				int norigvars;
				SCIP_Real generator_i;
				SCIP_Bool present;

				present = FALSE;
				npresentmastervars = GCGoriginalVarGetNMastervars(allorigvars[branchdata->S[p][0]]);
				presentmastervars = GCGoriginalVarGetMastervars(allorigvars[branchdata->S[p][0]]);
				norigvars = GCGmasterVarGetNOrigvars(mastervar);  
				origvars = GCGmasterVarGetOrigvars(mastervar);

				for(j=0; j<npresentmastervars; ++j)
				{
					if(mastervar == presentmastervars[j])
					{
						present = TRUE;
					}
				}
				if(present)
				{
					for(j=0; j<norigvars; ++j)
					{
						if(allorigvars[branchdata->S[p][0]] == origvars[j]) 
						{
							generator_i = GCGmasterVarGetOrigvals(mastervar)[j]; // ??? bei discr gleich dem generator ???
									break;
						}
					}
				}
				else
					generator_i = 0;
				if(branchdata->S[p][1] == 1 )
				{
					//if(GCGmasterVarGetOrigvals(mastervar)[branchdata->S[p][0]] < branchdata->S[p][2])
						if( generator_i < branchdata->S[p][2])
						{
							varinS = FALSE;
							break;
						}
				}
				else
				{
					if(generator_i >= branchdata->S[p][2])
					{
						varinS = FALSE;
						break;
					}
				}
			}
			if( varinS )
				SCIP_CALL( SCIPaddCoefLinear(masterscip, branchdata->mastercons, mastervar, 1.0) );

			parentcons = GCGconsMasterbranchGetParentcons(masterbranchcons);
		}

	}
	else
		SCIPinfoMessage(scip, NULL, "no masterbranchcons or mastervar found in SCIP <%s>\n", SCIPgetProbName(scip) );

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
