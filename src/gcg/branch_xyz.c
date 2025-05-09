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

/**@file    branch_xyz.c
 *
 * @brief   xyz branching rule (put your description here)
 * @author  Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/branch_xyz.h"
#include "gcg/type_branchgcg.h"


#define BRANCHRULE_NAME        "xyz branching rule"           /**< name of branching rule */
#define BRANCHRULE_DESC        "branching rule template"      /**< short description of branching rule */
#define BRANCHRULE_PRIORITY        0                          /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH        -1                         /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST    1.0                        /**< maximal relative distance from current node's
                                                                   dual bound to primal bound compared to best node's
                                                                   dual bound for applying branching */


/*
 * Data structures
 */

/* TODO: fill in the necessary branching rule data */

/** branching rule data */
struct SCIP_BranchruleData
{
};

/* TODO: fill in the necessary branching data */
struct GCG_BranchData
{
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/* TODO: Implement all necessary branching rule methods. The methods with an #ifdef SCIP_DISABLED_CODE ... #else #define ... are optional */


/** copy method for branchrule plugins (called when SCIP copies plugins) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHCOPY(branchCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyXyz NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHFREE(branchFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchFreeXyz NULL
#endif


/** deinitialization method of branching rule (called before transformed problem is freed) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXIT(branchExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitXyz NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHINITSOL(branchInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolXyz NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolXyz NULL
#endif


/** branching execution method for fractional LP solutions */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXECLP(branchExeclpXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExeclpXyz NULL
#endif


/** branching execution method for external candidates */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXECEXT(branchExecextXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextXyz NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXECPS(branchExecpsXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsXyz NULL
#endif

/*
 * GCG specific branching rule callbacks
 */

/** activation method for branchrule, called when a node in the master problem is activated,
 *  should perform changes to the current node's problem due to the branchdata
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchActiveMasterXyz NULL
#endif


/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchDeactiveMasterXyz NULL
#endif

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchPropMasterXyz NULL
#endif

/** method for branchrule, called when the master LP is solved at one node,
 *  can store pseudocosts for the branching decisions
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHMASTERSOLVED(branchMasterSolvedXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchMasterSolvedXyz NULL
#endif

/** frees branching data of an origbranch constraint (called when the origbranch constraint is deleted) */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchDataDeleteXyz NULL
#endif

/** called when a new mastervariable has been created and added to the master */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHNEWCOL(branchNewColXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchNewColXyz NULL
#endif

/** call to get a possible extended master cons generated by this branching rule */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHGETEXTENDEDMASTERCONS(branchGetExtendedMasterConsXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchGetExtendedmasterconsXyz NULL
#endif

/** call to determine the coefficient of a column solution in the extended master cons */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHGETEXTENDEDMASTERCONSCOEFF(branchGetExtendedmasterconsCoeffXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchGetExtendedmasterconsCoeffXyz NULL
#endif


/** initialization method of branching rule (called after problem was transformed) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHINIT(branchInitXyz)
{  /*lint --e{715}*/
   GCG_BRANCHRULE* gcgbranchrule = NULL;

   /* inform relaxator of GCG about the branching rule */
   SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, &gcgbranchrule, branchActiveMasterOrig,
         branchDeactiveMasterOrig, branchPropMasterOrig, branchMasterSolvedOrig, branchDataDeleteOrig,
         branchNewColXyz, branchGetExtendedMasterConsXyz, branchGetExtendedmasterconsCoeff) );

   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitXyz NULL
#endif


/*
 * branching rule specific interface methods
 */

/** creates the xyz branching rule and includes it in SCIP */
SCIP_RETCODE GCGincludeBranchruleXyz(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create xyz branching rule data */
   branchruledata = NULL;
   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST,
         branchCopyXyz, branchFreeXyz, branchInitXyz, branchExitXyz, branchInitsolXyz, branchExitsolXyz,
         branchExeclpXyz, branchExecextXyz, branchExecpsXyz,
         branchruledata) );

   /* add xyz branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
