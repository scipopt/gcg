/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    branch_compbnd.c
 *
 * @brief   branching rule based on vanderbeck's component bound branching
 * @author  Til Mohr
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "branch_compbnd.h"
#include "type_branchgcg.h"


#define BRANCHRULE_NAME            "compbnd"                      /**< name of branching rule */
#define BRANCHRULE_DESC            "component bound branching"    /**< short description of branching rule */
#define BRANCHRULE_PRIORITY        0                              /**< priority of this branching rule */
#define BRANCHRULE_MAXDEPTH        -1                             /**< maximal depth level of the branching rule */
#define BRANCHRULE_MAXBOUNDDIST    1.0                            /**< maximal relative distance from current node's
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
SCIP_DECL_BRANCHCOPY(branchCopyCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyCompBnd NULL
#endif

/** destructor of branching rule to free user data (called when SCIP is exiting) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHFREE(branchFreeCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchFreeCompBnd NULL
#endif


/** initialization method of branching rule (called after problem was transformed) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHINIT(branchInitCompBnd)
{  /*lint --e{715}*/

   /* inform relaxator of GCG about the branching rule */
   SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, branchActiveMasterOrig,
         branchDeactiveMasterOrig, branchPropMasterOrig, branchMasterSolvedOrig, branchDataDeleteOrig) );

   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitCompBnd NULL
#endif


/** deinitialization method of branching rule (called before transformed problem is freed) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXIT(branchExitCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitCompBnd NULL
#endif


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHINITSOL(branchInitsolCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchInitsolCompBnd NULL
#endif


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXITSOL(branchExitsolCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExitsolCompBnd NULL
#endif


/** branching execution method for fractional LP solutions */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXECLP(branchExeclpCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExeclpCompBnd NULL
#endif


/** branching execution method for external candidates */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXECEXT(branchExecextCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecextCompBnd NULL
#endif


/** branching execution method for not completely fixed pseudo solutions */
#ifdef SCIP_DISABLED_CODE
static
SCIP_DECL_BRANCHEXECPS(branchExecpsCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchExecpsCompBnd NULL
#endif

/*
 * GCG specific branching rule callbacks
 */

/** activation method for branchrule, called when a node in the master problem is activated,
 *  should perform changes to the current node's problem due to the branchdata
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchActiveMasterCompBnd NULL
#endif


/** deactivation method for branchrule, called when a node in the master problem is deactivated,
 *  should undo changes to the current node's problem due to the branchdata
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchDeactiveMasterCompBnd NULL
#endif

/** propagation method for branchrule, called when a node in the master problem is propagated,
 *  should perform propagation at the current node due to the branchdata
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchPropMasterCompBnd NULL
#endif

/** method for branchrule, called when the master LP is solved at one node,
 *  can store pseudocosts for the branching decisions
 */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHMASTERSOLVED(branchMasterSolvedCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchMasterSolvedCompBnd NULL
#endif

/** frees branching data of an origbranch constraint (called when the origbranch constraint is deleted) */
#ifdef SCIP_DISABLED_CODE
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteCompBnd)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of compbnd branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchDataDeleteCompBnd NULL
#endif

/*
 * branching rule specific interface methods
 */

/** creates the compbnd branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleCompBnd(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create compbnd branching rule data */
   branchruledata = NULL;
   /* TODO: (optional) create branching rule specific data here */

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
         BRANCHRULE_MAXBOUNDDIST,
         branchCopyCompBnd, branchFreeCompBnd, branchInitCompBnd, branchExitCompBnd, branchInitsolCompBnd, branchExitsolCompBnd,
         branchExeclpCompBnd, branchExecextCompBnd, branchExecpsCompBnd,
         branchruledata) );

   /* add compbnd branching rule parameters */
   /* TODO: (optional) add branching rule specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
