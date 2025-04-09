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

/**@file   bendersplugins.c
 * @brief  SCIP plugins for the master problem running in Benders' decomposition mode
 * @author Stephen J. Maher
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/bendersplugins.h"
#include "scip/scipdefplugins.h"
#include "scip/debug.h"

/* GCG specific stuff */
#include "gcg/dialog_master.h"
#include "gcg/disp_master.h"

/** includes default GCG Benders' decomposition plugins */
SCIP_RETCODE GCGincludeBendersPlugins(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* masterprob = GCGgetBendersMasterprob(gcg);
   /* including the dialog used for the master problem */
   SCIP_CALL( GCGincludeDialogMaster(gcg) );
   /* including the SCIP default plugins */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(masterprob) ); /* nonlinear must be before linear, quadratic, abspower, and and due to constraint upgrading */
   SCIP_CALL( SCIPincludeNlhdlrQuadratic(masterprob) ); /* quadratic must be before linear due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrLinear(masterprob) ); /* linear must be before its specializations due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrAnd(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrBenders(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrBenderslp(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrBounddisjunction(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrCardinality(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrConjunction(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrCountsols(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrCumulative(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrDisjunction(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrIndicator(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrIntegral(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrKnapsack(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrLinking(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrLogicor(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrOr(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrOrbisack(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrOrbitope(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrPseudoboolean(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrSetppc(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrSOS1(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrSOS2(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrSuperindicator(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrSymresack(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrVarbound(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrXor(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrComponents(masterprob) );
   SCIP_CALL( SCIPincludeReaderBnd(masterprob) );
   SCIP_CALL( SCIPincludeReaderCcg(masterprob) );
   SCIP_CALL( SCIPincludeReaderCip(masterprob) );
   SCIP_CALL( SCIPincludeReaderCnf(masterprob) );
   SCIP_CALL( SCIPincludeReaderCor(masterprob) );
   SCIP_CALL( SCIPincludeReaderDiff(masterprob) );
   SCIP_CALL( SCIPincludeReaderFix(masterprob) );
   SCIP_CALL( SCIPincludeReaderFzn(masterprob) );
   SCIP_CALL( SCIPincludeReaderGms(masterprob) );
   SCIP_CALL( SCIPincludeReaderLp(masterprob) );
   SCIP_CALL( SCIPincludeReaderMps(masterprob) );
   SCIP_CALL( SCIPincludeReaderMst(masterprob) );
   SCIP_CALL( SCIPincludeReaderOpb(masterprob) );
   SCIP_CALL( SCIPincludeReaderOsil(masterprob) );
   SCIP_CALL( SCIPincludeReaderPip(masterprob) );
   SCIP_CALL( SCIPincludeReaderPpm(masterprob) );
   SCIP_CALL( SCIPincludeReaderPbm(masterprob) );
   SCIP_CALL( SCIPincludeReaderRlp(masterprob) );
   SCIP_CALL( SCIPincludeReaderSmps(masterprob) );
   SCIP_CALL( SCIPincludeReaderSol(masterprob) );
   SCIP_CALL( SCIPincludeReaderSto(masterprob) );
   SCIP_CALL( SCIPincludeReaderTim(masterprob) );
   SCIP_CALL( SCIPincludeReaderWbo(masterprob) );
   SCIP_CALL( SCIPincludeReaderZpl(masterprob) );
   SCIP_CALL( SCIPincludePresolBoundshift(masterprob) );
   SCIP_CALL( SCIPincludePresolConvertinttobin(masterprob) );
   SCIP_CALL( SCIPincludePresolDomcol(masterprob) );
   SCIP_CALL( SCIPincludePresolDualagg(masterprob) );
   SCIP_CALL( SCIPincludePresolDualcomp(masterprob) );
   SCIP_CALL( SCIPincludePresolDualinfer(masterprob) );
   SCIP_CALL( SCIPincludePresolGateextraction(masterprob) );
   SCIP_CALL( SCIPincludePresolImplics(masterprob) );
   SCIP_CALL( SCIPincludePresolInttobinary(masterprob) );
   SCIP_CALL( SCIPincludePresolQPKKTref(masterprob) );
   SCIP_CALL( SCIPincludePresolRedvub(masterprob) );
   SCIP_CALL( SCIPincludePresolTrivial(masterprob) );
   SCIP_CALL( SCIPincludePresolTworowbnd(masterprob) );
   SCIP_CALL( SCIPincludePresolSparsify(masterprob) );
   SCIP_CALL( SCIPincludePresolStuffing(masterprob) );
   SCIP_CALL( SCIPincludeNodeselBfs(masterprob) );
   SCIP_CALL( SCIPincludeNodeselBreadthfirst(masterprob) );
   SCIP_CALL( SCIPincludeNodeselDfs(masterprob) );
   SCIP_CALL( SCIPincludeNodeselEstimate(masterprob) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(masterprob) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(masterprob) );
   SCIP_CALL( SCIPincludeNodeselUct(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleAllfullstrong(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleCloud(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleDistribution(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleFullstrong(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleInference(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleLeastinf(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleMostinf(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleMultAggr(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleNodereopt(masterprob) );
   SCIP_CALL( SCIPincludeBranchrulePscost(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleRandom(masterprob) );
   SCIP_CALL( SCIPincludeBranchruleRelpscost(masterprob) );
   SCIP_CALL( SCIPincludeEventHdlrSolvingphase(masterprob) );
   SCIP_CALL( SCIPincludeComprLargestrepr(masterprob) );
   SCIP_CALL( SCIPincludeComprWeakcompr(masterprob) );
   SCIP_CALL( SCIPincludeHeurActconsdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurBound(masterprob) );
   SCIP_CALL( SCIPincludeHeurClique(masterprob) );
   SCIP_CALL( SCIPincludeHeurCoefdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurCompletesol(masterprob) );
   SCIP_CALL( SCIPincludeHeurCrossover(masterprob) );
   SCIP_CALL( SCIPincludeHeurDins(masterprob) );
   SCIP_CALL( SCIPincludeHeurDistributiondiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurDualval(masterprob) );
   SCIP_CALL( SCIPincludeHeurFeaspump(masterprob) );
   SCIP_CALL( SCIPincludeHeurFixandinfer(masterprob) );
   SCIP_CALL( SCIPincludeHeurFracdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurGins(masterprob) );
   SCIP_CALL( SCIPincludeHeurGuideddiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurZeroobj(masterprob) );
   SCIP_CALL( SCIPincludeHeurIndicator(masterprob) );
   SCIP_CALL( SCIPincludeHeurIntdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurIntshifting(masterprob) );
   SCIP_CALL( SCIPincludeHeurLinesearchdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurLocalbranching(masterprob) );
   SCIP_CALL( SCIPincludeHeurLocks(masterprob) );
   SCIP_CALL( SCIPincludeHeurLpface(masterprob) );
   SCIP_CALL( SCIPincludeHeurAlns(masterprob) );
   SCIP_CALL( SCIPincludeHeurNlpdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurMutation(masterprob) );
   SCIP_CALL( SCIPincludeHeurMultistart(masterprob) );
   SCIP_CALL( SCIPincludeHeurMpec(masterprob) );
   SCIP_CALL( SCIPincludeHeurObjpscostdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurOctane(masterprob) );
   SCIP_CALL( SCIPincludeHeurOfins(masterprob) );
   SCIP_CALL( SCIPincludeHeurOneopt(masterprob) );
   SCIP_CALL( SCIPincludeHeurProximity(masterprob) );
   SCIP_CALL( SCIPincludeHeurPscostdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurRandrounding(masterprob) );
   SCIP_CALL( SCIPincludeHeurRens(masterprob) );
   SCIP_CALL( SCIPincludeHeurReoptsols(masterprob) );
   SCIP_CALL( SCIPincludeHeurRepair(masterprob) );
   SCIP_CALL( SCIPincludeHeurRins(masterprob) );
   SCIP_CALL( SCIPincludeHeurRootsoldiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurRounding(masterprob) );
   SCIP_CALL( SCIPincludeHeurShiftandpropagate(masterprob) );
   SCIP_CALL( SCIPincludeHeurShifting(masterprob) );
   SCIP_CALL( SCIPincludeHeurSimplerounding(masterprob) );
   SCIP_CALL( SCIPincludeHeurSubNlp(masterprob) );
   SCIP_CALL( SCIPincludeHeurTrivial(masterprob) );
   SCIP_CALL( SCIPincludeHeurTrivialnegation(masterprob) );
   SCIP_CALL( SCIPincludeHeurTrySol(masterprob) );
   SCIP_CALL( SCIPincludeHeurTwoopt(masterprob) );
   SCIP_CALL( SCIPincludeHeurUndercover(masterprob) );
   SCIP_CALL( SCIPincludeHeurVbounds(masterprob) );
   SCIP_CALL( SCIPincludeHeurVeclendiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurZirounding(masterprob) );
   SCIP_CALL( SCIPincludePropDualfix(masterprob) );
   SCIP_CALL( SCIPincludePropGenvbounds(masterprob) );
   SCIP_CALL( SCIPincludePropObbt(masterprob) );
   SCIP_CALL( SCIPincludePropNlobbt(masterprob) );
   SCIP_CALL( SCIPincludePropProbing(masterprob) );
   SCIP_CALL( SCIPincludePropPseudoobj(masterprob) );
   SCIP_CALL( SCIPincludePropRedcost(masterprob) );
   SCIP_CALL( SCIPincludePropRootredcost(masterprob) );
   SCIP_CALL( SCIPincludePropVbounds(masterprob) );
   SCIP_CALL( SCIPincludeSepaCGMIP(masterprob) );
   SCIP_CALL( SCIPincludeSepaClique(masterprob) );
   SCIP_CALL( SCIPincludeSepaClosecuts(masterprob) );
   SCIP_CALL( SCIPincludeSepaAggregation(masterprob) );
   SCIP_CALL( SCIPincludeSepaConvexproj(masterprob) );
   SCIP_CALL( SCIPincludeSepaDisjunctive(masterprob) );
   SCIP_CALL( SCIPincludeSepaEccuts(masterprob) );
   SCIP_CALL( SCIPincludeSepaGauge(masterprob) );
   SCIP_CALL( SCIPincludeSepaGomory(masterprob) );
   SCIP_CALL( SCIPincludeSepaImpliedbounds(masterprob) );
   SCIP_CALL( SCIPincludeSepaIntobj(masterprob) );
   SCIP_CALL( SCIPincludeSepaMcf(masterprob) );
   SCIP_CALL( SCIPincludeSepaOddcycle(masterprob) );
   SCIP_CALL( SCIPincludeSepaRapidlearning(masterprob) );
   SCIP_CALL( SCIPincludeSepaZerohalf(masterprob) );
   SCIP_CALL( SCIPincludeEventHdlrSofttimelimit(masterprob) );
   SCIP_CALL( SCIPincludeConcurrentScipSolvers(masterprob) );
   SCIP_CALL( SCIPincludeBendersDefault(masterprob) );

   SCIP_CALL( SCIPdebugIncludeProp(masterprob) ); /*lint !e506 !e774*/

   /* including the display used for the master problem */
   SCIP_CALL( GCGincludeDispMaster(gcg, masterprob) );
   SCIP_CALL( SCIPincludeTableDefault(masterprob) );

   /* Benders' decomposition specific plugins */


   return SCIP_OKAY;
}
