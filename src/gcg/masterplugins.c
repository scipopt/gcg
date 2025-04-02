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

/**@file   masterplugins.c
 * @brief  SCIP plugins for generic column generation
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <scip/def.h>
#define USEHEURS 1
#define USESEPA 0
#define USEPROP 1

#include "gcg/masterplugins.h"

#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_conjunction.h"
#include "scip/cutsel_hybrid.h"
#include "scip/cons_integral.h"
#include "scip/cons_indicator.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_or.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"

#if USEHEURS
#include "scip/heur_actconsdiving.h"
#include "scip/heur_adaptivediving.h"
#include "scip/heur_bound.h"
#include "scip/heur_clique.h"
#include "scip/heur_coefdiving.h"
#include "scip/heur_completesol.h"
#include "scip/heur_conflictdiving.h"
#include "scip/heur_crossover.h"
#include "scip/heur_dins.h"
#include "scip/heur_distributiondiving.h"
#include "scip/heur_dks.h"
#include "scip/heur_dps.h"
#include "scip/heur_dualval.h"
#include "scip/heur_farkasdiving.h"
#include "scip/heur_feaspump.h"
#include "scip/heur_fixandinfer.h"
#include "scip/heur_fracdiving.h"
#include "scip/heur_gins.h"
#include "scip/heur_guideddiving.h"
#include "scip/heur_indicator.h"
#include "scip/heur_indicatordiving.h"
#include "scip/heur_intdiving.h"
#include "scip/heur_intshifting.h"
#include "scip/heur_linesearchdiving.h"
#include "scip/heur_localbranching.h"
#include "scip/heur_locks.h"
#include "scip/heur_lpface.h"
#include "scip/heur_alns.h"
#include "scip/heur_multistart.h"
#include "scip/heur_mutation.h"
#include "scip/heur_mpec.h"
#include "scip/heur_nlpdiving.h"
#include "scip/heur_objpscostdiving.h"
#include "scip/heur_octane.h"
#include "scip/heur_ofins.h"
#include "scip/heur_oneopt.h"
#include "scip/heur_padm.h"
#include "scip/heur_pscostdiving.h"
#include "scip/heur_proximity.h"
#include "scip/heur_randrounding.h"
#include "scip/heur_rens.h"
#include "scip/heur_reoptsols.h"
#include "scip/heur_repair.h"
#include "scip/heur_rins.h"
#include "scip/heur_rootsoldiving.h"
#include "scip/heur_rounding.h"
#include "scip/heur_scheduler.h"
#include "scip/heur_shiftandpropagate.h"
#include "scip/heur_shifting.h"
#include "scip/heur_simplerounding.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trivial.h"
#include "scip/heur_trivialnegation.h"
#include "scip/heur_trustregion.h"
#include "scip/heur_trysol.h"
#include "scip/heur_twoopt.h"
#include "scip/heur_undercover.h"
#include "scip/heur_vbounds.h"
#include "scip/heur_veclendiving.h"
#include "scip/heur_zeroobj.h"
#include "scip/heur_zirounding.h"
#endif

#include "scip/presol_implics.h"
#include "scip/presol_inttobinary.h"
#include "gcg/presol_roundbound.h"
#include "scip/presol_boundshift.h"

#if USEPROP
#include "scip/prop_dualfix.h"
#include "scip/prop_genvbounds.h"
#include "scip/prop_probing.h"
#include "scip/prop_pseudoobj.h"
#include "scip/prop_rootredcost.h"
#include "scip/prop_redcost.h"
#include "scip/prop_vbounds.h"
#endif

#if USESEPA
#include "scip/sepa_clique.h"
#include "scip/sepa_cmir.h"
#include "scip/sepa_flowcover.h"
#include "scip/sepa_gomory.h"
#include "scip/sepa_impliedbounds.h"
#include "scip/sepa_intobj.h"
#include "scip/sepa_mcf.h"
#include "scip/sepa_oddcycle.h"
#include "scip/sepa_strongcg.h"
#endif

/* Jonas' stuff */
#include "gcg/sepa_basis.h"

#include "scip/reader_cip.h"
#include "scip/reader_lp.h"
#include "scip/scipshell.h"

/* GCG specific stuff */
#include "gcg/pricer_gcg.h"
#include "gcg/nodesel_master.h"
#include "gcg/cons_masterbranch.h"
#include "gcg/cons_integralorig.h"
#include "gcg/sepa_original.h"
#include "gcg/branch_ryanfoster.h"
#include "gcg/branch_orig.h"
#include "gcg/branch_relpsprob.h"
#include "gcg/branch_generic.h"
#include "gcg/branch_bpstrong.h"
#include "gcg/branch_compbnd.h"
#include "scip/debug.h"
#include "gcg/dialog_master.h"
#include "gcg/disp_master.h"
#include "gcg/solver_knapsack.h"
#include "gcg/solver_mip.h"
#include "gcg/event_bestsol.h"
#include "gcg/event_relaxsol.h"
#include "gcg/event_solvingstats.h"
#include "gcg/event_display.h"

/* Erik's GCG solver */
#include "gcg/solver_gcg.h"

/* Christian's heuristics */
#include "gcg/heur_greedycolsel.h"
#include "gcg/heur_masterdiving.h"
#include "gcg/heur_mastercoefdiving.h"
#include "gcg/heur_masterfracdiving.h"
#include "gcg/heur_masterlinesdiving.h"
#include "gcg/heur_mastervecldiving.h"
#include "gcg/heur_relaxcolsel.h"
#include "gcg/heur_restmaster.h"
#include "gcg/heur_setcover.h"

#include "heur_ipcolgen.h"

#ifdef WITH_CLIQUER
#include "gcg/solver_cliquer.h"
#endif

#ifdef WITH_CPLEXSOLVER
#include "gcg/solver_cplex.h"
#endif

#ifdef WITH_HIGHS
#include "gcg/solver_highs.h"
#endif

#include "scip/table_default.h"

/** includes default GCG master plugins */
SCIP_RETCODE GCGincludeMasterPlugins(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP* masterprob = GCGgetDwMasterprob(gcg);
   SCIP_CALL( GCGincludeDialogMaster(gcg) );
   SCIP_CALL( SCIPincludeConshdlrLinear(masterprob) ); /* linear must be first due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrAnd(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrBounddisjunction(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrConjunction(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrIndicator(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrIntegral(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrKnapsack(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrLogicor(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrOr(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrSetppc(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrVarbound(masterprob) );
   SCIP_CALL( SCIPincludeConshdlrXor(masterprob) );

   SCIP_CALL( SCIPincludeReaderCip(masterprob) );
   SCIP_CALL( SCIPincludeReaderLp(masterprob) );

   SCIP_CALL( SCIPincludePresolBoundshift(masterprob) );
   SCIP_CALL( SCIPincludePresolImplics(masterprob) );
   SCIP_CALL( SCIPincludePresolInttobinary(masterprob) );
   SCIP_CALL( GCGincludePresolRoundbound(masterprob) );

#if USEPROP
   SCIP_CALL( SCIPincludePropDualfix(masterprob) );
   SCIP_CALL( SCIPincludePropGenvbounds(masterprob) );
   SCIP_CALL( SCIPincludePropProbing(masterprob) );
   SCIP_CALL( SCIPincludePropPseudoobj(masterprob) );
   SCIP_CALL( SCIPincludePropRootredcost(masterprob) );
   SCIP_CALL( SCIPincludePropRedcost(masterprob) );
   SCIP_CALL( SCIPincludePropVbounds(masterprob) );
#endif

   SCIP_CALL( GCGincludeNodeselMaster(gcg) );
   SCIP_CALL( GCGincludeConshdlrIntegralOrig(gcg) );
   SCIP_CALL( GCGincludeBranchruleRyanfoster(gcg) );
   SCIP_CALL( GCGincludeBranchruleOrig(gcg) );
   SCIP_CALL( GCGincludeBranchruleRelpsprob(gcg) );
   SCIP_CALL( GCGincludeBranchruleGeneric(gcg) );
   SCIP_CALL( GCGincludeBranchruleBPStrong(gcg) );
   SCIP_CALL( GCGincludeBranchruleCompBnd(gcg) );

#if USEHEURS
   SCIP_CALL( SCIPincludeHeurActconsdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurAdaptivediving(masterprob) );
   SCIP_CALL( SCIPincludeHeurBound(masterprob) );
   SCIP_CALL( SCIPincludeHeurClique(masterprob) );
   SCIP_CALL( SCIPincludeHeurCoefdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurCompletesol(masterprob) );
   SCIP_CALL( SCIPincludeHeurConflictdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurCrossover(masterprob) );
   SCIP_CALL( SCIPincludeHeurDins(masterprob) );
   SCIP_CALL( SCIPincludeHeurDistributiondiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurDKS(masterprob) );
   SCIP_CALL( SCIPincludeHeurDps(masterprob) );
   SCIP_CALL( SCIPincludeHeurDualval(masterprob) );
   SCIP_CALL( SCIPincludeHeurFarkasdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurFeaspump(masterprob) );
   SCIP_CALL( SCIPincludeHeurFixandinfer(masterprob) );
   SCIP_CALL( SCIPincludeHeurFracdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurGins(masterprob) );
   SCIP_CALL( SCIPincludeHeurGuideddiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurZeroobj(masterprob) );
   SCIP_CALL( SCIPincludeHeurIndicator(masterprob) );
   SCIP_CALL( SCIPincludeHeurIndicatordiving(masterprob) );
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
   SCIP_CALL( SCIPincludeHeurPADM(masterprob) );
   SCIP_CALL( SCIPincludeHeurProximity(masterprob) );
   SCIP_CALL( SCIPincludeHeurPscostdiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurRandrounding(masterprob) );
   SCIP_CALL( SCIPincludeHeurRens(masterprob) );
   SCIP_CALL( SCIPincludeHeurReoptsols(masterprob) );
   SCIP_CALL( SCIPincludeHeurRepair(masterprob) );
   SCIP_CALL( SCIPincludeHeurRins(masterprob) );
   SCIP_CALL( SCIPincludeHeurRootsoldiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurRounding(masterprob) );
   SCIP_CALL( SCIPincludeHeurScheduler(masterprob) );
   SCIP_CALL( SCIPincludeHeurShiftandpropagate(masterprob) );
   SCIP_CALL( SCIPincludeHeurShifting(masterprob) );
   SCIP_CALL( SCIPincludeHeurSubNlp(masterprob) );
   SCIP_CALL( SCIPincludeHeurTrivial(masterprob) );
   SCIP_CALL( SCIPincludeHeurTrivialnegation(masterprob) );
   SCIP_CALL( SCIPincludeHeurTrustregion(masterprob) );
   SCIP_CALL( SCIPincludeHeurTrySol(masterprob) );
   SCIP_CALL( SCIPincludeHeurTwoopt(masterprob) );
   SCIP_CALL( SCIPincludeHeurUndercover(masterprob) );
   SCIP_CALL( SCIPincludeHeurVbounds(masterprob) );
   SCIP_CALL( SCIPincludeHeurVeclendiving(masterprob) );
   SCIP_CALL( SCIPincludeHeurZirounding(masterprob) );

   SCIP_CALL( SCIPincludeHeurSimplerounding(masterprob) );

   /* Christian's heuristics */
   SCIP_CALL( GCGincludeHeurGreedycolsel(gcg) );
   SCIP_CALL( GCGincludeEventHdlrMasterdiving(gcg) );
   SCIP_CALL( GCGincludeHeurMastercoefdiving(gcg) );
   SCIP_CALL( GCGincludeHeurMasterfracdiving(gcg) );
   SCIP_CALL( GCGincludeHeurMasterlinesdiving(gcg) );
   SCIP_CALL( GCGincludeHeurMastervecldiving(gcg) );
   SCIP_CALL( GCGincludeHeurRelaxcolsel(gcg) );
   SCIP_CALL( GCGincludeHeurRestmaster(gcg) );
   SCIP_CALL( GCGincludeHeurSetcover(gcg) );

   SCIP_CALL( SCIPincludeHeurIPcolgen(gcg) );
#endif

#if USESEPA
   SCIP_CALL( SCIPincludeSepaClique(masterprob) );
   SCIP_CALL( SCIPincludeSepaCmir(masterprob) );
   SCIP_CALL( SCIPincludeSepaFlowcover(masterprob) );
   SCIP_CALL( SCIPincludeSepaGomory(masterprob) );
   SCIP_CALL( SCIPincludeSepaImpliedbounds(masterprob) );
   SCIP_CALL( SCIPincludeSepaIntobj(masterprob) );
   SCIP_CALL( SCIPincludeSepaMcf(masterprob) );
   SCIP_CALL( SCIPincludeSepaOddcycle(masterprob) );
   SCIP_CALL( SCIPincludeSepaRedcost(masterprob) );
   SCIP_CALL( SCIPincludeSepaZerohalf(masterprob) );
#endif
   SCIP_CALL( GCGincludeSepaOriginal(gcg) );
   SCIP_CALL( SCIPincludeCutselHybrid(masterprob) );
   SCIP_CALL( GCGincludeDispMaster(gcg, masterprob) );
   SCIP_CALL( SCIPdebugIncludeProp(masterprob) ); /*lint !e506 !e774*/
   SCIP_CALL( SCIPincludeTableDefault(masterprob) );

   /* Jonas' stuff */
   SCIP_CALL( GCGincludeSepaBasis(gcg) );

   SCIP_CALL( GCGincludeSolverKnapsack(gcg) );
   SCIP_CALL( GCGincludeSolverMip(gcg) );
   SCIP_CALL( GCGincludeSolverGcg(gcg) );

#ifdef WITH_CLIQUER
   SCIP_CALL( GCGincludeSolverCliquer(gcg) );
#endif

#ifdef WITH_CPLEXSOLVER
   SCIP_CALL( GCGincludeSolverCplex(gcg) );
#endif

#ifdef WITH_HIGHS
   SCIP_CALL( GCGincludeSolverHighs(gcg) );
#endif

   /* include masterbranch constraint handler */
   SCIP_CALL( GCGincludeConshdlrMasterbranch(gcg) );

   SCIP_CALL( GCGincludeEventHdlrBestsol(masterprob) );
   SCIP_CALL( GCGincludeEventHdlrRelaxsol(gcg) );
   SCIP_CALL( GCGincludeEventHdlrSolvingstats(gcg) );
   SCIP_CALL( GCGincludeEventHdlrDisplay(gcg) );

   return SCIP_OKAY;
}
