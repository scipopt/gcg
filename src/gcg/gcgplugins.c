/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   gcgplugins.c
 * @brief  SCIP plugins for generic column generation
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Michael Bastubbe
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/gcgplugins.h"
#include "scip/debug.h"

#define USEHEURS 1
#define USESEPA 1
#define USEPROP 1

/* include header files here, such that the user only has to include
 * gcgplugins.h
 */
#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_fixedvar.h"
#include "scip/cons_indicator.h"
#include "scip/cons_integral.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_linking.h"
#include "scip/cons_logicor.h"
#include "scip/cons_or.h"
#include "scip/cons_orbitope.h"
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

#include "scip/nodesel_bfs.h"
#include "scip/nodesel_breadthfirst.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/nodesel_uct.h"

#include "scip/event_solvingphase.h"

#include "scip/presol_boundshift.h"
#include "scip/presol_convertinttobin.h"
#include "scip/presol_domcol.h"
#include "scip/presol_dualagg.h"
#include "scip/presol_dualcomp.h"
#include "scip/presol_dualinfer.h"
#include "scip/presol_dualsparsify.h"
#include "scip/presol_implics.h"
#include "scip/presol_implint.h"
#include "scip/presol_inttobinary.h"
#include "scip/presol_redvub.h"
#include "scip/presol_sparsify.h"
#include "scip/presol_stuffing.h"
#include "scip/presol_trivial.h"
#include "scip/presol_tworowbnd.h"

#if USEPROP
#include "scip/prop_dualfix.h"
#include "scip/prop_genvbounds.h"
#include "scip/prop_nlobbt.h"
#include "scip/prop_obbt.h"
#include "scip/prop_probing.h"
#include "scip/prop_pseudoobj.h"
#include "scip/prop_rootredcost.h"
#include "scip/prop_redcost.h"
#include "scip/prop_symmetry.h"
#include "scip/prop_vbounds.h"
#endif

#include "scip/reader_bnd.h"
#include "scip/reader_ccg.h"
#include "scip/reader_cip.h"
#include "scip/reader_cnf.h"
#include "scip/reader_fix.h"
#include "scip/reader_fzn.h"
#include "scip/reader_gms.h"
#include "scip/reader_lp.h"
#include "scip/reader_mps.h"
#include "scip/reader_opb.h"
#include "scip/reader_osil.h"
#include "scip/reader_pip.h"
#include "scip/reader_pbm.h"
#include "scip/reader_rlp.h"
#include "scip/reader_sol.h"
#include "scip/reader_wbo.h"
#include "scip/reader_zpl.h"

#if USESEPA
#include "scip/sepa_eccuts.h"
#include "scip/sepa_cgmip.h"
#include "scip/sepa_clique.h"
#include "scip/sepa_closecuts.h"
#include "scip/sepa_aggregation.h"
#include "scip/sepa_convexproj.h"
#include "scip/sepa_disjunctive.h"
#include "scip/sepa_gomory.h"
#include "scip/sepa_impliedbounds.h"
#include "scip/sepa_interminor.h"
#include "scip/sepa_intobj.h"
#include "scip/sepa_lagromory.h"
#include "scip/sepa_mcf.h"
#include "scip/sepa_minor.h"
#include "scip/sepa_mixing.h"
#include "scip/sepa_oddcycle.h"
#include "scip/sepa_rapidlearning.h"
#include "scip/sepa_zerohalf.h"
#endif

#include "scip/cutsel_hybrid.h"

#include "scip/scipshell.h"
#include "gcg/reader_blk.h"
#include "gcg/reader_dec.h"
#include "gcg/pricer_gcg.h"
#include "gcg/relax_gcg.h"
#include "gcg/branch_empty.h"

#if WITH_JANSSON
#include "gcg/reader_jdec.h"
#endif

#include "gcg/cons_origbranch.h"
#include "gcg/disp_gcg.h"
#include "gcg/dialog_gcg.h"
#include "gcg/reader_ref.h"
#include "gcg/event_bestsol.h"
#include "gcg/event_mastersol.h"

/* visualization */
#include "gcg/reader_gp.h"
#include "gcg/reader_tex.h"
#include "gcg/reader_cls.h"

/* detection */
#include "gcg/cons_decomp.h"
#include "gcg/dec_constype.h"
#include "gcg/dec_consclass.h"
#include "gcg/dec_densemasterconss.h"
#include "gcg/dec_stairheur.h"
#include "gcg/dec_compgreedily.h"
#include "gcg/dec_staircase_lsp.h"
#include "gcg/dec_postprocess.h"
#include "gcg/dec_mastersetpack.h"
#include "gcg/dec_mastersetpart.h"
#include "gcg/dec_mastersetcover.h"
#include "gcg/dec_neighborhoodmaster.h"
#if WITH_HMETIS
#include "gcg/dec_hrcgpartition.h"
#include "gcg/dec_hrgpartition.h"
#include "gcg/dec_hcgpartition.h"
#endif
#include "gcg/dec_connectedbase.h"
#include "gcg/dec_connected_noNewLinkingVars.h"
#include "gcg/dec_generalmastersetcover.h"
#include "gcg/dec_generalmastersetpack.h"
#include "gcg/dec_generalmastersetpart.h"
#include "gcg/dec_staircase_lsp.h"
#include "gcg/dec_varclass.h"

#ifndef NO_AUT_LIB
#include "gcg/dec_isomorph.h"
#endif

/* Christian's heuristics */
#include "gcg/heur_gcgcoefdiving.h"
#include "gcg/heur_gcgdins.h"
#include "gcg/heur_gcgfeaspump.h"
#include "gcg/heur_gcgfracdiving.h"
#include "gcg/heur_gcgguideddiving.h"
#include "gcg/heur_gcglinesdiving.h"
#include "gcg/heur_gcgpscostdiving.h"
#include "gcg/heur_gcgrens.h"
#include "gcg/heur_gcgrins.h"
#include "gcg/heur_gcgrounding.h"
#include "gcg/heur_gcgshifting.h"
#include "gcg/heur_gcgsimplerounding.h"
#include "gcg/heur_gcgveclendiving.h"
#include "gcg/heur_gcgzirounding.h"
#include "gcg/heur_origdiving.h"
#include "gcg/heur_xpcrossover.h"
#include "gcg/heur_xprins.h"

/* Friedrike's detection stuff */
#include "gcg/scip_misc.h"
#include "scip/table_default.h"

/* Igor's detection with clustering */
#include "gcg/dec_dbscan.h"
#include "gcg/dec_mst.h"

/* classifiers */
#include "gcg/clscons_miplibconstypes.h"
#include "gcg/clscons_nnonzeros.h"
#include "gcg/clscons_consnamelevenshtein.h"
#include "gcg/clscons_consnamenonumbers.h"
#include "gcg/clscons_scipconstypes.h"
#include "gcg/clscons_gamsdomain.h"
#include "gcg/clscons_gamssymbol.h"

#include "gcg/clsvar_gamsdomain.h"
#include "gcg/clsvar_gamssymbol.h"
#include "gcg/clsvar_objvalues.h"
#include "gcg/clsvar_scipvartypes.h"
#include "gcg/clsvar_objvaluesigns.h"

/* scores */
#include "gcg/score_bender.h"
#include "gcg/score_border.h"
#include "gcg/score_classic.h"
#include "gcg/score_fawh.h"
#include "gcg/score_forswh.h"
#include "gcg/score_maxwhite.h"
#include "gcg/score_spfawh.h"
#include "gcg/score_spfwh.h"
#include "gcg/score_strong.h"


/** includes default plugins for generic column generation into SCIP */
SCIP_RETCODE GCGincludeGcgPlugins(
   GCG*                 gcg               /**< GCG data structure */
   )
{
   SCIP* scip = GCGgetOrigprob(gcg);
   SCIP_CALL( GCGincludeDialogGcg(gcg) );

   SCIP_CALL( SCIPincludeConshdlrAnd(scip) );
   SCIP_CALL( SCIPincludeConshdlrBounddisjunction(scip) );
   SCIP_CALL( SCIPincludeConshdlrFixedvar(scip) );
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) ); /* linear must be first due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrIndicator(scip) );
   SCIP_CALL( SCIPincludeConshdlrIntegral(scip) );
   SCIP_CALL( SCIPincludeConshdlrLinking(scip) );
   SCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );
   SCIP_CALL( SCIPincludeConshdlrLogicor(scip) );
   SCIP_CALL( SCIPincludeConshdlrOr(scip) );
   SCIP_CALL( SCIPincludeConshdlrOrbitope(scip) );
   SCIP_CALL( SCIPincludeConshdlrSetppc(scip) );
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );
   SCIP_CALL( SCIPincludeConshdlrXor(scip) );

   SCIP_CALL( SCIPincludeReaderBnd(scip) );
   SCIP_CALL( SCIPincludeReaderCcg(scip) );
   SCIP_CALL( SCIPincludeReaderCip(scip) );
   SCIP_CALL( SCIPincludeReaderCnf(scip) );
   SCIP_CALL( SCIPincludeReaderFix(scip) );
   SCIP_CALL( SCIPincludeReaderFzn(scip) );
   SCIP_CALL( SCIPincludeReaderGms(scip) );
   SCIP_CALL( SCIPincludeReaderLp(scip) );
   SCIP_CALL( SCIPincludeReaderMps(scip) );
   SCIP_CALL( SCIPincludeReaderOpb(scip) );
   SCIP_CALL( SCIPincludeReaderOsil(scip) );
   SCIP_CALL( SCIPincludeReaderPip(scip) );
   SCIP_CALL( SCIPincludeReaderPbm(scip) );
   SCIP_CALL( SCIPincludeReaderRlp(scip) );
   SCIP_CALL( SCIPincludeReaderSol(scip) );
   SCIP_CALL( SCIPincludeReaderWbo(scip) );
   SCIP_CALL( SCIPincludeReaderZpl(scip) );

   SCIP_CALL( SCIPincludePresolBoundshift(scip) );
   SCIP_CALL( SCIPincludePresolConvertinttobin(scip) );
   SCIP_CALL( SCIPincludePresolDomcol(scip) );
   SCIP_CALL( SCIPincludePresolDualagg(scip) );
   SCIP_CALL( SCIPincludePresolDualcomp(scip) );
   SCIP_CALL( SCIPincludePresolDualinfer(scip) );
   // SCIP_CALL( SCIPincludePresolGateextraction(scip) );  @todo: GCG cannot handle this presolver currently
   SCIP_CALL( SCIPincludePresolImplics(scip) );
   SCIP_CALL( SCIPincludePresolImplint(scip) );
   SCIP_CALL( SCIPincludePresolInttobinary(scip) );
   SCIP_CALL( SCIPincludePresolRedvub(scip) );
   SCIP_CALL( SCIPincludePresolTrivial(scip) );
   SCIP_CALL( SCIPincludePresolTworowbnd(scip) );
   SCIP_CALL( SCIPincludePresolSparsify(scip) );
   SCIP_CALL( SCIPincludePresolDualsparsify(scip) );
   SCIP_CALL( SCIPincludePresolStuffing(scip) );
   // @todo: Papilo

   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselBreadthfirst(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   SCIP_CALL( SCIPincludeNodeselEstimate(scip) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(scip) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );
   SCIP_CALL( SCIPincludeNodeselUct(scip) );

   SCIP_CALL( SCIPincludeEventHdlrSolvingphase(scip) );

#if USEPROP
   SCIP_CALL( SCIPincludePropDualfix(scip) );
   SCIP_CALL( SCIPincludePropGenvbounds(scip) );
   SCIP_CALL( SCIPincludePropObbt(scip) );
   SCIP_CALL( SCIPincludePropNlobbt(scip) );
   SCIP_CALL( SCIPincludePropProbing(scip) );
   SCIP_CALL( SCIPincludePropPseudoobj(scip) );
   SCIP_CALL( SCIPincludePropRedcost(scip) );
   SCIP_CALL( SCIPincludePropRootredcost(scip) );
   // SCIP_CALL( SCIPincludePropSymmetry(scip) );
   SCIP_CALL( SCIPincludePropVbounds(scip) );
#endif


#if USEHEURS
   SCIP_CALL( SCIPincludeHeurActconsdiving(scip) );
   SCIP_CALL( SCIPincludeHeurAdaptivediving(scip) );
   SCIP_CALL( SCIPincludeHeurBound(scip) );
   SCIP_CALL( SCIPincludeHeurClique(scip) );
   SCIP_CALL( SCIPincludeHeurCoefdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCompletesol(scip) );
   SCIP_CALL( SCIPincludeHeurConflictdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCrossover(scip) );
   SCIP_CALL( SCIPincludeHeurDins(scip) );
   SCIP_CALL( SCIPincludeHeurDistributiondiving(scip) );
   SCIP_CALL( SCIPincludeHeurDKS(scip) );
   SCIP_CALL( SCIPincludeHeurDps(scip) );
   SCIP_CALL( SCIPincludeHeurDualval(scip) );
   SCIP_CALL( SCIPincludeHeurFarkasdiving(scip) );
   SCIP_CALL( SCIPincludeHeurFeaspump(scip) );
   SCIP_CALL( SCIPincludeHeurFixandinfer(scip) );
   SCIP_CALL( SCIPincludeHeurFracdiving(scip) );
   SCIP_CALL( SCIPincludeHeurGins(scip) );
   SCIP_CALL( SCIPincludeHeurGuideddiving(scip) );
   SCIP_CALL( SCIPincludeHeurZeroobj(scip) );
   SCIP_CALL( SCIPincludeHeurIndicator(scip) );
   SCIP_CALL( SCIPincludeHeurIndicatordiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntdiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntshifting(scip) );
   SCIP_CALL( SCIPincludeHeurLinesearchdiving(scip) );
   SCIP_CALL( SCIPincludeHeurLocalbranching(scip) );
   SCIP_CALL( SCIPincludeHeurLocks(scip) );
   SCIP_CALL( SCIPincludeHeurLpface(scip) );
   SCIP_CALL( SCIPincludeHeurAlns(scip) );
   SCIP_CALL( SCIPincludeHeurNlpdiving(scip) );
   SCIP_CALL( SCIPincludeHeurMutation(scip) );
   SCIP_CALL( SCIPincludeHeurMultistart(scip) );
   SCIP_CALL( SCIPincludeHeurMpec(scip) );
   SCIP_CALL( SCIPincludeHeurObjpscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurOctane(scip) );
   SCIP_CALL( SCIPincludeHeurOfins(scip) );
   SCIP_CALL( SCIPincludeHeurOneopt(scip) );
   SCIP_CALL( SCIPincludeHeurPADM(scip) );
   SCIP_CALL( SCIPincludeHeurProximity(scip) );
   SCIP_CALL( SCIPincludeHeurPscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurRandrounding(scip) );
   SCIP_CALL( SCIPincludeHeurRens(scip) );
   SCIP_CALL( SCIPincludeHeurReoptsols(scip) );
   SCIP_CALL( SCIPincludeHeurRepair(scip) );
   SCIP_CALL( SCIPincludeHeurRins(scip) );
   SCIP_CALL( SCIPincludeHeurRootsoldiving(scip) );
   SCIP_CALL( SCIPincludeHeurRounding(scip) );
   SCIP_CALL( SCIPincludeHeurScheduler(scip) );
   SCIP_CALL( SCIPincludeHeurShiftandpropagate(scip) );
   SCIP_CALL( SCIPincludeHeurShifting(scip) );
   SCIP_CALL( SCIPincludeHeurSubNlp(scip) );
   SCIP_CALL( SCIPincludeHeurTrivial(scip) );
   SCIP_CALL( SCIPincludeHeurTrivialnegation(scip) );
   SCIP_CALL( SCIPincludeHeurTrustregion(scip) );
   SCIP_CALL( SCIPincludeHeurTrySol(scip) );
   SCIP_CALL( SCIPincludeHeurTwoopt(scip) );
   SCIP_CALL( SCIPincludeHeurUndercover(scip) );
   SCIP_CALL( SCIPincludeHeurVbounds(scip) );
   SCIP_CALL( SCIPincludeHeurVeclendiving(scip) );
   SCIP_CALL( SCIPincludeHeurZirounding(scip) );
#endif
   SCIP_CALL( SCIPincludeHeurSimplerounding(scip) );

#if USESEPA
   SCIP_CALL( SCIPincludeSepaCGMIP(scip) );
   SCIP_CALL( SCIPincludeSepaClique(scip) );
   SCIP_CALL( SCIPincludeSepaClosecuts(scip) );
   SCIP_CALL( SCIPincludeSepaAggregation(scip) );
   SCIP_CALL( SCIPincludeSepaConvexproj(scip) );
   SCIP_CALL( SCIPincludeSepaDisjunctive(scip) );
   SCIP_CALL( SCIPincludeSepaEccuts(scip) );
   SCIP_CALL( SCIPincludeSepaGomory(scip) );
   SCIP_CALL( SCIPincludeSepaImpliedbounds(scip) );
   SCIP_CALL( SCIPincludeSepaInterminor(scip) );
   SCIP_CALL( SCIPincludeSepaIntobj(scip) );
   SCIP_CALL( SCIPincludeSepaLagromory(scip) );
   SCIP_CALL( SCIPincludeSepaMcf(scip) );
   SCIP_CALL( SCIPincludeSepaMinor(scip) );
   SCIP_CALL( SCIPincludeSepaMixing(scip) );
   SCIP_CALL( SCIPincludeSepaOddcycle(scip) );
   SCIP_CALL( SCIPincludeSepaRapidlearning(scip) );
   SCIP_CALL( SCIPincludeSepaZerohalf(scip) );
#endif

   SCIP_CALL( SCIPincludeCutselHybrid(scip) );

   SCIP_CALL( GCGincludeRelaxGcg(gcg) );
   SCIP_CALL( GCGincludeReaderBlk(gcg) );
   SCIP_CALL( GCGincludeReaderDec(gcg) );
   SCIP_CALL( GCGincludeReaderRef(gcg) );
#if WITH_JANSSON
   SCIP_CALL( GCGincludeReaderJDec(gcg) );
#endif
   SCIP_CALL( GCGincludeBranchruleEmpty(gcg) );

   SCIP_CALL( GCGincludeConshdlrOrigbranch(gcg) );
   SCIP_CALL( GCGincludeEventHdlrBestsol(scip) );
   SCIP_CALL( GCGincludeEventHdlrMastersol(gcg) );

   /* Visualizations */
   SCIP_CALL( GCGincludeReaderGp(gcg) );
   SCIP_CALL( GCGincludeReaderTex(gcg) );
   SCIP_CALL( GCGincludeReaderCls(gcg) );

   /* Detectors and decompositions */
   SCIP_CALL( GCGincludeConshdlrDecomp(gcg) );
   SCIP_CALL( GCGincludeDetectorConstype(gcg) );
   SCIP_CALL( GCGincludeDetectorPostprocess(gcg) );
   SCIP_CALL( GCGincludeDetectorConsclass(gcg) );
   SCIP_CALL( GCGincludeDetectorDensemasterconss(gcg) );
   SCIP_CALL( GCGincludeDetectorNeighborhoodmaster(gcg) );
   SCIP_CALL( GCGincludeDetectorStairheur(gcg) );
   SCIP_CALL( GCGincludeDetectorStaircaseLsp(gcg) );
#ifdef WITH_GSL
   SCIP_CALL( GCGincludeDetectorDBSCAN(gcg) );
   SCIP_CALL( GCGincludeDetectorMST(gcg) );
#endif
   SCIP_CALL( GCGincludeDetectorCompgreedily(gcg) );
   SCIP_CALL( GCGincludeDetectorMastersetcover(gcg) );
   SCIP_CALL( GCGincludeDetectorMastersetpack(gcg) );
   SCIP_CALL( GCGincludeDetectorMastersetpart(gcg) );
#if WITH_HMETIS
   SCIP_CALL( GCGincludeDetectorHcgpartition(gcg) );
   SCIP_CALL( GCGincludeDetectorHrgpartition(gcg) );
   SCIP_CALL( GCGincludeDetectorHrcgpartition(gcg) );
#endif
   SCIP_CALL( GCGincludeDetectorConnectedbase(gcg) );
   SCIP_CALL( GCGincludeDetectorConnected_noNewLinkingVars(gcg) );
   SCIP_CALL( GCGincludeDetectorGeneralmastersetpack(gcg) );
   SCIP_CALL( GCGincludeDetectorGeneralmastersetpart(gcg) );
   SCIP_CALL( GCGincludeDetectorGeneralmastersetcover(gcg) );
   SCIP_CALL( GCGincludeDetectorVarclass(gcg) );

   #ifndef NO_AUT_LIB
      SCIP_CALL( GCGincludeDetectorIsomorphism(gcg) );
   #endif


   /* Classifiers */
   SCIP_CALL( GCGincludeConsClassifierNNonzeros(gcg) );
   SCIP_CALL( GCGincludeConsClassifierScipConstypes(gcg) );
   SCIP_CALL( GCGincludeConsClassifierMiplibConstypes(gcg) );
   SCIP_CALL( GCGincludeConsClassifierConsnameLevenshtein(gcg) );
   SCIP_CALL( GCGincludeConsClassifierForConsnamesDigitFreeIdentical(gcg) );
   SCIP_CALL( GCGincludeConsClassifierGamsdomain(gcg) );
   SCIP_CALL( GCGincludeConsClassifierGamssymbol(gcg) );

   SCIP_CALL( GCGincludeVarClassifierGamsdomain(gcg) );
   SCIP_CALL( GCGincludeVarClassifierGamssymbol(gcg) );
   SCIP_CALL( GCGincludeVarClassifierScipVartypes(gcg) );
   SCIP_CALL( GCGincludeVarClassifierObjValues(gcg) );
   SCIP_CALL( GCGincludeVarClassifierObjValueSigns(gcg) );


   /* Scores */
   SCIP_CALL( GCGincludeScoreBender(gcg) );
   SCIP_CALL( GCGincludeScoreBorder(gcg) );
   SCIP_CALL( GCGincludeScoreClassic(gcg) );
   SCIP_CALL( GCGincludeScoreFawh(gcg) );
   SCIP_CALL( GCGincludeScoreForswh(gcg) );
   SCIP_CALL( GCGincludeScoreMaxwhite(gcg) );
   SCIP_CALL( GCGincludeScoreSpfawh(gcg) );
   SCIP_CALL( GCGincludeScoreSpfwh(gcg) );
   SCIP_CALL( GCGincludeScoreStrongDecomp(gcg) );


   /* Christian's heuristics */
   SCIP_CALL( GCGincludeEventHdlrOrigdiving(gcg) );
   SCIP_CALL( GCGincludeHeurGcgcoefdiving(gcg) );
   SCIP_CALL( GCGincludeHeurGcgfracdiving(gcg) );
   SCIP_CALL( GCGincludeHeurGcgguideddiving(gcg) );
   SCIP_CALL( GCGincludeHeurGcglinesdiving(gcg) );
   SCIP_CALL( GCGincludeHeurGcgpscostdiving(gcg) );
   SCIP_CALL( GCGincludeHeurGcgveclendiving(gcg) );
   SCIP_CALL( GCGincludeHeurGcgdins(gcg) );
   SCIP_CALL( GCGincludeHeurGcgfeaspump(gcg) );
   SCIP_CALL( GCGincludeHeurGcgrens(gcg) );
   SCIP_CALL( GCGincludeHeurGcgrins(gcg) );
   SCIP_CALL( GCGincludeHeurGcgrounding(gcg) );
   SCIP_CALL( GCGincludeHeurGcgshifting(gcg) );
   SCIP_CALL( GCGincludeHeurGcgsimplerounding(gcg) );
   SCIP_CALL( GCGincludeHeurGcgzirounding(gcg) );
   SCIP_CALL( GCGincludeHeurXpcrossover(gcg) );
   SCIP_CALL( GCGincludeHeurXprins(gcg) );

   /* Jonas' stuff */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable conflict analysis since adding constraints after structure detection may destroy symmetries */
   SCIP_CALL( SCIPsetBoolParam(scip, "conflict/enable", FALSE) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/clique/freq", -1) );
   SCIP_CALL( SCIPfixParam(scip, "heuristics/clique/freq") );
   SCIP_CALL( SCIPfixParam(scip, "conflict/enable") );

   SCIP_CALL( GCGincludeDispGcg(gcg) );
   SCIP_CALL( GCGincludeDialogsGraph(gcg) );
   SCIP_CALL( SCIPincludeTableDefault(scip) );

   return SCIP_OKAY;
}
