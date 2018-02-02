/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   gcgplugins.c
 * @brief  SCIP plugins for generic column generation
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcgplugins.h"
#include "scip/debug.h"

#define USEHEURS 1
#define USESEPA 1
#define USEPROP 1

/* include header files here, such that the user only has to include
 * gcgplugins.h
 */
#include "scip/cons_integral.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"

#if USEHEURS
#include "scip/heur_actconsdiving.h"
#include "scip/heur_clique.h"
#include "scip/heur_coefdiving.h"
#include "scip/heur_crossover.h"
#include "scip/heur_dins.h"
#include "scip/heur_dualval.h"
#include "scip/heur_feaspump.h"
#include "scip/heur_fixandinfer.h"
#include "scip/heur_fracdiving.h"
#include "scip/heur_guideddiving.h"
#include "scip/heur_intdiving.h"
#include "scip/heur_intshifting.h"
#include "scip/heur_linesearchdiving.h"
#include "scip/heur_localbranching.h"
#include "scip/heur_mutation.h"
#include "scip/heur_nlpdiving.h"
#include "scip/heur_objpscostdiving.h"
#include "scip/heur_octane.h"
#include "scip/heur_oneopt.h"
#include "scip/heur_proximity.h"
#include "scip/heur_pscostdiving.h"
#include "scip/heur_randrounding.h"
#include "scip/heur_rens.h"
#include "scip/heur_rins.h"
#include "scip/heur_rootsoldiving.h"
#include "scip/heur_rounding.h"
#include "scip/heur_shiftandpropagate.h"
#include "scip/heur_shifting.h"
#include "scip/heur_simplerounding.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trivial.h"
#include "scip/heur_trysol.h"
#include "scip/heur_twoopt.h"
#include "scip/heur_undercover.h"
#include "scip/heur_vbounds.h"
#include "scip/heur_veclendiving.h"
#include "scip/heur_zirounding.h"
#include "scip/heur_zeroobj.h"
#endif

#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"

#include "scip/presol_implics.h"
#include "scip/presol_inttobinary.h"
#include "scip/presol_trivial.h"
#include "scip/presol_boundshift.h"
#include "scip/presol_domcol.h"
#include "scip/presol_convertinttobin.h"

#if USEPROP
#include "scip/prop_dualfix.h"
#include "scip/prop_probing.h"
#include "scip/prop_pseudoobj.h"
#include "scip/prop_rootredcost.h"
#include "scip/prop_redcost.h"
#include "scip/prop_genvbounds.h"
#include "scip/prop_vbounds.h"
#include "scip/prop_obbt.h"
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
#include "scip/sepa_clique.h"
#include "scip/sepa_gomory.h"
#include "scip/sepa_impliedbounds.h"
#include "scip/sepa_intobj.h"
#include "scip/sepa_mcf.h"
#include "scip/sepa_oddcycle.h"
#include "scip/sepa_strongcg.h"
#include "scip/sepa_zerohalf.h"

/** added by Jonas */
#include "scip/sepa_closecuts.h"
#include "scip/sepa_rapidlearning.h"
#endif



#include "scip/scipshell.h"
#include "reader_blk.h"
#include "reader_dec.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "branch_empty.h"

#include "cons_origbranch.h"
#include "disp_gcg.h"
#include "dialog_gcg.h"
#include "reader_ref.h"
#include "event_bestsol.h"
#include "event_mastersol.h"

/* Martin's detection stuff */
#include "reader_gp.h"
#include "cons_decomp.h"
#include "dec_connected.h"

#ifndef NBLISS
#include "dec_isomorph.h"
#endif

#include "dec_arrowheur.h"
#include "dec_stairheur.h"
#include "dec_staircase.h"
#include "dec_random.h"
#include "dec_colors.h"

/* Christian's heuristics */
#include "heur_gcgcoefdiving.h"
#include "heur_gcgdins.h"
#include "heur_gcgfeaspump.h"
#include "heur_gcgfracdiving.h"
#include "heur_gcgguideddiving.h"
#include "heur_gcglinesdiving.h"
#include "heur_gcgpscostdiving.h"
#include "heur_gcgrens.h"
#include "heur_gcgrins.h"
#include "heur_gcgrounding.h"
#include "heur_gcgshifting.h"
#include "heur_gcgsimplerounding.h"
#include "heur_gcgveclendiving.h"
#include "heur_gcgzirounding.h"
#include "heur_origdiving.h"
#include "heur_xpcrossover.h"
#include "heur_xprins.h"

/* Friedrike's detection stuff */
#include "dec_cutpacking.h"
#include "scip_misc.h"
#include "scip/table_default.h"


/** includes default plugins for generic column generation into SCIP */
SCIP_RETCODE SCIPincludeGcgPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) ); /* linear must be first due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrIntegral(scip) );
   SCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );
   SCIP_CALL( SCIPincludeConshdlrLogicor(scip) );
   SCIP_CALL( SCIPincludeConshdlrSetppc(scip) );
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );

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
   SCIP_CALL( SCIPincludePresolImplics(scip) );
   SCIP_CALL( SCIPincludePresolInttobinary(scip) );
   SCIP_CALL( SCIPincludePresolTrivial(scip) );
   SCIP_CALL( SCIPincludePresolDomcol(scip) );
   SCIP_CALL( SCIPincludePresolConvertinttobin(scip) );

   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   SCIP_CALL( SCIPincludeNodeselEstimate(scip) );
   SCIP_CALL( SCIPincludeNodeselHybridestim(scip) );
   SCIP_CALL( SCIPincludeNodeselRestartdfs(scip) );

#if USEPROP
   SCIP_CALL( SCIPincludePropDualfix(scip) );
   SCIP_CALL( SCIPincludePropPseudoobj(scip) );
   SCIP_CALL( SCIPincludePropRootredcost(scip) );
   SCIP_CALL( SCIPincludePropGenvbounds(scip) );
   SCIP_CALL( SCIPincludePropProbing(scip) );
   SCIP_CALL( SCIPincludePropRedcost(scip) );
   SCIP_CALL( SCIPincludePropVbounds(scip) );
   SCIP_CALL( SCIPincludePropObbt(scip) );
#endif


#if USEHEURS
   SCIP_CALL( SCIPincludeHeurActconsdiving(scip) );
   SCIP_CALL( SCIPincludeHeurClique(scip) );
   SCIP_CALL( SCIPincludeHeurCoefdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCrossover(scip) );
   SCIP_CALL( SCIPincludeHeurDins(scip) );
   SCIP_CALL( SCIPincludeHeurDualval(scip) );
   SCIP_CALL( SCIPincludeHeurFeaspump(scip) );
   SCIP_CALL( SCIPincludeHeurFixandinfer(scip) );
   SCIP_CALL( SCIPincludeHeurFracdiving(scip) );
   SCIP_CALL( SCIPincludeHeurGuideddiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntdiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntshifting(scip) );
   SCIP_CALL( SCIPincludeHeurLinesearchdiving(scip) );
   SCIP_CALL( SCIPincludeHeurLocalbranching(scip) );
   SCIP_CALL( SCIPincludeHeurMutation(scip) );
   SCIP_CALL( SCIPincludeHeurNlpdiving(scip) );
   SCIP_CALL( SCIPincludeHeurObjpscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurOctane(scip) );
   SCIP_CALL( SCIPincludeHeurOneopt(scip) );
   SCIP_CALL( SCIPincludeHeurProximity(scip) );
   SCIP_CALL( SCIPincludeHeurPscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurRandrounding(scip) );
   SCIP_CALL( SCIPincludeHeurRens(scip) );
   SCIP_CALL( SCIPincludeHeurRins(scip) );
   SCIP_CALL( SCIPincludeHeurRootsoldiving(scip) );
   SCIP_CALL( SCIPincludeHeurRounding(scip) );
   SCIP_CALL( SCIPincludeHeurShiftandpropagate(scip) );
   SCIP_CALL( SCIPincludeHeurShifting(scip) );
   SCIP_CALL( SCIPincludeHeurSubNlp(scip) );
   SCIP_CALL( SCIPincludeHeurTrivial(scip) );
   SCIP_CALL( SCIPincludeHeurTrySol(scip) );
   SCIP_CALL( SCIPincludeHeurTwoopt(scip) );
   SCIP_CALL( SCIPincludeHeurUndercover(scip) );
   SCIP_CALL( SCIPincludeHeurVbounds(scip) );
   SCIP_CALL( SCIPincludeHeurVeclendiving(scip) );
   SCIP_CALL( SCIPincludeHeurZirounding(scip) );
   SCIP_CALL( SCIPincludeHeurZeroobj(scip) );
#endif
   SCIP_CALL( SCIPincludeHeurSimplerounding(scip) );

#if USESEPA
   SCIP_CALL( SCIPincludeSepaClique(scip) );
   SCIP_CALL( SCIPincludeSepaGomory(scip) );
   SCIP_CALL( SCIPincludeSepaImpliedbounds(scip) );
   SCIP_CALL( SCIPincludeSepaIntobj(scip) );
   SCIP_CALL( SCIPincludeSepaMcf(scip) );
   SCIP_CALL( SCIPincludeSepaOddcycle(scip) );
   SCIP_CALL( SCIPincludeSepaStrongcg(scip) );
   SCIP_CALL( SCIPincludeSepaZerohalf(scip) );

   /* added by Jonas */
   SCIP_CALL( SCIPincludeSepaClosecuts(scip) );
   SCIP_CALL( SCIPincludeSepaRapidlearning(scip) );
#endif

   SCIP_CALL( SCIPincludeRelaxGcg(scip) );
   SCIP_CALL( SCIPincludeReaderBlk(scip) );
   SCIP_CALL( SCIPincludeReaderDec(scip) );
   SCIP_CALL( SCIPincludeReaderRef(scip) );
   SCIP_CALL( SCIPincludeBranchruleEmpty(scip) );

   SCIP_CALL( SCIPincludeConshdlrOrigbranch(scip) );
   SCIP_CALL( SCIPincludeEventHdlrBestsol(scip) );
   SCIP_CALL( SCIPincludeEventHdlrMastersol(scip) );

   /* Detectors and decompositions */
   SCIP_CALL( SCIPincludeReaderGp(scip) );
   SCIP_CALL( SCIPincludeConshdlrDecomp(scip) );
   SCIP_CALL( SCIPincludeDetectorConnected(scip) );
   SCIP_CALL( SCIPincludeDetectorArrowheur(scip) );
   SCIP_CALL( SCIPincludeDetectorStairheur(scip) );
   SCIP_CALL( SCIPincludeDetectorStaircase(scip) );
   SCIP_CALL( SCIPincludeDetectorRandom(scip) );
   SCIP_CALL( SCIPincludeDetectorColors(scip) );
   SCIP_CALL( SCIPincludeDetectorCutpacking(scip) );


#ifndef NBLISS
   SCIP_CALL( SCIPincludeDetectorIsomorphism(scip) );
#endif

   /* Christian's heuristics */
   SCIP_CALL( SCIPincludeEventHdlrOrigdiving(scip) );
   SCIP_CALL( GCGincludeHeurGcgcoefdiving(scip) );
   SCIP_CALL( GCGincludeHeurGcgfracdiving(scip) );
   SCIP_CALL( GCGincludeHeurGcgguideddiving(scip) );
   SCIP_CALL( GCGincludeHeurGcglinesdiving(scip) );
   SCIP_CALL( GCGincludeHeurGcgpscostdiving(scip) );
   SCIP_CALL( GCGincludeHeurGcgveclendiving(scip) );
   SCIP_CALL( SCIPincludeHeurGcgdins(scip) );
   SCIP_CALL( SCIPincludeHeurGcgfeaspump(scip) );
   SCIP_CALL( SCIPincludeHeurGcgrens(scip) );
   SCIP_CALL( SCIPincludeHeurGcgrins(scip) );
   SCIP_CALL( SCIPincludeHeurGcgrounding(scip) );
   SCIP_CALL( SCIPincludeHeurGcgshifting(scip) );
   SCIP_CALL( SCIPincludeHeurGcgsimplerounding(scip) );
   SCIP_CALL( SCIPincludeHeurGcgzirounding(scip) );
   SCIP_CALL( SCIPincludeHeurXpcrossover(scip) );
   SCIP_CALL( SCIPincludeHeurXprins(scip) );

   /* Jonas' stuff */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIP_CALL( SCIPincludeDispGcg(scip) );
   SCIP_CALL( SCIPincludeDialogGcg(scip) );
   SCIP_CALL( GCGincludeDialogsGraph(scip) );
   SCIP_CALL( SCIPincludeTableDefault(scip) );

   return SCIP_OKAY;
}
