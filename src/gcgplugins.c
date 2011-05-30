/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"

/**@file   gcgplugins.c
 * @brief  SCIP plugins for generic column generation
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcgplugins.h"
#include "scip/debug.h"

#define USEHEURS 1
#define USESEPA 1


/** includes default plugins for generic column generation into SCIP */
SCIP_RETCODE SCIPincludeGcgPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) ); /* linear must be first due to constraint upgrading */
   SCIP_CALL( SCIPincludeConshdlrAnd(scip) );
   SCIP_CALL( SCIPincludeConshdlrBounddisjunction(scip) );
   SCIP_CALL( SCIPincludeConshdlrConjunction(scip) );
   SCIP_CALL( SCIPincludeConshdlrIntegral(scip) );
   SCIP_CALL( SCIPincludeConshdlrKnapsack(scip) );
   SCIP_CALL( SCIPincludeConshdlrLogicor(scip) );
   SCIP_CALL( SCIPincludeConshdlrOr(scip) );
   SCIP_CALL( SCIPincludeConshdlrSetppc(scip) );
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );
   SCIP_CALL( SCIPincludeConshdlrXor(scip) );

   SCIP_CALL( SCIPincludeReaderCip(scip) );
   SCIP_CALL( SCIPincludeReaderLp(scip) );
   SCIP_CALL( SCIPincludeReaderSol(scip) );
   SCIP_CALL( SCIPincludeReaderMps(scip) );

   SCIP_CALL( SCIPincludePresolBoundshift(scip) );
   SCIP_CALL( SCIPincludePresolDualfix(scip) );
   SCIP_CALL( SCIPincludePresolImplics(scip) );
   SCIP_CALL( SCIPincludePresolInttobinary(scip) );
   SCIP_CALL( SCIPincludePresolProbing(scip) );
   SCIP_CALL( SCIPincludePresolTrivial(scip) );

   //SCIP_CALL( SCIPincludeHeurGcgfeaspump(scip) );
   //SCIP_CALL( SCIPincludeHeurGcgrounding(scip) );

   SCIP_CALL( SCIPincludeNodeselBfs(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );

#if USEHEURS
   SCIP_CALL( SCIPincludeHeurActconsdiving(scip) );
   SCIP_CALL( SCIPincludeHeurClique(scip) );
   SCIP_CALL( SCIPincludeHeurCoefdiving(scip) );
   SCIP_CALL( SCIPincludeHeurCrossover(scip) );
   SCIP_CALL( SCIPincludeHeurDins(scip) );
   SCIP_CALL( SCIPincludeHeurFeaspump(scip) );
   SCIP_CALL( SCIPincludeHeurFixandinfer(scip) );
   SCIP_CALL( SCIPincludeHeurFracdiving(scip) );
   SCIP_CALL( SCIPincludeHeurGuideddiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntdiving(scip) );
   SCIP_CALL( SCIPincludeHeurIntshifting(scip) );
   SCIP_CALL( SCIPincludeHeurLinesearchdiving(scip) );
   SCIP_CALL( SCIPincludeHeurLocalbranching(scip) );
   SCIP_CALL( SCIPincludeHeurMutation(scip) );
   SCIP_CALL( SCIPincludeHeurObjpscostdiving(scip) );
   SCIP_CALL( SCIPincludeHeurOctane(scip) );
   SCIP_CALL( SCIPincludeHeurOneopt(scip) );
   SCIP_CALL( SCIPincludeHeurPscostdiving(scip) );
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
#endif
   SCIP_CALL( SCIPincludeHeurSimplerounding(scip) );


   SCIP_CALL( SCIPincludePropPseudoobj(scip) );
   SCIP_CALL( SCIPincludePropRootredcost(scip) );

#if USESEPA
   SCIP_CALL( SCIPincludeSepaClique(scip) );
   SCIP_CALL( SCIPincludeSepaCmir(scip) );
   SCIP_CALL( SCIPincludeSepaFlowcover(scip) );
   SCIP_CALL( SCIPincludeSepaGomory(scip) );
   SCIP_CALL( SCIPincludeSepaImpliedbounds(scip) );
   SCIP_CALL( SCIPincludeSepaIntobj(scip) );
   SCIP_CALL( SCIPincludeSepaMcf(scip) );
   SCIP_CALL( SCIPincludeSepaRedcost(scip) );
   SCIP_CALL( SCIPincludeSepaStrongcg(scip) );
   SCIP_CALL( SCIPincludeSepaZerohalf(scip) );
#endif

   SCIP_CALL( SCIPincludeRelaxGcg(scip) );
   SCIP_CALL( SCIPincludeReaderBlk(scip) );
   SCIP_CALL( SCIPincludeReaderRef(scip) );
   SCIP_CALL( SCIPincludeBranchruleOrig(scip) );
   SCIP_CALL( SCIPincludeBranchruleRyanfoster(scip) );
   SCIP_CALL( SCIPincludeBranchruleRelpsprob(scip) );
   SCIP_CALL( SCIPincludeConshdlrOrigbranch(scip) );
   SCIP_CALL( SCIPincludeDispGcg(scip) );

   //SCIP_CALL( SCIPincludeDispDefault(scip) );
   //SCIP_CALL( SCIPincludeDialogDefault(scip) );
   SCIP_CALL( SCIPincludeDialogGcg(scip) );
   //SCIP_CALL( SCIPdebugIncludeProp(scip) ); 

   /* martin's decomp stuff */
   SCIP_CALL(SCIPincludeReaderGp(scip));
   SCIP_CALL(SCIPincludeConshdlrDecomp(scip));

   /* Christian's heuristics */
   SCIP_CALL( SCIPincludeHeurGcgcoefdiving(scip) );
   SCIP_CALL( SCIPincludeHeurGcgfracdiving(scip) );
   SCIP_CALL( SCIPincludeHeurGcgrens(scip) );
   //SCIP_CALL( SCIPincludeHeurGcgrounding(scip) );
   SCIP_CALL( SCIPincludeHeurGcgshifting(scip) );
   SCIP_CALL( SCIPincludeHeurGcgsimplerounding(scip) );
   SCIP_CALL( SCIPincludeHeurGcgzirounding(scip) );



   return SCIP_OKAY;
}
