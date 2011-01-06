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

/**@file   gcgplugins.h
 * @brief  SCIP plugins for coloring
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIPGCGPLUGINS_H__
#define __SCIP_SCIPGCGPLUGINS_H__


#include "scip/scip.h"

/* include header files here, such that the user only has to include
 * gcgplugins.h
 */
#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"
#include "scip/cons_and.h"
#include "scip/cons_bounddisjunction.h"
#include "scip/cons_conjunction.h"
#include "scip/cons_integral.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_or.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"
//#include "scip/dialog_default.h"
#include "scip/heur_actconsdiving.h"
#include "scip/heur_coefdiving.h"
#include "scip/heur_crossover.h"
#include "scip/heur_feaspump.h"
#include "scip/heur_fixandinfer.h"
#include "scip/heur_fracdiving.h"
#include "scip/heur_guideddiving.h"
#include "scip/heur_intdiving.h"
#include "scip/heur_intshifting.h"
#include "scip/heur_linesearchdiving.h"
#include "scip/heur_localbranching.h"
#include "scip/heur_mutation.h"
#include "scip/heur_objpscostdiving.h"
#include "scip/heur_octane.h"
#include "scip/heur_oneopt.h"
#include "scip/heur_pscostdiving.h"
#include "scip/heur_rens.h"
#include "scip/heur_rins.h"
#include "scip/heur_rootsoldiving.h"
#include "scip/heur_rounding.h"
#include "scip/heur_shifting.h"
#include "scip/heur_simplerounding.h"
#include "scip/heur_veclendiving.h"

#include "heur_gcgfracdiving.h"

#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"
#include "scip/presol_dualfix.h"
#include "scip/presol_implics.h"
#include "scip/presol_inttobinary.h"
#include "scip/presol_probing.h"
#include "scip/presol_trivial.h"
#include "scip/presol_boundshift.h"
#include "scip/prop_pseudoobj.h"
#include "scip/prop_rootredcost.h"
#include "scip/reader_cnf.h"
#include "scip/reader_fix.h"
#include "scip/reader_lp.h"
#include "scip/reader_mps.h"
#include "scip/reader_sol.h"
#include "scip/reader_zpl.h"
#include "scip/reader_cip.h"
#include "scip/sepa_clique.h"
#include "scip/sepa_cmir.h"
#include "scip/sepa_flowcover.h"
#include "scip/sepa_gomory.h"
#include "scip/sepa_impliedbounds.h"
#include "scip/sepa_intobj.h"
#include "scip/sepa_mcf.h"
#include "scip/sepa_redcost.h"
#include "scip/sepa_strongcg.h"
#include "scip/sepa_zerohalf.h"
#include "scip/scipshell.h"
#include "scip/disp_default.h"
#include "reader_blk.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "branch_orig.h"
#include "branch_ryanfoster.h"
#include "cons_origbranch.h"
#include "disp_gcg.h"
#include "dialog_gcg.h"
#include "branch_relpsprob.h"



/** includes default SCIP plugins into SCIP */
extern
SCIP_RETCODE SCIPincludeGcgPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   );


#endif
