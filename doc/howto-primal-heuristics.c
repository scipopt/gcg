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

/**@file   howto_primal_heuristics.c
 * @brief  howto document page
 * @author   Gerald Gamrath
 * @author   Martin Bergner
 * @author   Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add primal heuristics
 *
 * For general information on how to add your own primal heuristics to GCG, first check the SCIP documentation.
 * However, one has to take into account some peculiarities when implementing heuristics that are included
 * in the original SCIP instance, i.e. work on the original variables.
 *
 * @section HEUR_LPSOLUTIONS Access to LP feasible solutions (on the original variables)
 *
 * Many MIP heuristics make use of an LP feasible solution. In SCIP, such a solution is obtained by solving the LP relaxation.
 * In GCG however, no LP relaxation is solved by default. A linearly feasible solution on the original variables comes from the
 * GCG relaxator plug-in; it is a solution of the master LP that has been translated back into the original variables. To access
 * it, one should use the method GCGrelaxGetCurrentOrigSol() in relax_gcg.h.
 * Its fractional variables can be accessed by SCIPgetExternBranchCands() (rather than SCIPgetLPBranchCands() which is used
 * in SCIP heuristics).
 * \n
 * Note also that heuristics using LP solutions should use another timing than SCIP heuristics. Heuristics that are called after
 * solving a node's relaxation typically have the timing SCIP_HEURTIMING_AFTERLPNODE.
 * By default, no LPs are solved on the original problem. A heuristic relying on a linearly feasible solution should therefore
 * have the timing SCIP_HEURTIMING_AFTERNODE to ensure that the heuristic is called at all. One then must ensure that the node's
 * relaxation has indeed been solved to optimality and that the relaxation solution is valid. This can be done by placing
 * \code
 * /* do not execute the heuristic on invalid relaxation solutions
 *  * (which is the case if the node has been cut off)
 *  */
 *  if( !SCIPisRelaxSolValid(scip) )
 *     return SCIP_OKAY;
 *
 *  /* only call heuristic, if an optimal LP solution is at hand */
 *  if( SCIPgetStage(GCGrelaxGetMasterprob(scip)) > SCIP_STAGE_SOLVING || SCIPgetLPSolstat(GCGrelaxGetMasterprob(scip)) != SCIP_LPSOLSTAT_OPTIMAL )
 *     return SCIP_OKAY;
 * \endcode
 * at the beginning of the HEUREXEC callback.
 *
 * @section HEUR_DIVING Diving on original variables
 *
 * A common class of heuristics are diving heuristics; they solve LPs with modified bounds to perform a depth-first
 * search on the Branch-and-bound tree. For this purpose, a probing mode and a diving mode have been implemented in SCIP,
 * which can be invoked by SCIPstartProbing(scip) and SCIPstartDive(scip), respectively. In this mode, temporary bound changes
 * on variables can be made, and modified LPs can be solved.
 * \n
 * In GCG, a special probing mode has been implemented for the original instance. This mode serves for performing changes on
 * the original instance but using the master LP instead of the original LP. It is invoked by GCGrelaxStartProbing() and terminated
 * by GCGrelaxEndProbing() and features the methods GCGrelaxPerformProbing() and GCGrelaxPerformProbingWithPricing(), which will
 * propagate any bound changes on the original instance to the extended instance and solve the resulting modified master LP, either
 * without or with pricing new variables in. See e.g. heur_gcgcoefdiving.c for an example on how to use them.
 *
 * @section HEURCOPY The HEURCOPY callback
 *
 * The HEURCOPY callback is executed when a SCIP instance is copied, e.g. to
 * solve a sub-SCIP. By
 * defining this callback as
 * <code>NULL</code> the user disables the execution of the specified
 * heuristic for all copied SCIP instances. This may deteriorate the performance
 * of primal heuristics using sub-SCIPs.
 * \n
 * For heuristics that are included in the original instance and make use of the extended instance as well (in
 * particular, most of the heur_gcg* and heur_xp* plug-ins), this callback should be set to null. This is because
 * sub-SCIPs are solved by SCIP rather than GCG and therefore do not know any master problem; including a GCG
 * specific heuristic into them would cause errors.
 *
 */
 
