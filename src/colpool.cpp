/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2014 Operations Research, RWTH Aachen University       */
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

/**@file   colpool.cpp
 * @brief  column pool
 * @author Jonas Witt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */
#include <cassert>
#include <cstring>

/*lint -e64 disable useless and wrong lint warning */

#ifdef __INTEL_COMPILER
#ifndef _OPENMP
#pragma warning disable 3180  /* disable wrong and useless omp warnings */
#endif
#endif
#include "scip/scip.h"
#include "gcg.h"

#include "scip/cons_linear.h"
#include "scip/cons_knapsack.h"

#include "pricer_gcg.h"
#include "objpricer_gcg.h"
#include "sepa_master.h"

#include "relax_gcg.h"
#include "struct_solver.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "cons_masterbranch.h"
#include "objscip/objscip.h"
#include "objpricer_gcg.h"
#include "class_pricingtype.h"
#include "class_stabilization.h"
#include "branch_generic.h"
#include "event_display.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace scip;




