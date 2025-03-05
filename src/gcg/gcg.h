/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   gcg.h
 * @ingroup PUBLICCOREAPI
 * @brief  GCG interface methods
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_STATISTIC */

#ifndef GCG_H_
#define GCG_H_

#include "scip/scip.h"
#include "gcg/def.h"

#include "gcg/type_branchgcg.h"
#include "gcg/type_classifier.h"
#include "gcg/type_colpool.h"
#include "gcg/type_consclassifier.h"
#include "gcg/type_decomp.h"
#include "gcg/type_detector.h"
#include "gcg/type_extendedmasterconsdata.h"
#include "gcg/type_gcg.h"
#include "gcg/type_gcgcol.h"
#include "gcg/type_gcgpqueue.h"
#include "gcg/type_masterdiving.h"
#include "gcg/type_origdiving.h"
#include "gcg/type_parameter.h"
#include "gcg/type_pricingjob.h"
#include "gcg/type_pricingprob.h"
#include "gcg/type_pricingstatus.h"
#include "gcg/type_score.h"
#include "gcg/type_solver.h"
#include "gcg/type_varclassifier.h"

/* include public interfaces, s.t. the user only needs to include gcg.h */
#include "gcg/pub_clscons.h"
#include "gcg/pub_clsvar.h"
#include "gcg/pub_colpool.h"
#include "gcg/pub_decomp.h"
#include "gcg/pub_extendedmasterconsdata.h"
#include "gcg/pub_gcg.h"
#include "gcg/pub_gcgcol.h"
#include "gcg/pub_gcgheur.h"
#include "gcg/pub_gcgpqueue.h"
#include "gcg/pub_gcgsepa.h"
#include "gcg/pub_gcgvar.h"
#include "gcg/pub_pricingjob.h"
#include "gcg/pub_pricingprob.h"
#include "gcg/pub_score.h"
#include "gcg/pub_solver.h"

#include "gcg/relax_gcg.h"
#include "gcg/gcg_general.h"

#endif /* GCG_H_ */
