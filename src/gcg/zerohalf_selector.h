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

/**@file   sepa_zerohalf.h
 * @ingroup SEPARATORS
 * @brief  {0,1/2}-cuts separator
 * @author Manuel Kutschka
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_ZEROHALF_SELECTOR_H__
#define GCG_ZEROHALF_SELECTOR_H__


#include <scip/def.h>
#include <scip/type_retcode.h>
#include <scip/def.h>
#include <scip/type_retcode.h>
#include <scip/type_scip.h>

#include "gcg/gcg.h"

#ifdef __cplusplus
extern "C" {
#endif

struct GCG_ZeroHalfData
{
   SCIP_Real             minviol;            /**< minimal violation to generate zerohalfcut for */
   SCIP_Real             maxslack;           /**< maximal slack of rows to be used in aggregation */
   SCIP_Real             maxslackroot;       /**< maximal slack of rows to be used in aggregation in the root node */
   SCIP_Real             maxrowdensity;      /**< maximal density of row to be used in aggregation */
   SCIP_Bool             infeasible;         /**< infeasibility was detected after adding a zerohalf cut */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   int                   maxrounds;          /**< maximal number of zerohalf separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of zerohalf separation rounds in the root node (-1: unlimited) */
   int                   densityoffset;      /**< additional number of variables allowed in row on top of density */
   int                   nreductions;        /**< number of reductions to the mod 2 system found so far */
   SCIP_CONS**           origmasterconss;    /**< array of orig master conss (is only allocated if not equal to data of relaxator) */
   SCIP_CONS**           masterconss;        /**< array of master conss (is only allocated if not equal to data of relaxator) */
   int                   nmasterconss;       /**< number of master conss */
};
typedef struct GCG_ZeroHalfData GCG_ZEROHALFDATA;

struct GCG_CutIndices
{
   int*  indices;       /**< indices of the the constraint used to create the cut */
   int   nindices;      /**< number of constraints used to create the cut */
};
typedef struct GCG_CutIndices GCG_CUTINDICES;

struct RowIndex
{
   unsigned int          type:2;             /**< type of row index; 0 means lp row using the right hand side,
                                              *   1 means lp row using the left hand side, and 2 means a
                                              *   transformed integral row */
   unsigned int          index:30;           /**< lp position of original row, or index of transformed integral row */
};
typedef struct RowIndex ROWINDEX;

/** perform the zerohalf cut separation */
SCIP_RETCODE GCGselectConstraintsZeroHalf(
   GCG*                  gcg,
   SCIP_SOL*             sol,
   SCIP_Bool             allowlocal,
   int                   depth,               /* current depth */
   GCG_ZEROHALFDATA*     zhdata,
   int                   ncalls,
   int                   maxcuts,
   GCG_CUTINDICES***     cutindices,
   int*                  ncutindices
);

/** create an instance of cut indices */
SCIP_RETCODE GCGcreateCutIndicesFromRowindex(
   SCIP*                scip,          /**< SCIP data structure (master problem) */
   GCG_CUTINDICES**     cutindices,    /**< pointer to store cutindices */
   int                  nindices,      /**< number of indices */
   ROWINDEX*            rowindex       /**< indices */
);

/** create an instance of cut indices */
SCIP_RETCODE GCGcreateCutIndicesFromArray(
   SCIP*                scip,          /**< SCIP data structure (master problem) */
   GCG_CUTINDICES**     cutindices,    /**< pointer to store cutindices */
   int                  nindices,      /**< number of indices */
   int*                 indices        /**< indices */
);

/** free an instance of cut indices */
SCIP_RETCODE GCGfreeCutIndices(
   SCIP*             scip,               /**< SCIP data structure (master problem) */
   GCG_CUTINDICES**  cutindices          /**< pointer to cutindices */
);

#ifdef __cplusplus
}
#endif

#endif
