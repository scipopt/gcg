/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file   struct_mastersepacutdata.h
 * @ingroup DATASTRUCTURES
 * @brief  data structures for GCG separator cuts
 * @author Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_STRUCT_MASTERSEPACUTDATA_H__
#define GCG_STRUCT_MASTERSEPACUTDATA_H__

#include <scip/def.h>

#include "gcg/type_extendedmasterconsdata.h"
#include "gcg/type_mastersepacut.h"
#include "gcg/type_gcgvarhistory.h"
#include "gcg/type_sepagcg.h"

#ifdef __cplusplus
extern "C" {
#endif


/** additional data for subset row cuts */
struct GCG_ChvatalGomoryCutData
{
   int*                    conssindices;           /**< indices of constraints used to create cut */
   SCIP_Real*              weights;                /**< weights used to create cut */
   int                     nconssindices;          /**< number of constraints used to create cut */

};

/** additional data for master separator cuts */
struct GCG_SeparatorMasterCutData
{
   union
   {
      GCG_CHVATALGOMORYCUTDATA    chvatalgomorycutdata;       /**< data for Chvatal-Gomory cut */
   } data;
};

/** master separator cut data structure */
struct GCG_SeparatorMasterCut
{
   GCG_SEPA*                        sepa;                   /**< GCG Master Separator which generated the mastercut */
   GCG_SEPARATORMASTERCUTTYPE       type;                   /**< type of mastercut */
   GCG_VARHISTORY*                  knownvarhistory;        /**< pointer to the history of priced variables */
   GCG_SEPARATORMASTERCUTDATA*      data;                   /**< additional data helpful to compute coefficients */
   int                              nuses;                  /**< number of times this separator mastercut is referenced */
};

#ifdef __cplusplus
}
#endif

#endif