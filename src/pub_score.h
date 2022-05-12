/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2022 Operations Research, RWTH Aachen University       */
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

/**@file   pub_score.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for score
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_PUB_SCORE_H__
#define __GCG_PUB_SCORE_H__


#include "scip/def.h"
#include "type_score.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicScoreMethods
 *
 * @{
 */

/** gets user data of score */
SCIP_EXPORT
DEC_SCOREDATA* GCGscoreGetData(
   DEC_SCORE*            score               /**< score */
   );

/** sets user data of score; user has to free old data in advance! */
SCIP_EXPORT
void GCGscoreSetData(
   DEC_SCORE*            score,              /**< score */
   DEC_SCOREDATA*        scoredata           /**< new score user data */
   );

/** gets name of score */
SCIP_EXPORT
const char* GCGscoreGetName(
   DEC_SCORE*            score               /**< score */
   );

/** gets shortname of score */
SCIP_EXPORT
const char* GCGscoreGetShortname(
   DEC_SCORE*            score               /**< score */
   );

/** gets description of score */
SCIP_EXPORT
const char* GCGscoreGetDesc(
   DEC_SCORE*            score               /**< score */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
