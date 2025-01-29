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

/**@file   pub_clsvar.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for variable classifiers
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __GCG_PUB_CLSVAR_H__
#define __GCG_PUB_CLSVAR_H__


#include "def.h"
#include "type_varclassifier.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicVarClassifierMethods
 *
 * @{
 */

/** gets user data of variable classifier */
GCG_EXPORT
GCG_CLASSIFIERDATA* GCGvarClassifierGetData(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** sets user data of variable classifier; user has to free old data in advance! */
GCG_EXPORT
void GCGvarClassifierSetData(
   GCG_VARCLASSIFIER* varClassifier,                /**< variable classifier data structure */
   GCG_CLASSIFIERDATA* clsdata                      /**< new variable classifier user data */
   );

/** gets name of variable classifier */
GCG_EXPORT
const char* GCGvarClassifierGetName(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** gets priority of variable classifier */
GCG_EXPORT
int GCGvarClassifierGetPriority(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** returns True iff variable classifier is enabled */
GCG_EXPORT
SCIP_Bool GCGvarClassifierIsEnabled(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** gets description of variable classifier */
GCG_EXPORT
const char* GCGvarClassifierGetDesc(
   GCG_VARCLASSIFIER* varClassifier                 /**< variable classifier data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
