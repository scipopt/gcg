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

/**@file   clscons.c
 * @ingroup DECOMP
 * @brief  interface for constraint classifier
 * @author Jurgen Lentz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"

#include "struct_consclassifier.h"
#include "pub_clscons.h"


/** gets user data of constraint classifier */
GCG_CLASSIFIERDATA* GCGconsClassifierGetData(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->clsdata;
}

/** sets user data of constraint classifier; user has to free old data in advance! */
void GCGconsClassifierSetData(
   GCG_CONSCLASSIFIER* consClassifier,              /**< constraint classifier data structure */
   GCG_CLASSIFIERDATA* clsdata                      /**< new constraint classifier user data */
   )
{
   assert(consClassifier != NULL);

   consClassifier->clsdata = clsdata;
}

/** gets name of constraint classifier */
const char* GCGconsClassifierGetName(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->name;
}

/** gets priority of constraint classifier */
int GCGconsClassifierGetPriority(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->priority;
}

/** returns True iff constraint classifier is enabled */
SCIP_Bool GCGconsClassifierIsEnabled(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->enabled;
}

/** gets description of constraint classifier */
const char* GCGconsClassifierGetDesc(
   GCG_CONSCLASSIFIER* consClassifier  /**< constraint classifier data structure */
   )
{
   assert(consClassifier != NULL);

   return consClassifier->description;
}
