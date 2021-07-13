/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2021 Operations Research, RWTH Aachen University       */
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

/**@file   type_varclassifier.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for variable classifier in GCG projects
 * @author William Ma
 */

#ifndef GCG_TYPE_VARCLASSIFIER_H__
#define GCG_TYPE_VARCLASSIFIER_H__

#include <scip/def.h>
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"

typedef struct DEC_VarClassifier DEC_VARCLASSIFIER;

/** destructor of classifier to free classifier data (called when GCG is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - classifier      : classifier data structure
 */
#define DEC_DECL_FREEVARCLASSIFIER(x) SCIP_RETCODE x (SCIP* scip, DEC_VARCLASSIFIER* classifier)

/**
 * Tries to classify variables with data of the according detprobdata and store the classification in the detprobdata
 *
 * input:
 *  - scip                 : SCIP data structure
 *  - transformed          : should use data from transformed detprobdata or not
 */
#define DEC_DECL_VARCLASSIFY(x) SCIP_RETCODE x (SCIP* scip, SCIP_Bool transformed)

#endif //GCG_TYPE_VARCLASSIFIER_H__