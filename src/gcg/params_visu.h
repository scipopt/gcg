/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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

/**@file   reader_tex.h
 * @brief  parameter settings for visualization readers
 * @author Hanna Franzen

 * This file provides universally used parameters for visualizations.

 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_PARAMS_VISU_H__
#define GCG_PARAMS_VISU_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "scip/scip.h"
#include "gcg/gcg.h"

/** Colorscheme selection for the visualizations */
enum Colorscheme
{
   COLORSCHEME_DEFAULT  = 0,     /**< default colors (supposedly eye-friendly) */
   COLORSCHEME_GREY     = 1,     /**< on a range from black to white */
   COLORSCHEME_MANUAL   = 2      /**< take user-defined input */
};

typedef enum Colorscheme VISU_COLORSCHEME; /**< visualization colorscheme type */

/** includes the visualization parameters into GCG
 * @returns SCIP return code */
GCG_EXPORT
SCIP_RETCODE GCGcreateParamsVisu(
   GCG*  gcg,                 /**< GCG data structure */
   GCG_PARAMDATA** paramdata  /**< input empty paramdata, oputput new set of param data */
   );

/** gets whether draftmode is on
 * draftmode lets visualizations omit nonzeros
 * @returns true if draftmode is on */
GCG_EXPORT
SCIP_Bool GCGvisuGetDraftmode(
   GCG*  gcg   /**< GCG data structure */
   );

/** sets draftmode
 * draftmode lets visualizations omit nonzeros
 */
GCG_EXPORT
void GCGvisuSetDraftmode(
   GCG*  gcg,        /**< GCG data structure */
   SCIP_Bool setmode /**< true iff draftmode should be on */
   );

/** gets the colorscheme for visualizations
 * @returns current colorscheme */
GCG_EXPORT
VISU_COLORSCHEME GCGvisuGetColorscheme(
   GCG*  gcg   /**< GCG data structure */
   );

/** sets colorscheme for visualizations
 * 
 */
GCG_EXPORT
void GCGvisuSetColorscheme(
   GCG*  gcg,                    /**< GCG data structure */
   VISU_COLORSCHEME newscheme    /**< new colorscheme */
   );

/** sets color for mastercons block in current color scheme
 *
 */
GCG_EXPORT
void GCGvisuSetColorManMasterconss(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for mastervar block in current color scheme
 * 
 */
GCG_EXPORT
void GCGvisuSetColorManMastervars(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for linking blocks in current color scheme
 * 
 */
GCG_EXPORT
void GCGvisuSetColorManLinking(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for stairlinking blocks in current color scheme
 * 
 */
GCG_EXPORT
void GCGvisuSetColorManStairlinking(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for normal decomp blocks in current color scheme
 * 
 */
GCG_EXPORT
void GCGvisuSetColorManBlock(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for open blocks in current color scheme
 * 
 */
GCG_EXPORT
void GCGvisuSetColorManOpen(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for non-zero points in current color scheme
 * 
 */
GCG_EXPORT
void GCGvisuSetColorManNonzero(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for lines in current color scheme
 * 
 */
GCG_EXPORT
void GCGvisuSetColorManLine(
   GCG*  gcg,           /**< GCG data structure */
   const char* newcolor /**< new color */
   );

/** gets color for mastercons block in current color scheme
 * @returns mastercons color */
GCG_EXPORT
const char* GCGvisuGetColorMasterconss(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets color for mastervar block in current color scheme
 * @returns mastervar color */
GCG_EXPORT
const char* GCGvisuGetColorMastervars(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets color for linking blocks in current color scheme
 * @returns linking color */
GCG_EXPORT
const char* GCGvisuGetColorLinking(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets color for stairlinking blocks in current color scheme
 * @returns stairlinking color */
GCG_EXPORT
const char* GCGvisuGetColorStairlinking(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets color for normal decomp blocks in current color scheme
 * @returns block color */
GCG_EXPORT
const char* GCGvisuGetColorBlock(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets color for open blocks in current color scheme
 * @returns open color */
GCG_EXPORT
const char* GCGvisuGetColorOpen(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets color for non-zero points in current color scheme
 * @returns non-zero color */
GCG_EXPORT
const char* GCGvisuGetColorNonzero(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets color for lines in current color scheme
 * @returns line color */
GCG_EXPORT
const char* GCGvisuGetColorLine(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets appropriate radius for nonzeros
 * needs highest indices of both axes for scaling
 * @returns radius */
GCG_EXPORT
float GCGvisuGetNonzeroRadius(
   GCG*  gcg,           /**< GCG data structure */
   int maxindx,         /**< highest index x-axis */
   int maxindy,         /**< highest index y-axis */
   float scalingfactor  /**< percentage to scale radius, 1 if no scaling */
   );


/** if true gp reader should be used for sub-visualizations, otherwise tex reader
 * @returns true if gp reader should be used, false if tex reader should be used */
GCG_EXPORT
SCIP_Bool GCGgetUseGp(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets the name of the pdf reader that should be used
 * @returns name of pdf reader */
GCG_EXPORT
const char* GCGVisuGetPdfReader(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets the max number of decomps to be included in reports
 * @returns max number of decomps */
GCG_EXPORT
int GCGreportGetMaxNDecomps(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets what type of decomps to show in reports (where 0 corresponds to 'show all')
 * @returns type of decomps */
GCG_EXPORT
GCG_DECTYPE GCGreportGetDecompTypeToShow(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets whether a titlepage should be included in reports
 * @returns true iff title page should be generated */
GCG_EXPORT
SCIP_Bool GCGreportGetShowTitlepage(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets whether a table of contents should be included in reports
 * @returns true iff table of contents should be generated */
GCG_EXPORT
SCIP_Bool GCGreportGetShowToc(
   GCG*  gcg   /**< GCG data structure */
   );

/** gets whether statistics should be included for each decomp in reports
 * @returns true iff statistics for each decomp should be generated */
GCG_EXPORT
SCIP_Bool GCGreportGetShowStatistics(
   GCG*  gcg   /**< GCG data structure */
   );

/** frees all visualization parameters
 * 
 */
GCG_EXPORT
void GCGVisuFreeParams(
   GCG*  gcg,                 /**< GCG data structure */
   GCG_PARAMDATA* paramdata   /**< input empty paramdata, oputput new set of param data */
   );

#ifdef __cplusplus
}
#endif

#endif
