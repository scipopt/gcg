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

#include "def.h"
#include "type_decomp.h"
#include "type_parameter.h"

#include "scip/scip.h"

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
SCIP_RETCODE SCIPcreateParamsVisu(
   SCIP* scip,                /**< SCIP data structure */
   GCG_PARAMDATA** paramdata  /**< input empty paramdata, oputput new set of param data */
   );

/** gets whether draftmode is on
 * draftmode lets visualizations omit nonzeros
 * @returns true if draftmode is on */
GCG_EXPORT
SCIP_Bool SCIPvisuGetDraftmode(
   SCIP* scip  /**< SCIP data structure */
   );

/** sets draftmode
 * draftmode lets visualizations omit nonzeros
 */
GCG_EXPORT
void SCIPvisuSetDraftmode(
   SCIP* scip,       /**< SCIP data structure */
   SCIP_Bool setmode /**< true iff draftmode should be on */
   );

/** gets the colorscheme for visualizations
 * @returns current colorscheme */
GCG_EXPORT
VISU_COLORSCHEME SCIPvisuGetColorscheme(
   SCIP* scip  /**< SCIP data structure */
   );

/** sets colorscheme for visualizations
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorscheme(
   SCIP* scip,                   /**< SCIP data structure */
   VISU_COLORSCHEME newscheme    /**< new colorscheme */
   );

/** sets color for mastercons block in current color scheme
 *
 */
GCG_EXPORT
void SCIPvisuSetColorManMasterconss(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for mastervar block in current color scheme
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorManMastervars(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for linking blocks in current color scheme
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorManLinking(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for stairlinking blocks in current color scheme
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorManStairlinking(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for normal decomp blocks in current color scheme
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorManBlock(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for open blocks in current color scheme
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorManOpen(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for non-zero points in current color scheme
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorManNonzero(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** sets manual color for lines in current color scheme
 * 
 */
GCG_EXPORT
void SCIPvisuSetColorManLine(
   SCIP* scip,          /**< SCIP data structure */
   const char* newcolor /**< new color */
   );

/** gets color for mastercons block in current color scheme
 * @returns mastercons color */
GCG_EXPORT
const char* SCIPvisuGetColorMasterconss(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets color for mastervar block in current color scheme
 * @returns mastervar color */
GCG_EXPORT
const char* SCIPvisuGetColorMastervars(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets color for linking blocks in current color scheme
 * @returns linking color */
GCG_EXPORT
const char* SCIPvisuGetColorLinking(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets color for stairlinking blocks in current color scheme
 * @returns stairlinking color */
GCG_EXPORT
const char* SCIPvisuGetColorStairlinking(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets color for normal decomp blocks in current color scheme
 * @returns block color */
GCG_EXPORT
const char* SCIPvisuGetColorBlock(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets color for open blocks in current color scheme
 * @returns open color */
GCG_EXPORT
const char* SCIPvisuGetColorOpen(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets color for non-zero points in current color scheme
 * @returns non-zero color */
GCG_EXPORT
const char* SCIPvisuGetColorNonzero(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets color for lines in current color scheme
 * @returns line color */
GCG_EXPORT
const char* SCIPvisuGetColorLine(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets appropriate radius for nonzeros
 * needs highest indices of both axes for scaling
 * @returns radius */
GCG_EXPORT
float SCIPvisuGetNonzeroRadius(
   SCIP* scip,          /**< SCIP data structure */
   int maxindx,         /**< highest index x-axis */
   int maxindy,         /**< highest index y-axis */
   float scalingfactor  /**< percentage to scale radius, 1 if no scaling */
   );


/** if true gp reader should be used for sub-visualizations, otherwise tex reader
 * @returns true if gp reader should be used, false if tex reader should be used */
GCG_EXPORT
SCIP_Bool GCGgetUseGp(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets the name of the pdf reader that should be used
 * @returns name of pdf reader */
GCG_EXPORT
const char* GCGVisuGetPdfReader(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets the max number of decomps to be included in reports
 * @returns max number of decomps */
GCG_EXPORT
int GCGreportGetMaxNDecomps(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets what type of decomps to show in reports (where 0 corresponds to 'show all')
 * @returns type of decomps */
GCG_EXPORT
GCG_DECTYPE GCGreportGetDecompTypeToShow(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets whether a titlepage should be included in reports
 * @returns true iff title page should be generated */
GCG_EXPORT
SCIP_Bool GCGreportGetShowTitlepage(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets whether a table of contents should be included in reports
 * @returns true iff table of contents should be generated */
GCG_EXPORT
SCIP_Bool GCGreportGetShowToc(
   SCIP* scip  /**< SCIP data structure */
   );

/** gets whether statistics should be included for each decomp in reports
 * @returns true iff statistics for each decomp should be generated */
GCG_EXPORT
SCIP_Bool GCGreportGetShowStatistics(
   SCIP* scip  /**< SCIP data structure */
   );

/** frees all visualization parameters
 * 
 */
GCG_EXPORT
void GCGVisuFreeParams(
   SCIP* scip,                /**< SCIP data structure */
   GCG_PARAMDATA* paramdata   /**< input empty paramdata, oputput new set of param data */
   );

#ifdef __cplusplus
}
#endif

#endif
