/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2019 Operations Research, RWTH Aachen University       */
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

#include "type_decomp.h"

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
extern SCIP_RETCODE SCIPincludeParamsVisu(
   SCIP* scip     /**< SCIP data structure */
   );

/** gets whether draftmode is on
 * draftmode lets visualizations omit nonzeros
 * @returns true if draftmode is on */
extern SCIP_Bool SCIPvisuGetDraftmode(void);

/** sets draftmode
 * draftmode lets visualizations omit nonzeros
 * @returns nothing */
extern void SCIPvisuSetDraftmode(
   SCIP_Bool setmode /**< true iff draftmode should be on */
   );

/** gets the colorscheme for visualizations
 * @returns current colorscheme */
extern VISU_COLORSCHEME SCIPvisuGetColorscheme(void);

/** sets colorscheme for visualizations
 * @returns nothing */
extern void SCIPvisuSetColorscheme(
   VISU_COLORSCHEME newscheme    /**< new colorscheme */
   );

/** sets color for mastercons block in current color scheme
 *  @returns nothing*/
extern void SCIPvisuSetColorManMasterconss(
   char* newcolor       /**< new color */
   );

/** sets manual color for mastervar block in current color scheme
 * @returns nothing */
extern void SCIPvisuSetColorManMastervars(
   char* newcolor       /**< new color */
   );

/** sets manual color for linking blocks in current color scheme
 * @returns nothing */
extern void SCIPvisuSetColorManLinking(
   char* newcolor       /**< new color */
   );

/** sets manual color for stairlinking blocks in current color scheme
 * @returns nothing */
extern void SCIPvisuSetColorManStairlinking(
   char* newcolor       /**< new color */
   );

/** sets manual color for normal decomp blocks in current color scheme
 * @returns nothing */
extern void SCIPvisuSetColorManBlock(
   char* newcolor       /**< new color */
   );

/** sets manual color for open blocks in current color scheme
 * @returns nothing */
extern void SCIPvisuSetColorManOpen(
   char* newcolor       /**< new color */
   );

/** sets manual color for non-zero points in current color scheme
 * @returns nothing */
extern void SCIPvisuSetColorManNonzero(
   char* newcolor       /**< new color */
   );

/** sets manual color for lines in current color scheme
 * @returns nothing */
extern void SCIPvisuSetColorManLine(
   char* newcolor       /**< new color */
   );

/** gets color for mastercons block in current color scheme
 * @returns mastercons color */
extern char* SCIPvisuGetColorMasterconss(void);

/** gets color for mastervar block in current color scheme
 * @returns mastervar color */
extern char* SCIPvisuGetColorMastervars(void);

/** gets color for linking blocks in current color scheme
 * @returns linking color */
extern char* SCIPvisuGetColorLinking(void);

/** gets color for stairlinking blocks in current color scheme
 * @returns stairlinking color */
extern char* SCIPvisuGetColorStairlinking(void);

/** gets color for normal decomp blocks in current color scheme
 * @returns block color */
extern char* SCIPvisuGetColorBlock(void);

/** gets color for open blocks in current color scheme
 * @returns open color */
extern char* SCIPvisuGetColorOpen(void);

/** gets color for non-zero points in current color scheme
 * @returns non-zero color */
extern char* SCIPvisuGetColorNonzero(void);

/** gets color for lines in current color scheme
 * @returns line color */
extern char* SCIPvisuGetColorLine(void);

/** gets appropriate radius for nonzeros
 * needs highest indices of both axes for scaling
 * @returns radius */
extern float SCIPvisuGetNonzeroRadius(
   int maxindx,         /**< highest index x-axis */
   int maxindy,         /**< highest index y-axis */
   float scalingfactor  /**< percentage to scale radius, 1 if no scaling */
   );

/** if true gp reader should be used for sub-visualizations, otherwise tex reader
 * @returns true if gp reader should be used, false if tex reader should be used */
extern SCIP_Bool GCGgetUseGp(void);

/** gets the name of the pdf reader that should be used
 * @returns name of pdf reader */
extern char* GCGVisuGetPdfReader(void);

/** gets the max number of decomps to be included in reports
 * @returns max number of decomps */
extern int GCGreportGetMaxNDecomps(void);

/** gets what type of decomps to show in reports (where 0 corresponds to 'show all')
 * @returns type of decomps */
extern DEC_DECTYPE GCGreportGetDecompTypeToShow(void);

/** gets whether a titlepage should be included in reports
 * @returns true iff title page should be generated */
extern SCIP_Bool GCGreportGetShowTitlepage(void);

/** gets whether a table of contents should be included in reports
 * @returns true iff table of contents should be generated */
extern SCIP_Bool GCGreportGetShowToc(void);

/** gets whether statistics should be included for each decomp in reports
 * @returns true iff statistics for each decomp should be generated */
extern SCIP_Bool GCGreportGetShowStatistics(void);

/** gets the max number of finished decomps to be included in family tree
 * @returns max number of finished decomps */
extern int GCGfamtreeGetMaxNDecomps(void);

/*@todo include this somewhere where it is automatically called at quitting GCG */
/** frees all visualization parameters
 * @returns nothing */
extern void GCGVisuFreeParams(
   SCIP* scip     /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
