/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2016 Operations Research, RWTH Aachen University       */
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
 * @brief  tex file reader for writing decomposition details to LaTeX files
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

#define COLOR_WHITE     "#FFFFFF"
#define COLOR_BLUE      "#00549F"
#define COLOR_LBLUE     "#8EBAE5"
#define COLOR_PURPLE    "#7A6FAC"
#define COLOR_VIOLET    "#612158"
#define COLOR_CARMINE   "#A11035"
#define COLOR_RED       "#CC071E"
#define COLOR_MAGENTA   "#E30066"
#define COLOR_ORANGE    "#F6A800"
#define COLOR_YELLOW    "#FFED00"
#define COLOR_GRASS     "#BDAB27"
#define COLOR_GREEN     "#57AB27"
#define COLOR_CYAN      "#0098A1"
#define COLOR_TEAL      "#006165"
#define COLOR_BLACK     "#000000"

#define DEFAULT_COLOR_MASTERVARS   COLOR_WHITE   /* for mastervars (in block area) */
#define DEFAULT_COLOR_MASTERCONS   COLOR_BLUE    /* for mastercons */
#define DEFAULT_COLOR_LINKING      COLOR_PURPLE
#define DEFAULT_COLOR_STAIRLINKING COLOR_MAGENTA
#define DEFAULT_COLOR_BLOCK        COLOR_TEAL
#define DEFAULT_COLOR_OPEN         COLOR_GREEN   /* for open (not assigned) elements */
#define DEFAULT_COLOR_NONZERO      COLOR_BLACK
#define DEFAULT_COLOR_LINE         COLOR_BLACK   /* for outlines of blocks */

enum Colorscheme
{
   COLORSCHEME_DEFAULT  = 0,     /**< default colors (supposedly eye-friendly) */
   COLORSCHEME_GREY     = 1,     /**< on a range from black to white */
   COLORSCHEME_MANUAL   = 2      /**< take user-defined input */
};

typedef enum Colorscheme VISU_COLORSCHEME; /**< visualization colorscheme type */

/** gets if draftmode is on
 * draftmode lets visualizations omit nonzeros */
extern
SCIP_Bool getDraftmode(void);

/** gets the colorscheme for visualizations */
extern
VISU_COLORSCHEME getColorscheme(void);

#ifdef __cplusplus
}
#endif

#endif
