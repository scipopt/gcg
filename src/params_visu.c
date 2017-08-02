/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
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

/**@file    params_visu.c
 * @brief   parameter-related stuff for visualization
 * @author  Hanna Franzen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "params_visu.h"

/* global visualization parameters */

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

/*@todo defines for black and white color scheme */

SCIP_Bool visudraftmode;
VISU_COLORSCHEME visucolorscheme;

char* mancolormastervars;
char* mancolormastercons;
char* mancolorlinking;
char* mancolorstairlinking;
char* mancolorblock;
char* mancoloropen;
char* mancolornonzero;
char* mancolorline;

/*@todo vars for black and white scheme */

/** includes the visualization parameters into GCG */
SCIP_RETCODE SCIPincludeParamsVisu(
   SCIP* scip /**< SCIP data structure */
   )
{
   visudraftmode = FALSE;
   visucolorscheme = COLORSCHEME_DEFAULT;

   mancolormastervars = DEFAULT_COLOR_MASTERVARS;
   mancolormastercons = DEFAULT_COLOR_MASTERCONS;
   mancolorlinking = DEFAULT_COLOR_LINKING;
   mancolorstairlinking = DEFAULT_COLOR_STAIRLINKING;
   mancolorblock = DEFAULT_COLOR_BLOCK;
   mancoloropen = DEFAULT_COLOR_OPEN;
   mancolornonzero = DEFAULT_COLOR_NONZERO;
   mancolorline = DEFAULT_COLOR_LINE;

   /*@todo initialize black and white scheme*/
}

/* getter & setter */

/** gets if draftmode is on
 * draftmode lets visualizations omit nonzeros */
SCIP_Bool SCIPvisuGetDraftmode()
{
   return visudraftmode;
}

/** sets draftmode
 * draftmode lets visualizations omit nonzeros */
void SCIPvisuSetDraftmode(SCIP_Bool setmode)
{
   visudraftmode = setmode;
}

/** gets the colorscheme for visualizations */
VISU_COLORSCHEME SCIPvisuGetColorscheme()
{
   return visucolorscheme;
}

/** sets colorscheme for visualizations */
void SCIPvisuSetColorscheme(VISU_COLORSCHEME newscheme)
{
   visucolorscheme = newscheme;
}

/*@todo (getter &) setter for manual color scheme*/
