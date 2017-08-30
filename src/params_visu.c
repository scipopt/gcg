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
#define DEFAULT_COLOR_MASTERCONSS  COLOR_BLUE    /* for masterconss */
#define DEFAULT_COLOR_LINKING      COLOR_PURPLE
#define DEFAULT_COLOR_STAIRLINKING COLOR_MAGENTA
#define DEFAULT_COLOR_BLOCK        COLOR_TEAL
#define DEFAULT_COLOR_OPEN         COLOR_GREEN   /* for open (not assigned) elements */
#define DEFAULT_COLOR_NONZERO      COLOR_BLACK
#define DEFAULT_COLOR_LINE         COLOR_BLACK   /* for outlines of blocks */

#define GREY_COLOR_MASTERVARS   "#323232"
#define GREY_COLOR_MASTERCONS   "#999999"
#define GREY_COLOR_LINKING      "#666666"
#define GREY_COLOR_STAIRLINKING "#191919"
#define GREY_COLOR_BLOCK        "#4C4C4C"
#define GREY_COLOR_OPEN         "#7F7F7F"
#define GREY_COLOR_NONZERO      COLOR_BLACK
#define GREY_COLOR_LINE         COLOR_BLACK

#define DEFAULT_VISU_RADIUS 5    /* possible scale: 1-10 */


SCIP_Bool visudraftmode;
VISU_COLORSCHEME visucolorscheme;

char* mancolormastervars;
char* mancolormasterconss;
char* mancolorlinking;
char* mancolorstairlinking;
char* mancolorblock;
char* mancoloropen;
char* mancolornonzero;
char* mancolorline;

char* greycolormastervars;
char* greycolormasterconss;
char* greycolorlinking;
char* greycolorstairlinking;
char* greycolorblock;
char* greycoloropen;
char* greycolornonzero;
char* greycolorline;

int visuradius;


/** includes the visualization parameters into GCG */
SCIP_RETCODE SCIPincludeParamsVisu(
   SCIP* scip /**< SCIP data structure */
   )
{
   /* set defaults */
   visudraftmode = FALSE;
   visucolorscheme = COLORSCHEME_DEFAULT;

   mancolormastervars =    (char*) DEFAULT_COLOR_MASTERVARS;
   mancolormasterconss =   (char*) DEFAULT_COLOR_MASTERCONSS;
   mancolorlinking =       (char*) DEFAULT_COLOR_LINKING;
   mancolorstairlinking =  (char*) DEFAULT_COLOR_STAIRLINKING;
   mancolorblock =         (char*) DEFAULT_COLOR_BLOCK;
   mancoloropen =          (char*) DEFAULT_COLOR_OPEN;
   mancolornonzero =       (char*) DEFAULT_COLOR_NONZERO;
   mancolorline =          (char*) DEFAULT_COLOR_LINE;

   greycolormastervars =   (char*) GREY_COLOR_MASTERVARS;
   greycolormasterconss =  (char*) GREY_COLOR_MASTERCONS;
   greycolorlinking =      (char*) GREY_COLOR_LINKING;
   greycolorstairlinking = (char*) GREY_COLOR_STAIRLINKING;
   greycolorblock =        (char*) GREY_COLOR_BLOCK;
   greycoloropen =         (char*) GREY_COLOR_OPEN;
   greycolornonzero =      (char*) GREY_COLOR_NONZERO;
   greycolorline =         (char*) GREY_COLOR_LINE;

   visuradius = DEFAULT_VISU_RADIUS;

   /* add general parameters */

   SCIP_CALL( SCIPaddBoolParam(scip,
      "visualization/draftmode", "if true no nonzeros are shown (may improve performance)",
      &visudraftmode, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "visualization/colorscheme", "type number: 0=default, 1=black and white, 2=manual",
      (int*) &visucolorscheme, FALSE, 0, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "visualization/nonzeroradius", "integer value to scale dots from 1-10, default: 5",
      &visuradius, FALSE, DEFAULT_VISU_RADIUS, 1, 10, NULL, NULL) );

   /* add parameters for manual colors */

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/colormastervars", "color for master variables in hex code (e.g. #000000)",
      &mancolormastervars, FALSE, DEFAULT_COLOR_MASTERVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/colormasterconss", "color for master constraints in hex code (e.g. #000000)",
      &mancolormasterconss, FALSE, DEFAULT_COLOR_MASTERCONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/colorlinking", "color for linking variables in hex code (e.g. #000000)",
      &mancolorlinking, FALSE, DEFAULT_COLOR_LINKING, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/colorstairlinking", "color for stairlinking variables in hex code (e.g. #000000)",
      &mancolorstairlinking, FALSE, DEFAULT_COLOR_STAIRLINKING, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/colorblock", "color for found blocks in hex code (e.g. #000000)",
      &mancolorblock, FALSE, DEFAULT_COLOR_BLOCK, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/coloropen", "color for open areas in hex code (e.g. #000000)",
      &mancoloropen, FALSE, DEFAULT_COLOR_OPEN, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/colornonzeros", "color for nonzeros in hex code (e.g. #000000)",
      &mancolornonzero, FALSE, DEFAULT_COLOR_NONZERO, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visualization/colors/colorlines", "color for lines in hex code (e.g. #000000)",
      &mancolorline, FALSE, DEFAULT_COLOR_LINE, NULL, NULL) );

   return SCIP_OKAY;
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


/*@todo setter for manual color scheme*/


/** gets color for mastercon block in current color scheme */
char* SCIPvisuGetColorMasterconss()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycolormasterconss;
      break;
   case COLORSCHEME_MANUAL:
      return mancolormasterconss;
   default:
      return (char*) DEFAULT_COLOR_MASTERCONSS;
   }
}

/** gets color for mastervar block in current color scheme */
char* SCIPvisuGetColorMastervars()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycolormastervars;
      break;
   case COLORSCHEME_MANUAL:
      return mancolormastervars;
   default:
      return (char*) DEFAULT_COLOR_MASTERVARS;
   }
}

/** gets color for linking blocks in current color scheme */
char* SCIPvisuGetColorLinking()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycolorlinking;
      break;
   case COLORSCHEME_MANUAL:
      return mancolorlinking;
   default:
      return (char*) DEFAULT_COLOR_LINKING;
   }
}

/** gets color for stairlinking blocks in current color scheme */
char* SCIPvisuGetColorStairlinking()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycolorstairlinking;
      break;
   case COLORSCHEME_MANUAL:
      return mancolorstairlinking;
   default:
      return (char*) DEFAULT_COLOR_STAIRLINKING;
   }
}

/** gets color for normal decomp blocks in current color scheme */
char* SCIPvisuGetColorBlock()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycolorblock;
      break;
   case COLORSCHEME_MANUAL:
      return mancolorblock;
   default:
      return (char*) DEFAULT_COLOR_BLOCK;
   }
}

/** gets color for open blocks in current color scheme */
char* SCIPvisuGetColorOpen()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycoloropen;
      break;
   case COLORSCHEME_MANUAL:
      return mancoloropen;
   default:
      return (char*) DEFAULT_COLOR_OPEN;
   }
}

/** gets color for non-zero points in current color scheme */
char* SCIPvisuGetColorNonzero()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycolornonzero;
      break;
   case COLORSCHEME_MANUAL:
      return mancolornonzero;
   default:
      return (char*) DEFAULT_COLOR_NONZERO;
   }
}

/** gets color for lines in current color scheme */
char* SCIPvisuGetColorLine()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return greycolorline;
      break;
   case COLORSCHEME_MANUAL:
      return mancolorline;
   default:
      return (char*) DEFAULT_COLOR_LINE;
   }
}

/** gets appropriate radius for nonzeros
 * needs highest indices of both axes */
float SCIPvisuGetNonzeroRadius(
   int maxindx,     /**< highest index x-axis */
   int maxindy,    /**< highest index y-axis */
   float scalingfactor /**< percentage to scale radius, 1 if no scaling */
   )
{
   int maxind = 0;

   /* the max indices must be at least one to be compatible with division */
   maxindx = maxindx < 1 ? 1 : maxindx;
   maxindy = maxindy < 1 ? 1 : maxindy;

   /* determine the highest index */
   maxind = maxindx>maxindy?maxindx:maxindy;

   /* scale by coordinate system size and given factor */
   return (visuradius / maxind) * scalingfactor;
}
