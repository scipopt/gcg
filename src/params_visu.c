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

#define DEFAULT_VISU_RADIUS 5

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

char* greycolormastervars;
char* greycolormastercons;
char* greycolorlinking;
char* greycolorstairlinking;
char* greycolorblock;
char* greycoloropen;
char* greycolornonzero;
char* greycolorline;

int visuradius;

/*@todo vars for black and white scheme */

/** includes the visualization parameters into GCG */
SCIP_RETCODE SCIPincludeParamsVisu(
   SCIP* scip /**< SCIP data structure */
   )
{
   visudraftmode = FALSE;
   visucolorscheme = COLORSCHEME_DEFAULT;

   mancolormastervars =    (char*) DEFAULT_COLOR_MASTERVARS;
   mancolormastercons =    (char*) DEFAULT_COLOR_MASTERCONS;
   mancolorlinking =       (char*) DEFAULT_COLOR_LINKING;
   mancolorstairlinking =  (char*) DEFAULT_COLOR_STAIRLINKING;
   mancolorblock =         (char*) DEFAULT_COLOR_BLOCK;
   mancoloropen =          (char*) DEFAULT_COLOR_OPEN;
   mancolornonzero =       (char*) DEFAULT_COLOR_NONZERO;
   mancolorline =          (char*) DEFAULT_COLOR_LINE;

   /*@todo initialize black and white scheme*/

   visuradius = DEFAULT_VISU_RADIUS;

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
      return greycolormastercons;
      break;
   case COLORSCHEME_MANUAL:
      return mancolormastercons;
   default:
      return (char*) DEFAULT_COLOR_MASTERCONS;
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
