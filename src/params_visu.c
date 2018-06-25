/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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
 * @author  Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "params_visu.h"
#include "type_decomp.h"
#include "cons_decomp.h"

#include <limits.h>


/* color defaults to build default color layout with */
#define COLOR_WHITE     "#FFFFFF"   /**< standard white */
#define COLOR_BLUE1     "#ACBCE9"   /**< very light blue */
#define COLOR_BLUE2     "#718CDB"   /**< light blue */
#define COLOR_BLUE3     "#3C64DD"   /**< middle blue */
#define COLOR_BLUE4     "#1340C7"   /**< dark blue */
#define COLOR_BLUE5     "#1F377D"   /**< very dark blue */
#define COLOR_ORANGE1   "#FFD88F"   /**< very light orange */
#define COLOR_ORANGE2   "#FFCB69"   /**< light orange */
#define COLOR_ORANGE3   "#FFB72D"   /**< orange */
#define COLOR_BROWN1    "#B38208"   /**< light brown */
#define COLOR_BROWN2    "#886100"   /**< brown */
#define COLOR_BROWN3    "#443000"   /**< dark brown */
#define COLOR_BLACK     "#000000"   /**< standard black */

/* default colors (use defines above for changes) */
#define DEFAULT_COLOR_MASTERVARS   COLOR_BLUE4     /**< for mastervars (in block area) */
#define DEFAULT_COLOR_MASTERCONSS  COLOR_BLUE4     /**< for masterconss */
#define DEFAULT_COLOR_LINKING      COLOR_ORANGE3   /**< for linking areas */
#define DEFAULT_COLOR_STAIRLINKING COLOR_BROWN2    /**< for stairlinking areas */
#define DEFAULT_COLOR_BLOCK        COLOR_BLUE2     /**< for finished blocks */
#define DEFAULT_COLOR_OPEN         COLOR_ORANGE1   /**< for open (not assigned) elements */
#define DEFAULT_COLOR_NONZERO      COLOR_BLACK     /**< for nonzero dots */
#define DEFAULT_COLOR_LINE         COLOR_BLACK     /**< for outlines of blocks */

/* 8 shades of grey */
#define GREY_COLOR_MASTERVARS   "#323232"    /**< for mastervars (in block area) */
#define GREY_COLOR_MASTERCONS   "#666666"    /**< for masterconss */
#define GREY_COLOR_LINKING      "#4C4C4C"    /**< for linking areas */
#define GREY_COLOR_STAIRLINKING "#191919"    /**< for stairlinking areas */
#define GREY_COLOR_BLOCK        "#d3d3d3"    /**< for finished blocks */
#define GREY_COLOR_OPEN         "#7F7F7F"    /**< for open (not assigned) elements */
#define GREY_COLOR_NONZERO      COLOR_BLACK  /**< for nonzero dots */
#define GREY_COLOR_LINE         COLOR_BLACK  /**< for outlines of blocks */

/* visualization imaging defaults */
#define DEFAULT_VISU_DRAFTMODE   FALSE                /**< if true no nonzeros are shown in visualizations */
#define DEFAULT_VISU_COLORSCHEME COLORSCHEME_DEFAULT  /**< is of type VISU_COLORSCHEME */
#define DEFAULT_VISU_RADIUS      2                    /**< possible scale: 1-10 */
#define DEFAULT_VISU_USEGP       FALSE                /**< if true gnuplot is used for visualizations,
                                                       * otherwise LaTeX/Tikz */

/* pdf reader default */
#define DEFAULT_PDFREADER        "evince"             /**< name of pdf reader, must be callable by system */

/* report parameter defaults */
#define DEFAULT_REPORT_MAXNDECOMPS     20       /**< maximum number of decomps to be shown in report */
#define DEFAULT_REPORT_SHOWTYPE        0        /**< what type of decomps to show
                                                 * (DEC_DECTYPE, but 0 corresponds to 'show all') */
#define DEFAULT_REPORT_SHOWTITLEPAGE   TRUE     /**< if true a titlepage is included */
#define DEFAULT_REPORT_SHOWTOC         TRUE     /**< if true a table of contents is included */
#define DEFAULT_REPORT_SHOWSTATISTICS  TRUE     /**< if true statistics are included for each decomp */

/* familytree parameter defaults */
#define DEFAULT_FAMTREE_MAXNDECOMPS    5        /**< maximum number of finished decompositions in family tree */


struct GCG_VisualizationData
{
   SCIP_Bool         visudraftmode;    /**< true if no nonzeros should be shown */
   VISU_COLORSCHEME  visucolorscheme;  /**< stores the current color scheme */
   int               visuradius;       /**< radius for nonzeros */
   SCIP_Bool         visuusegp;        /**< if true gnuplot is used for visualizations, otherwise LaTeX/Tikz */

   char* mancolormastervars;           /**< manual color for master variables */
   char* mancolormasterconss;          /**< manual color for master constraints */
   char* mancolorlinking;              /**< manual color for linking */
   char* mancolorstairlinking;         /**< manual color for stairlinking */
   char* mancolorblock;                /**< manual color for blocks */
   char* mancoloropen;                 /**< manual color for nonassigned areas */
   char* mancolornonzero;              /**< manual color for nonzeros */
   char* mancolorline;                 /**< manual color for lines */

   char* greycolormastervars;          /**< black and white color for master variables */
   char* greycolormasterconss;         /**< black and white color for master constraints */
   char* greycolorlinking;             /**< black and white color for linking */
   char* greycolorstairlinking;        /**< black and white color for stairlinking */
   char* greycolorblock;               /**< black and white color for blocks */
   char* greycoloropen;                /**< black and white color for nonassigned areas */
   char* greycolornonzero;             /**< black and white color for nonzeros */
   char* greycolorline;                /**< black and white color for lines */

   char* pdfreader;                    /**< name of pdfreader to open files with */

   int         rep_maxndecomps;        /**< maximum number of decomps to be shown in report */
   DEC_DECTYPE rep_showtype;           /**< what type of decomps to show (where 0 corresponds to 'show all') */
   SCIP_Bool   rep_showtitle;          /**< if true a titlepage is included */
   SCIP_Bool   rep_showtoc;            /**< if true a table of contents is included */
   SCIP_Bool   rep_statistics;         /**< if true statistics are included for each decomp */

   int         fam_maxndecomps;        /**< maximum number of finished decompositions in family tree */
   int         nmaxdecompstowrite;     /**< maximum number of decompositions to write */
};

/** visualization parameter data */
struct GCG_VisualizationData* visudata;


/** includes the visualization parameters into GCG */
SCIP_RETCODE SCIPincludeParamsVisu(
   SCIP* scip /**< SCIP data structure */
   )
{
   visudata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &visudata) );

   /* init string params with NULL pointer */
   visudata->pdfreader = NULL;
   visudata->mancolormastervars = NULL;
   visudata->mancolormasterconss = NULL;
   visudata->mancolorlinking = NULL;
   visudata->mancolorstairlinking = NULL;
   visudata->mancolorblock = NULL;
   visudata->mancoloropen = NULL;
   visudata->mancolornonzero = NULL;
   visudata->mancolorline = NULL;

   /* add general parameters */

   SCIP_CALL( SCIPaddBoolParam(scip,
      "visual/draftmode", "if true no nonzeros are shown (may improve performance)",
      &visudata->visudraftmode, FALSE, DEFAULT_VISU_DRAFTMODE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "visual/colorscheme", "type number: 0=default, 1=black and white, 2=manual",
      (int*) &visudata->visucolorscheme, FALSE, DEFAULT_VISU_COLORSCHEME, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "visual/nonzeroradius", "integer value to scale points on range 1-10",
      &visudata->visuradius, FALSE, DEFAULT_VISU_RADIUS, 1, 10, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "visual/nmaxdecompstowrite", "maximum number of decompositions to write (-1: no limit)",
      &visudata->nmaxdecompstowrite, FALSE, -1, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/pdfreader", "pdf reader that opens visualizations in decomposition explorer",
      &visudata->pdfreader, FALSE,
      DEFAULT_PDFREADER,
      NULL, NULL) );

   /* add parameters for manual colors */

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/colormastervars", "color for master variables in hex code",
      &visudata->mancolormastervars, FALSE, DEFAULT_COLOR_MASTERVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/colormasterconss", "color for master constraints in hex code",
      &visudata->mancolormasterconss, FALSE, DEFAULT_COLOR_MASTERCONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/colorlinking", "color for linking variables in hex code",
      &visudata->mancolorlinking, FALSE, DEFAULT_COLOR_LINKING, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/colorstairlinking", "color for stairlinking variables in hex code",
      &visudata->mancolorstairlinking, FALSE, DEFAULT_COLOR_STAIRLINKING, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/colorblock", "color for found blocks in hex code",
      &visudata->mancolorblock, FALSE, DEFAULT_COLOR_BLOCK, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/coloropen", "color for open areas in hex code",
      &visudata->mancoloropen, FALSE, DEFAULT_COLOR_OPEN, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/colornonzeros", "color for nonzeros in hex code",
      &visudata->mancolornonzero, FALSE, DEFAULT_COLOR_NONZERO, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
      "visual/colors/colorlines", "color for lines in hex code",
      &visudata->mancolorline, FALSE, DEFAULT_COLOR_LINE, NULL, NULL) );

   /* add parameters for report */

   SCIP_CALL( SCIPaddIntParam(scip,
      "visual/report/maxndecomps", "maximum number of decompositions shown in report (best scores first)",
      &visudata->rep_maxndecomps, FALSE, DEFAULT_REPORT_MAXNDECOMPS, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "visual/report/showtype",
      "only decompositions of type: 0=all types, 1=arrowhead, 2=staircase, 3=diagonal, 4=bordered",
      (int*) &visudata->rep_showtype, FALSE, DEFAULT_REPORT_SHOWTYPE, 0, 4, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "visual/report/showtitle", "if true a title page is included",
      &visudata->rep_showtitle, FALSE, DEFAULT_REPORT_SHOWTITLEPAGE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "visual/report/showtoc", "if true a table of contents is included",
      &visudata->rep_showtoc, FALSE, DEFAULT_REPORT_SHOWTOC, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "visual/report/showstatistics", "if true statistics are included for each decomp",
      &visudata->rep_statistics, FALSE, DEFAULT_REPORT_SHOWSTATISTICS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "visual/report/usegp", "if true gnuplot is used for sub-visualizations in report, otherwise LaTeX/Tikz",
      &visudata->visuusegp, FALSE, DEFAULT_VISU_USEGP, NULL, NULL) );

   /* add parameters for family tree */

   SCIP_CALL( SCIPaddIntParam(scip,
      "visual/famtree/maxndecomps", "maximum number of finished decompositions in family tree",
      &visudata->fam_maxndecomps, FALSE, DEFAULT_FAMTREE_MAXNDECOMPS, 1, INT_MAX, NULL, NULL) );

   /* initialize black and white color scheme */

   visudata->greycolormastervars = (char*) GREY_COLOR_MASTERVARS;
   visudata->greycolormasterconss = (char*) GREY_COLOR_MASTERCONS;
   visudata->greycolorlinking = (char*) GREY_COLOR_LINKING;
   visudata->greycolorstairlinking = (char*) GREY_COLOR_STAIRLINKING;
   visudata->greycolorblock = (char*) GREY_COLOR_BLOCK;
   visudata->greycoloropen = (char*) GREY_COLOR_OPEN;
   visudata->greycolornonzero = (char*) GREY_COLOR_NONZERO;
   visudata->greycolorline = (char*) GREY_COLOR_LINE;

   return SCIP_OKAY;
}

/* getter & setter */


/** gets if draftmode is on
 * draftmode lets visualizations omit nonzeros */
SCIP_Bool SCIPvisuGetDraftmode()
{
   return visudata->visudraftmode;
}


/** sets draftmode
 * draftmode lets visualizations omit nonzeros */
void SCIPvisuSetDraftmode(
   SCIP_Bool setmode
   )
{
   visudata->visudraftmode = setmode;
}


/** gets the colorscheme for visualizations */
VISU_COLORSCHEME SCIPvisuGetColorscheme()
{
   return visudata->visucolorscheme;
}


/** sets colorscheme for visualizations */
void SCIPvisuSetColorscheme(
   VISU_COLORSCHEME newscheme
   )
{
   visudata->visucolorscheme = newscheme;
}


/** gets color for mastercon block in current color scheme */
char* SCIPvisuGetColorMasterconss()
{
   switch(SCIPvisuGetColorscheme())
   {
   case COLORSCHEME_GREY:
      return visudata->greycolormasterconss;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancolormasterconss;
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
      return visudata->greycolormastervars;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancolormastervars;
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
      return visudata->greycolorlinking;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancolorlinking;
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
      return visudata->greycolorstairlinking;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancolorstairlinking;
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
      return visudata->greycolorblock;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancolorblock;
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
      return visudata->greycoloropen;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancoloropen;
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
      return visudata->greycolornonzero;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancolornonzero;
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
      return visudata->greycolorline;
      break;
   case COLORSCHEME_MANUAL:
      return visudata->mancolorline;
   default:
      return (char*) DEFAULT_COLOR_LINE;
   }
}


/** sets color for mastercon block in current color scheme */
void SCIPvisuSetColorManMasterconss(
   char* newcolor       /**< new color */
   )
{
   visudata->mancolormasterconss = newcolor;
}


/** sets manual color for mastervar block in current color scheme */
void SCIPvisuSetColorManMastervars(
   char* newcolor       /**< new color */
   )
{
   visudata->mancolormastervars = newcolor;
}


/** sets manual color for linking blocks in current color scheme */
void SCIPvisuSetColorManLinking(
   char* newcolor       /**< new color */
   )
{
   visudata->mancolorlinking = newcolor;
}


/** sets manual color for stairlinking blocks in current color scheme */
void SCIPvisuSetColorManStairlinking(
   char* newcolor       /**< new color */
   )
{
   visudata->mancolorstairlinking = newcolor;
}


/** sets manual color for normal decomp blocks in current color scheme */
void SCIPvisuSetColorManBlock(
   char* newcolor       /**< new color */
   )
{
   visudata->mancolorblock = newcolor;
}


/** sets manual color for open blocks in current color scheme */
void SCIPvisuSetColorManOpen(
   char* newcolor       /**< new color */
   )
{
   visudata->mancoloropen = newcolor;
}


/** sets manual color for non-zero points in current color scheme */
void SCIPvisuSetColorManNonzero(
   char* newcolor       /**< new color */
   )
{
   visudata->mancolornonzero = newcolor;
}


/** sets manual color for lines in current color scheme */
void SCIPvisuSetColorManLine(
   char* newcolor       /**< new color */
   )
{
   visudata->mancolorline = newcolor;
}


/** gets appropriate radius for nonzeros
 * needs highest indices of both axes */
float SCIPvisuGetNonzeroRadius(
   int maxindx,         /**< highest index x-axis */
   int maxindy,         /**< highest index y-axis */
   float scalingfactor  /**< percentage to scale radius, 1 if no scaling */
   )
{
   int maxind = 0;

   /* the max indices must be at least one to be compatible with division */
   if(maxindx <= 0)
      maxindx = 1;

   if(maxindy <= 0)
      maxindy = 1;

   /* determine the highest index */
   if(maxindx > maxindy)
      maxind = maxindx;
   else
      maxind = maxindy;

   /* scale by coordinate system size and given factor */
   return ( (float) visudata->visuradius / (float) maxind) * scalingfactor;
}


/** if true gp reader should be used for sub-visualizations, otherwise tex reader */
SCIP_Bool GCGgetUseGp()
{
   return visudata->visuusegp;
}


/** gets the name of the pdf reader that should be used */
char* GCGVisuGetPdfReader()
{
   return visudata->pdfreader;
}


/** gets the max number of decomps to be included in reports */
int GCGreportGetMaxNDecomps()
{
   return visudata->rep_maxndecomps;
}


/** gets what type of decomps to show in reports (where 0 corresponds to 'show all') */
DEC_DECTYPE GCGreportGetDecompTypeToShow()
{
   return visudata->rep_showtype;
}


/** gets whether a titlepage should be included in reports */
SCIP_Bool GCGreportGetShowTitlepage()
{
   return visudata->rep_showtitle;
}


/** gets whether a table of contents should be included in reports */
SCIP_Bool GCGreportGetShowToc()
{
   return visudata->rep_showtoc;
}


/** gets whether statistics should be included for each decomp in reports */
SCIP_Bool GCGreportGetShowStatistics()
{
   return visudata->rep_statistics;
}


/** gets the max number of finished decomps to be included in family tree */
int GCGfamtreeGetMaxNDecomps()
{
   return visudata->fam_maxndecomps;
}


/** frees all visualization parameters */
void GCGVisuFreeParams(
   SCIP* scip     /**< SCIP data structure */
   )
{
   if ( visudata != NULL )
      SCIPfreeMemory(scip, &visudata);
}
