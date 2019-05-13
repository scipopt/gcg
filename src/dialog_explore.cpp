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

/**@file   dialog_explore.c
 * @brief  dialog menu for exploring decompositions
 * @author Michael Bastubbe
 * @author Hanna Franzen
 *
 * This file contains all dialog calls to build and use the explore menu.
 * The explore menu gives the user detailed information about all decompositions and a possibility to edit such.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <stddef.h>

#include <string>
#include <iostream>
#include <regex>
#include <map>

#include "class_seeed.h"
#include "cons_decomp.h"
#include "wrapper_seeed.h"

/* column headers */
#define DEFAULT_COLUMN_MIN_WIDTH  4 /**< min width of a column in the menu table */
#define DEFAULT_COLUMN_MAX_WIDTH 10 /**< max width of a column (also determines max width of column header abbreviation) */
#define DEFAULT_COLUMNS "nr id nbloc nmacon nlivar nmavar nstlva score history pre nopcon nopvar sel" /**< default column headers */

#define DESC_NR      "number of the decomposition (use this number for choosing the decomposition)" /**< description of number */
#define DESC_ID      "id of the decomposition (identifies the decomposition in reports/statistics/visualizations/etc.)" /**< description of seeed id */
#define DESC_NBLOC   "number of blocks" /**< description of block number */
#define DESC_NMACON  "number of master constraints" /**< description of master conss number */
#define DESC_NLIVAR  "number of linking variables" /**< description of linking vars number */
#define DESC_NMAVAR  "number of master variables (do not occur in blocks)" /**< description of master vars number */
#define DESC_NSTLVA  "number of stairlinking variables (disjoint from linking variables)" /**< description of stairlinking vars number */

#define DESC_SCORE   " " //@todo put this back in the actual scores, they should know their description 

#define DESC_HISTORY "list of detector chars worked on this decomposition" /**< description of detection history */
#define DESC_PRE     "is this decomposition for the presolved problem" /**< description of presolved bool */
#define DESC_NOPCON  "number of open constraints" /**< description of open conss number */
#define DESC_NOPVAR  "number of open variables" /**< description of open vars number */
#define DESC_USR     "whether this decomposition was given by the user" /**< description of user given bool */
#define DESC_SEL     "is this decomposition selected at the moment" /**< description of selected bool */

#define DEFAULT_MENULENGTH 10 /**< initial number of entries in menu */

namespace gcg
{


/** modifies menulength according to input and updates menu accordingly
 * @returns SCIP return code */
static
SCIP_RETCODE SCIPdialogSetNEntires(
   SCIP* scip,                   /**< SCIP data structure */
   SCIP_DIALOGHDLR* dialoghdlr,  /**< dialog handler for user input management */
   SCIP_DIALOG* dialog,          /**< dialog for user input management */
   int* menulength               /**< current menu length to be modified */
   )
{
   char* ntovisualize;
   SCIP_Bool endoffile;
   int newlength;
   int commandlen;

   SCIPdialogMessage(scip, NULL, "Please specify the amount of entries to be shown in this menu:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   newlength = -1;
   if( commandlen != 0 )
      newlength = atoi(ntovisualize);

   /* check whether there are decompositions,
    * (preventing "Why doesn't it show anything? Maybe the entry number is 0") */
   if( SCIPconshdlrDecompGetNSeeeds(scip) == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }

   if( commandlen == 0 || newlength < 1 )
   {
      SCIPdialogMessage( scip, NULL, "The input was not a valid number." );
      return SCIP_OKAY;
   }

   *menulength = newlength;

   return SCIP_OKAY;
}


/**
 * @brief method to update the list of incomplete decompositions
 *
 * this list changes due to new decompositions, modified, decompositions or changes of the score
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPdialogUpdateSeeedlist(
   SCIP* scip,       /**< SCIP data structure */
   int* startindex   /**< start index of menu */
   )
{
   assert( SCIPconshdlrDecompCheckConsistency(scip) );

   *startindex = 0;

   SCIPconshdlrdataDecompUnselectAll(scip);

   /* sort decomposition and finished seeeds according to max white score */
   /*@todo remove this when manual sorting in menu is implemented */
   SCIP_CALL( DECconshdlrDecompSortDecompositionsByScore(scip) );

   return SCIP_OKAY;
}

/**
 * Outputs the given char x times as SCIPdialogMessage
 * @returns SCIP status
 */
static
SCIP_RETCODE outputCharXTimes(
   SCIP* scip,          /**< SCIP data structure */
   const char letter,   /**< char to write */
   int x                /**< write char x times */
   )
{
   for(int i = 0; i < x; i++)
      SCIPdialogMessage(scip, NULL, "%c", letter);

   return SCIP_OKAY;
}

/** @brief show current menu containing seeed information
 *
 * Update length of seeed list in case it changed since the last command
 * and show the table of seeeds. 
 * @returns SCIP status
 */
static
SCIP_RETCODE SCIPdialogShowMenu(
   SCIP* scip,                         /**< SCIP data structure */
   std::vector<std::string> columns,   /**< list of column headers (abbreviations) */
   int* nseeeds,                       /**< max number of seeeds */
   const int startindex,               /**< index (in seeed list) of uppermost seeed in extract */
   int menulength,                     /**< number of menu entries */
   int** idlist,                       /**< current list of seeed ids */
   int*  listlength                    /**< length of idlist */
   )
{
   assert(scip != NULL);

   /* update size of seeed list in case it changed */
   if(*nseeeds < SCIPconshdlrDecompGetNSeeeds(scip))
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, idlist, *nseeeds, SCIPconshdlrDecompGetNSeeeds(scip)) );
      *nseeeds = SCIPconshdlrDecompGetNSeeeds(scip);
   }

   /* update list of seeeds */
   SCIPconshdlrDecompGetSeeedLeafList(scip, idlist, listlength);

   /* count corresponding seeeds for overview statistics */
   int ndetectedpresolved = 0;
   int ndetectedunpresolved = 0;

   for(int i = 0; i < *listlength; ++i)
   {
      Seeed_Wrapper sw;
      GCGgetSeeedFromID(scip, &(*idlist)[i], &sw);
      Seeed* seeed;
      seeed = sw.seeed;
      /* finished seeeds */
      if(seeed->isComplete())
      {
         /* from presolved problem */
         if(!seeed->isFromUnpresolved())
         {
            ++ndetectedpresolved;
         }
         /* from original problem */
         else
         {
            ++ndetectedunpresolved;
         }
      }
   }

   /* count width of menu table by summing width of column headers,
    * make header line, and border line for header (will be beneath the column headers) as follows:
    * a border line of the length of the column width as '-' for each column and a space between the columns,
    * e.g. header line "   nr   id nbloc nmacon  sel ",
    * e.g. underscores " ---- ---- ----- ------ ---- " */
   std::string headerline;
   std::string borderline;
   std::map<std::string, int> columnlength;
   int linelength = 0;

   /* line starts with a space */
   headerline = " ";
   borderline = " ";

   for(auto header : columns)
   {
      /* make sure the header name is unique and add a length for header */
      assert(columnlength.find(header) == columnlength.end());
      columnlength.insert(std::pair<std::string,int>(header, 0));
      /* if header is smaller than min column width, add spaces to header first */
      if(header.size() < DEFAULT_COLUMN_MIN_WIDTH)
      {
         for(int i = 0; i < (DEFAULT_COLUMN_MIN_WIDTH - (int) header.size()); i++)
         {
            headerline += " ";
            borderline += "-";
            columnlength.at(header)++;
         }
      }
      /* add header to headerline and add #chars of header as '-' to borderline*/
      headerline += header;
      for(int i = 0; i < (int) header.size(); i++)
         borderline += "-";
      columnlength.at(header) += (int) header.size();
      /* add space to both lines as column border */
      headerline += " ";
      borderline += " ";
      /* add columnlength (+1 for border space) to overall linelength */
      linelength += columnlength.at(header) + 1;
   }

   /* display overview statistics */
   SCIPdialogMessage(scip, NULL, "\n");
   outputCharXTimes(scip, '=', linelength);
   SCIPdialogMessage(scip, NULL, " \n");
   SCIPdialogMessage(scip, NULL, "Summary              presolved       original \n");
   SCIPdialogMessage(scip, NULL, "                     ---------       -------- \n");
   SCIPdialogMessage(scip, NULL, "detected             ");
   SCIPdialogMessage(scip, NULL, "%9d       ", ndetectedpresolved );
   SCIPdialogMessage(scip, NULL, "%8d\n", ndetectedunpresolved );
   outputCharXTimes(scip, '=', linelength);
   SCIPdialogMessage(scip, NULL, " \n");
   /* display header of table */
   SCIPdialogMessage(scip, NULL, "%s\n", headerline.c_str());
   SCIPdialogMessage(scip, NULL, "%s\n", borderline.c_str());

   /*@todo make this part generic */
   /* go through all seeeds that should currently be displayed,
    * so from startindex on menulength many entries if there are that much left in the id list */
   for(int i = startindex; i < startindex + menulength && i < *listlength; ++i)
   {
      Seeed_Wrapper sw;
      Seeed* seeed;
      SCIP_CALL( GCGgetSeeedFromID(scip, &(*idlist)[i], &sw) );
      seeed = sw.seeed;

      SCIPdialogMessage(scip, NULL, " %4d   ", i );
      SCIPdialogMessage(scip, NULL, "%4d  ", seeed->getID() );
      SCIPdialogMessage(scip, NULL, "%5d  ", seeed->getNBlocks() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMasterconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNLinkingvars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMastervars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNTotalStairlinkingvars() );
      if( seeed->isComplete() )
         SCIPdialogMessage(scip, NULL, "%.4f  ",  seeed->getScore(SCIPconshdlrDecompGetScoretype(scip)) );
      else
         SCIPdialogMessage(scip, NULL, "<=%.2f  ", seeed->getScore(SCIPconshdlrDecompGetScoretype(scip)) );
      SCIPdialogMessage(scip, NULL, "%7s  ", seeed->getDetectorChainString() );
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->isFromUnpresolved() ? "no" : "yes")  );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenvars() );
      SCIPdialogMessage(scip, NULL, "%3s  \n", (seeed->isSelected() ? "yes" : "no")  );
   }

   SCIPdialogMessage(scip, NULL, "=================================================================================================== \n");

   return SCIP_OKAY;
}


/** Shows information about the explore screen and its abbreviations
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogShowLegend(
   SCIP* scip,             /**< SCIP data structure */
   SCORETYPE scoretype     /**< current score type */
   )
{
   assert(scip != NULL);
   DEC_DETECTOR** detectors;

   std::string scorename = SCIPconshdlrDecompGetScoretypeShortName(scip, scoretype);
   std::string scoredescr = SCIPconshdlrDecompGetScoretypeDescription(scip, scoretype);

   SCIPdialogMessage(scip, NULL, "List of included detectors for decompositions histories: \n");

   SCIPdialogMessage(scip, NULL, "\n%30s    %4s\n", "detector" , "char");
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "--------" , "----");

   detectors = SCIPconshdlrDecompGetDetectors(scip);

   for( int det = 0; det < SCIPconshdlrDecompGetNDetectors(scip); ++det )
   {
      DEC_DETECTOR* detector;
      detector = detectors[det];

      SCIPdialogMessage(scip, NULL, "%30s    %4c\n", DECdetectorGetName(detector), DECdetectorGetChar(detector));
   }
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "given by user" , "U");
   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "=================================================================================================== \n");

   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "List of abbreviations of decomposition table \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "abbreviation", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "------------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nr", "number of the decomposition (use this number for choosing the decomposition)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "id", "id of the decomposition (identifies the decomposition in reports/statistics/visualizations/etc.)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nbloc", "number of blocks");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nmacon", "number of master constraints");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nlivar", "number of linking variables");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nmavar", "number of master variables (do not occur in blocks)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nstlva", "number of stairlinking variables (disjoint from linking variables)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", scorename.c_str(), scoredescr.c_str());
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "history", "list of detector chars worked on this decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "pre", "is this decomposition for the presolved problem");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopcon", "number of open constraints");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopvar", "number of open variables");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "sel", "is this decomposition selected at the moment");

   SCIPdialogMessage(scip, NULL, "\n=================================================================================================== \n");

   return SCIP_OKAY;
}

/** Shows help section of explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogShowHelp(
   SCIP* scip  /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "=================================================================================================== \n");
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "List of selection commands \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "command", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "-------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "select", "selects/unselects decomposition with given id");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "previous", "displays the preceding decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "next", "displays the subsequent decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "top", "displays the first decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "end", "displays the last decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "legend", "displays the legend for table header and history abbreviations");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "help", "displays this help");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "number_entries", "modifies the number of displayed decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "visualizes the specified decomposition (requires gnuplot)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "inspect", "displays detailed information for the specified decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "calc_strong", "calculates and displays the strong decomposition score for this decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "return to main menu");

   SCIPdialogMessage(scip, NULL, "\n=================================================================================================== \n");

   return SCIP_OKAY;
}


/** Shows a visualization of the seeed specified by the user via the dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogSelectVisualize(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   int**                   idlist,     /**< current list of seeed ids */
   int*                    listlength  /**< length of idlist */
   )
{
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;
   int commandlen;

   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the nr of the decomposition to be visualized:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0 )
      idtovisu = atoi(ntovisualize);

   /* check whether the seeed exists */
   if( commandlen == 0 || idtovisu < 0 || idtovisu >= *listlength )
   {
      SCIPdialogMessage( scip, NULL, "This nr is out of range." );
      return SCIP_OKAY;
   }

   /* get and show seeed */
   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &(*idlist)[idtovisu], &sw) );
   Seeed* seeed = sw.seeed;
   assert( seeed != NULL );

   seeed->showVisualisation();

   return SCIP_OKAY;
}


/**
 * Calculates and displays the strong decomposition score for this decomposition in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE SCIPdialogCalcStrongDecompositionScore(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   int**                   idlist,     /**< current list of seeed ids */
   int*                    listlength  /**< length of idlist */
   )
{
   char* ntocalcstrong;
   SCIP_Bool endoffile;
   int idtocalcstrong;
   int commandlen;

   assert( scip != NULL );

   /* read the id of the decomposition to calculate strong decomp score for */
   SCIPdialogMessage( scip, NULL,
      "Please specify the nr of the decomposition that should be evaluated by strong decomposition score:\n" );
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntocalcstrong, &endoffile ) );
   commandlen = strlen( ntocalcstrong );

   idtocalcstrong = -1;
   if( commandlen != 0 )
   {
      std::stringstream convert( ntocalcstrong );
      convert >> idtocalcstrong;

      if( idtocalcstrong == 0 && ntocalcstrong[0] != '0' )
      {
         idtocalcstrong = -1;
      }
   }

   /* check whether the seeed exists */
   if( commandlen == 0 || idtocalcstrong < 0 || idtocalcstrong >= *listlength )
   {
      SCIPdialogMessage( scip, NULL, "This nr is out of range." );
      return SCIP_OKAY;
   }

   /* call calculation strong decomp score method according to chosen parameters */
   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &(*idlist)[idtocalcstrong], &sw) );
   assert(sw.seeed != NULL);

   SCIP_Real score;
   Seeedpool* seeedpool = sw.seeed->getSeeedpool();
   seeedpool->calcStrongDecompositionScore(sw.seeed, &score);
   SCIPdialogMessage( scip, NULL, "Strong decomposition score of this decomposition is %f.", score) ;

   return SCIP_OKAY;
}


/**
 * Displays information about a seeed that is chosen by the user in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE SCIPdialogInspectSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   int**                   idlist,     /**< current list of seeed ids */
   int*                    listlength  /**< length of idlist */
   )
{
   char* ntoinspect;
   char* ndetaillevel;
   SCIP_Bool endoffile;
   int idtoinspect;
   int detaillevel;

   int commandlen;

   assert( scip != NULL );

   /* read the id of the decomposition to be inspected */
   SCIPdialogMessage( scip, NULL, "Please specify the nr of the decomposition to be inspected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntoinspect, &endoffile ) );
   commandlen = strlen( ntoinspect );

   idtoinspect = -1;
   if( commandlen != 0 )
      idtoinspect = atoi( ntoinspect );

   if(idtoinspect < 0 || idtoinspect >= *listlength){
      SCIPdialogMessage( scip, NULL, "This nr is out of range." );
      return SCIP_OKAY;
   }

   /* check whether ID is in valid range */
   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &(*idlist)[idtoinspect], &sw) );

   assert( sw.seeed != NULL );

   /* read the desired detail level; for wrong input, it is set to 1 by default */
   SCIPdialogMessage( scip, NULL,
      "Please specify the detail level:\n  0 - brief overview\n  1 - block and detector info (default)\n  2 - cons and var assignments\n" );
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ndetaillevel, &endoffile ) );
   commandlen = strlen( ndetaillevel );

   detaillevel = 1;
   if( commandlen != 0 )
   {
      std::stringstream convert( ndetaillevel );
      convert >> detaillevel;

      if( detaillevel < 0 || ( detaillevel == 0 && ndetaillevel[0] != '0' ) )
      {
         detaillevel = 1;
      }
   }

   /* call displayInfo method according to chosen parameters */
   sw.seeed->displayInfo( detaillevel );

   return SCIP_OKAY;
}


/** Lets the user select decompositions from the explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogSelect(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   int**                   idlist,     /**< current list of seeed ids */
   int*                    listlength  /**< length of idlist */
   )
{
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;

   int commandlen;

   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the nr of the decomposition to be selected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   if( commandlen == 0 || idtovisu < 0 || idtovisu >= *listlength )
   {
      SCIPdialogMessage( scip, NULL, "This nr is out of range, nothing was selected." );
      return SCIP_OKAY;
   }

   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &(*idlist)[idtovisu], &sw) );

   assert( sw.seeed != NULL );

   sw.seeed->setSelected(!sw.seeed->isSelected() );

   if( sw.seeed->isSelected() )
      std::cout << "is selected!" << sw.seeed->isSelected() << std::endl;

   return SCIP_OKAY;
}

static
SCIP_RETCODE SCIPdialogExecCommand(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog,
   std::vector<std::string> columns,
   char*                   command,
   SCIP_Bool               endoffile,
   int*                    startindex,
   int*                    menulength,
   SCIP_Bool*              finished,
   SCORETYPE*              scoretype,
   int*                    nseeeds,
   int**                   idlist,
   int*                    listlength
   )
{
   

   int commandlen = strlen(command);

      if( strncmp( command, "previous", commandlen) == 0 )
      {
         *startindex = *startindex - *menulength;
         if(*startindex < 0 )
            *startindex = 0;
      }
      else if( strncmp( command, "next", commandlen) == 0 )
      {
         *startindex = *startindex + *menulength;
         if( *startindex > *listlength - *menulength )
            *startindex = *listlength - *menulength;
      }
      else if( strncmp( command, "top", commandlen) == 0 )
      {
         *startindex = 0;
      }
      else if( strncmp( command, "end", commandlen) == 0 )
      {
         *startindex = *listlength - *menulength;
      }

      else if( strncmp( command, "quit", commandlen) == 0 )
      {
         *finished = TRUE;
         SCIP_CALL( SCIPconshdlrDecompChooseCandidatesFromSelected(scip) );
      }

      else if( strncmp( command, "legend", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogShowLegend(scip, *scoretype) );
      }

      else if( strncmp( command, "help", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogShowHelp(scip) );
      }

      else if( strncmp( command, "number_entries", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogSetNEntires(scip, dialoghdlr, dialog, menulength) );
      }

      else if( strncmp( command, "visualize", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogSelectVisualize(scip, dialoghdlr, dialog, idlist, listlength) );
      }

      else if( strncmp( command, "inspect", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogInspectSeeed( scip, dialoghdlr, dialog, idlist, listlength) );
      }

      else if( strncmp( command, "calc_strong", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogCalcStrongDecompositionScore(scip, dialoghdlr, dialog, idlist, listlength) );
      }

      else if( strncmp( command, "select", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogSelect(scip, dialoghdlr, dialog, idlist, listlength) );
      }

   return SCIP_OKAY;
}

extern "C" {

SCIP_RETCODE GCGdialogExecExplore(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   )
{
   /* set navigation defaults */
   int startindex = 0;
   int menulength = DEFAULT_MENULENGTH;

   /* check for available seeeds */
   int nseeeds = SCIPconshdlrDecompGetNSeeeds(scip);
   if(nseeeds == 0)
   {
      SCIPdialogMessage( scip, NULL, "There are no decompositions to explore yet, please detect first.\n" );
      return SCIP_OKAY;
   }

   /* get initial seeed list */
   SCIP_CALL( SCIPdialogUpdateSeeedlist(scip, &startindex) );
   int* idlist;
   int listlength;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idlist, nseeeds) );

   /* set initial columns */
   std::vector<std::string> columns;
   SCORETYPE scoretype = scoretype::MAX_WHITE;
   char* scorename;
   char columnstr[] = DEFAULT_COLUMNS;
   char* tempcolumns = strtok(columnstr, " ");
   while(tempcolumns != NULL)
   {
      /* get each column header of default */
      char newchar[DEFAULT_COLUMN_MAX_WIDTH]; // cutting string at max column width if longer
      strcpy(newchar, tempcolumns);

      /* "score" is a wildcard for the current score */
      if( strncmp(newchar, "score", strlen(newchar)) == 0 )
      {
         scoretype = SCIPconshdlrDecompGetScoretype(scip);
         scorename = SCIPconshdlrDecompGetScoretypeShortName(scip, scoretype);
         char name[DEFAULT_COLUMN_MAX_WIDTH];
         strcpy(name, scorename);
         columns.push_back(name);
      }
      else
      {
         /* if the name is not score, just use the current char */
         columns.push_back(newchar);
      }
      /* get the next item in the list */
      tempcolumns = strtok (NULL, " ");
   }
   /*@todo hand this vector, scorename, scoretype down and use it to make menu columns generic */

   /* while user has not aborted: show current list extract and catch commands */
   SCIP_Bool finished = false;
   char* command;
   SCIP_Bool endoffile;
   while( !finished )
   {
      SCIPdialogShowMenu(scip, columns, &nseeeds, startindex, menulength, &idlist, &listlength);

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Please enter command or decomposition id to select (or \"h\" for help) : \nGCG/explore> ", &command, &endoffile) );

      SCIPdialogExecCommand(scip, dialoghdlr, dialog, columns, command, endoffile, &startindex, &menulength, &finished, &scoretype, &nseeeds, &idlist, &listlength);
   }

   SCIPfreeBlockMemoryArray(scip, &idlist, nseeeds);
   return SCIP_OKAY;
}

} // extern "C"

} // namespace gcg
