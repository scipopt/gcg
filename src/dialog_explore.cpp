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

#include "class_seeed.h"
#include "cons_decomp.h"
#include "wrapper_seeed.h"

#define DEFAULT_MENULENGTH 10

namespace gcg
{

/**
 * Gets the shortname of the given scoretype
 *
 * @returns the shortname of the given Scoretype
 */
static
char* SCIPgetScoretypeShortName(
   SCIP*       scip,    /**< SCIP data structure */
   SCORETYPE   sctype   /**< scoretype */
   )
{
   char scoretypename[SCIP_MAXSTRLEN];
   char* copy;
   /* set detector chain info string */
   SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "") ;


   if( sctype == scoretype::MAX_WHITE )
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maxwhi") ;

   if( sctype == scoretype::CLASSIC )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "classi") ;

   if( sctype == scoretype::BORDER_AREA )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "border") ;

   if( sctype == scoretype::MAX_FORESSEEING_WHITE )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "forswh") ;

   if( sctype == scoretype::MAX_FORESEEING_AGG_WHITE )
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "fawh") ;


   if( sctype == scoretype::SETPART_FWHITE )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "spfwh ") ;

   if( sctype == scoretype::SETPART_AGG_FWHITE )
           SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "spfawh") ;


   if( sctype == scoretype::BENDERS )
              SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "bender") ;


   SCIP_CALL_ABORT ( SCIPduplicateBlockMemoryArray(scip, &copy, scoretypename, SCIP_MAXSTRLEN ) );

   return copy;

}

/*!
 * returns the description of the given scoretype
 *
 * @returns description of the scoretype
 */
static
char*  SCIPconshdlrDecompGetScoretypeDescription(
   SCIP*       scip,    /**< SCIP data structure */
   SCORETYPE   sctype   /**< scoretype */
      )
{
   char scoretypename[SCIP_MAXSTRLEN];
   char* copy;

   /* set detector chain info string */
   SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "") ;

   if( sctype == scoretype::MAX_WHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum white area score (i.e. maximize fraction of white area score; white area is nonblock and nonborder area, stairlinking variables count as linking)") ;

   if( sctype == scoretype::CLASSIC)
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "classical score") ;

   if( sctype == scoretype::BORDER_AREA)
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "minimum border score (i.e. minimizes fraction of border area score; )")  ;

   if( sctype == scoretype::MAX_FORESSEEING_WHITE)
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum foreseeing  white area score (i.e. maximize fraction of white area score considering problem with copied linking variables and corresponding master constraints; white area is nonblock and nonborder area, stairlinking variables count as linking)")  ;

   if( sctype == scoretype::MAX_FORESEEING_AGG_WHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum foreseeing  white area score with aggregation information(i.e. maximize fraction of white area score considering problem with copied linking variables and corresponding master constraints; white area is nonblock and nonborder area, stairlinking variables count as linking)")  ;


   if( sctype == scoretype::SETPART_FWHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "setpartitioning maximum foreseeing  white area score (i.e. convex combination of maximum foreseeing white area score and a boolean score rewarding a master containing only setppc and cardinality constraints )")  ;

   if( sctype == scoretype::SETPART_AGG_FWHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "setpartitioning maximum foreseeing white area score with aggregation information (i.e. convex combination of maximum foreseeing white area score and a boolean score rewarding a master containing only setppc and cardinality constraints )")  ;

   if( sctype == scoretype::BENDERS)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "experimental score to evaluate benders decompositions")  ;


   SCIP_CALL_ABORT ( SCIPduplicateBlockMemoryArray(scip, &copy, scoretypename, SCIP_MAXSTRLEN ) );

   return copy;

}

/** Shows header for seeed information in explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogShowListExtractHeader(
   SCIP*                   scip  /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   int ndetectedpresolved;
   int ndetectedunpresolved;
   int nuserpresolvedfull;
   int nuserpresolvedpartial;
   int nuserunpresolvedfull;
   int nuserunpresolvedpartial;
   char* scorename;
   size_t i;

   scorename = SCIPgetScoretypeShortName(scip, SCIPconshdlrDecompGetScoretype(scip) );

   ndetectedpresolved = 0;
   ndetectedunpresolved = 0;
   nuserpresolvedfull = 0;
   nuserpresolvedpartial = 0;
   nuserunpresolvedfull = 0;
   nuserunpresolvedpartial = 0;
   int* idlist; //@todo SCIP malloc for nseeeds
   int listlength;
   SCIPconshdlrDecompGetSelectedSeeeds(scip, &idlist, &listlength);

   /* count corresponding seeeds */
   for ( i = 0; i < listlength; ++i )
   {
      Seeed_Wrapper sw;
      GCGgetSeeedFromID(scip, &(idlist[i]), &sw);
      Seeed* seeed;
      seeed = sw.seeed;
      if( seeed->isComplete() && seeed->getUsergiven() == USERGIVEN::NOT && !seeed->isFromUnpresolved() )
         ++ndetectedpresolved;
      if( seeed->isComplete() && seeed->getUsergiven() == USERGIVEN::NOT && seeed->isFromUnpresolved() )
         ++ndetectedunpresolved;
      if( seeed->isComplete() && ( seeed->getUsergiven() == USERGIVEN::COMPLETE || seeed->getUsergiven() == USERGIVEN::COMPLETED_CONSTOMASTER) && !seeed->isFromUnpresolved() )
         ++nuserpresolvedfull;
      if( !seeed->isComplete() && seeed->getUsergiven() == USERGIVEN::PARTIAL && !seeed->isFromUnpresolved() )
         ++nuserpresolvedpartial;
      if( seeed->isComplete() && ( seeed->getUsergiven() == USERGIVEN::COMPLETE || seeed->getUsergiven() == USERGIVEN::COMPLETED_CONSTOMASTER) && seeed->isFromUnpresolved() )
         ++nuserunpresolvedfull;
      if( !seeed->isComplete() && seeed->getUsergiven() == USERGIVEN::PARTIAL && seeed->isFromUnpresolved() )
         ++nuserunpresolvedpartial;

   }

   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, "============================================================================================= ");
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, "Summary              presolved       original \n");
   SCIPdialogMessage(scip, NULL, "                     ---------       -------- \n");
   SCIPdialogMessage(scip, NULL, "detected             ");
   SCIPdialogMessage(scip, NULL, "%9d       ", ndetectedpresolved );
   SCIPdialogMessage(scip, NULL, "%8d\n", ndetectedunpresolved );
   SCIPdialogMessage(scip, NULL, "user given (partial) ");
   SCIPdialogMessage(scip, NULL, "%9d       ", nuserpresolvedpartial );
   SCIPdialogMessage(scip, NULL, "%8d\n", nuserunpresolvedpartial );
   SCIPdialogMessage(scip, NULL, "user given (full)    ");
   SCIPdialogMessage(scip, NULL, "%9d       ", nuserpresolvedfull );
   SCIPdialogMessage(scip, NULL, "%8d\n", nuserunpresolvedfull );

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");
   SCIPdialogMessage(scip, NULL, "   id   nbloc  nmacon  nlivar  nmavar  nstlva  %.6s  history  pre  nopcon  nopvar  usr  sel \n", scorename );
   SCIPdialogMessage(scip, NULL, " ----   -----  ------  ------  ------  ------  ------  -------  ---  ------  ------  ---  --- \n");

   SCIPfreeBlockMemoryArrayNull(scip, &scorename, SCIP_MAXSTRLEN);


   return SCIP_OKAY;
}


/** Shows detailed information about seeeds in explore menu
 *
 *@returns SCIP status
 * */
static
SCIP_RETCODE SCIPdialogShowListExtract(
   SCIP* scip,             /**< SCIP data structure */
   const int startindex    /**< index (in seeed list) of uppermost seeed in extract */
   )
{
   assert(scip != NULL);
   int i;

   int* idlist; //@todo malloc for nseeeds
   int listlength;
   SCIPconshdlrDecompGetSeeedLeafList(scip, &idlist, &listlength);

   for( i = startindex; i < (size_t) startindex + DEFAULT_MENULENGTH && i < (int) listlength; ++i)
   {
      Seeed_Wrapper sw;
      Seeed* seeed;
      SCIP_CALL( GCGgetSeeedFromID(scip, &(idlist[i]), &sw) );
      seeed = sw.seeed;

      assert( seeed->checkConsistency( ) );

      SCIPdialogMessage(scip, NULL, " %4d   ", i );
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
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->getUsergiven() == USERGIVEN::NOT ? "no" : "yes")   );
      SCIPdialogMessage(scip, NULL, "%3s  \n", (seeed->isSelected() ? "yes" : "no")  );
   }

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");

   return SCIP_OKAY;
}

/** Shows help for the user toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogShowToolboxInfo(
   SCIP* scip  /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Options to proceed: \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "option", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "conss", "assign unassigned constraints to master/blocks");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "vars", "assign unassigned variables to master(only)/linking/blocks");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "refine", "refine implicit constraint and variables assignments");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "finish", "choose a finishing detector that completes the decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "quit the modification process and returns to main menu");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "undo", "last modification is undone (atm only the last modification can be undone)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "shows a visualization of the current decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "propagate", "list all detectors that can propagate the current seeed and apply one to propagate it");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "finish", "list all detectors that can finish the current seeed and apply one to finish it");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "postprocess", "apply postprocessing to a finished seeed by selecting a suitable postprocessor");
   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");

   return SCIP_OKAY;
}

/** Shows information about the explore screen and its abbreviations
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogShowLegend(
   SCIP* scip  /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   DEC_DETECTOR** detectors;
   char * scorename;
   char * scoredescr;

   scorename = SCIPgetScoretypeShortName(scip, SCIPconshdlrDecompGetScoretype(scip) );
   scoredescr = SCIPconshdlrDecompGetScoretypeDescription(scip, SCIPconshdlrDecompGetScoretype(scip) );


   SCIPdialogMessage(scip, NULL, "List of included detectors for decompositions histories: \n" );

   SCIPdialogMessage(scip, NULL, "\n%30s    %4s\n", "detector" , "char"  );
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "--------" , "----"  );

   detectors = SCIPconshdlrDecompGetDetectors(scip);

   for( int det = 0; det < SCIPconshdlrDecompGetNDetectors(scip); ++det )
   {
      DEC_DETECTOR* detector;

      detector = detectors[det];

      SCIPdialogMessage(scip, NULL, "%30s    %4c\n", DECdetectorGetName(detector), DECdetectorGetChar(detector)  );
   }
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "given by user" , "U"  );

   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");

   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "List of abbreviations of decomposition table \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "abbreviation", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "------------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "id", "id of the decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nbloc", "number of blocks");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nmacon", "number of master constraints");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nlivar", "number of linking variables");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nmavar", "number of master variables (do not occur in blocks)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nstlva", "number of stairlinking variables (disjoint from linking variables)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", scorename, scoredescr);
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "history", "list of detector chars worked on this decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "pre", "is this decomposition for the presolved problem");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopcon", "number of open constraints");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopvar", "number of open variables");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "usr", "was this decomposition given by the user");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "sel", "is this decomposition selected at the moment");

   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");

   SCIPfreeBlockMemoryArrayNull(scip, &scorename, SCIP_MAXSTRLEN);
   SCIPfreeBlockMemoryArrayNull(scip, &scoredescr, SCIP_MAXSTRLEN);

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

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "List of selection commands \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "command", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "-------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "select", "selects/unselects decomposition with given id");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "modify", "modify an existing decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "create", "create a new decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "back", "displays the preceding decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "next", "displays the subsequent decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "top", "displays the first decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "end", "displays the last decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "legend", "displays the legend for table header and history abbreviations");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "help", "displays this help");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "dispNEntries", "modifies the number of displayed decompositions ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "finishes decomposition explorer and goes back to main menu");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "experimental feature: visualizes the specified decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "inspect", "displays detailed information for the specified decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "calc_strong", "calculates and displays the strong decomposition score for this decomposition");

   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");

   return SCIP_OKAY;
}


/** Shows a visualization of the seeed specified by the user via the dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogSelectVisualize(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;

   int commandlen;

   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be visualized:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0 )
      idtovisu = atoi(ntovisualize);

   /* check whether ID is in valid range */
   if( SCIPconshdlrDecompGetNSeeeds(scip) == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &idtovisu, &sw) );
   Seeed* seeed = sw.seeed;
   if( commandlen == 0 || seeed == NULL )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_OKAY;
   }

   seeed->showVisualisation();

   return SCIP_OKAY;
}


/**
 * Calculates and displays the strong decomposition score for this decomposition in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE SCIPdialogSelectCalcStrongDecompositionScore(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   char* ntocalcstrong;
   SCIP_Bool endoffile;
   int idtocalcstrong;
   int commandlen;

   assert( scip != NULL );

   /* read the id of the decomposition to be calculate strong decomp score */
   SCIPdialogMessage( scip, NULL,
      "Please specify the id of the decomposition that should be evaluated by strong decomposition score:\n" );
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntocalcstrong, &endoffile ) );
   commandlen = strlen( ntocalcstrong );

   idtocalcstrong = -1;
   if( commandlen != 0 )
   {
      std::stringstream convert( ntocalcstrong );
      convert >> idtocalcstrong;

      if ( idtocalcstrong == 0 && ntocalcstrong[0] != '0' )
      {
         idtocalcstrong = -1;
      }
   }

   /* call calculation strong decomp score method according to chosen parameters */
   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &idtocalcstrong, &sw) );
   if(sw.seeed != NULL)
   {
      SCIP_Real score;
      Seeedpool* seeedpool = sw.seeed->getSeeedpool();
      seeedpool->calcStrongDecompositionScore(sw.seeed, &score);
      SCIPdialogMessage( scip, NULL, "Strong decomposition score of this decomposition is %f.", score) ;
   }
   else
   {
      SCIPdialogMessage( scip, NULL, "This is not an existing id." );
   }

   return SCIP_OKAY;
}


/**
 * Displays information about a seeed that is chosen by the user in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE SCIPdialogSelectInspect(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
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
   SCIPdialogMessage( scip, NULL, "Please specify the id of the decomposition to be inspected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntoinspect, &endoffile ) );
   commandlen = strlen( ntoinspect );

   idtoinspect = -1;
   if( commandlen != 0 )
      idtoinspect = atoi( ntoinspect );

   /* check whether ID is in valid range */
   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &idtoinspect, &sw) );

   if( sw.seeed == NULL )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_PARAMETERWRONGVAL;
   }

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

      if ( detaillevel < 0 || ( detaillevel == 0 && ndetaillevel[0] != '0' ) )
      {
         detaillevel = 1;
      }
   }

   /* call displayInfo method according to chosen parameters */
   sw.seeed->displayInfo( detaillevel );

   return SCIP_OKAY;
}


/**
 * @brief sets a variable by name to the linking variables in the current user seeed
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPtoolboxUserSeeedSetVarToLinking(
   SCIP*                 scip,                /**< SCIP data structure */
   Seeed*                userseeed,           /**< user seeed to edit */
   const char*           varname              /**< name of the variable */
   )
{
   Seeedpool* currseeedpool;
   int varindex;

   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = userseeed->isFromUnpresolved()
      ? SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip)
      : SCIPconshdlrDecompGetSeeedpool(scip);

   varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname ) );

   userseeed->bookAsLinkingVar(varindex);

   return SCIP_OKAY;
}

/**
 * @brief sets the number of blocks
 *
 * set the number of blocks in the current user seeed
 * @returns SCIP return code
 */
extern
SCIP_RETCODE SCIPtoolboxUserSeeedSetnumberOfBlocks(
   SCIP*                 scip,                /**< SCIP data structure */
   Seeed*                userseeed,           /**< user seeed to edit */
   int                   nblocks              /**< number of blocks */
   )
{
   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   userseeed->setNBlocks(nblocks);

   return SCIP_OKAY;
}


/**
 * @brief sets a variable by name to the master in the current user seeed
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPtoolboxUserSeeedSetVarToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   Seeed*                userseeed,           /**< user seeed to edit */
   const char*           varname              /**< name of the variable */
   )
{
   Seeedpool* currseeedpool;
   int varindex;

   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = userseeed->isFromUnpresolved()
      ? SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip)
      : SCIPconshdlrDecompGetSeeedpool(scip);

   varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname) );

   userseeed->bookAsMasterVar(varindex);

   return SCIP_OKAY;
}


/**
 * @brief sets a constraint by name to a block in the current user seeed
 * @param scip SCIP data structure
 * @param consname name of the constraint that should be set to a block
 * @param blockid index of the block the constraint should be assigned to
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPtoolboxUserSeeedSetConsToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   Seeed*                userseeed,           /**< user seeed to edit */
   const char*           consname,            /**< name of the constraint */
   int                   blockid              /**< block index ( counting from 0) */
   )
{
   SCIP_CONS* cons;
   Seeedpool* currseeedpool;
   int consindex;

   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = userseeed->isFromUnpresolved()
         ? SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip)
         : SCIPconshdlrDecompGetSeeedpool(scip);

   cons = userseeed->isFromUnpresolved() ? (SCIPfindOrigCons(scip, consname) == NULL
      ? SCIPfindCons(scip, consname): SCIPfindOrigCons(scip, consname)) : SCIPfindCons(scip, consname);
   consindex = currseeedpool->getIndexForCons( cons ) ;

   if( blockid >= userseeed->getNBlocks() )
         userseeed->setNBlocks(blockid+1);
   userseeed->bookAsBlockCons(consindex, blockid);

   return SCIP_OKAY;
}


/**
 * @brief sets a variable by name to a block in the current user seeed
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPtoolboxUserSeeedSetVarToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   Seeed*                userseeed,           /**< user seeed to edit */
   const char*           varname,             /**< name of the variable */
   int                   blockid              /**< block index ( counting from 0) */
   )
{
   Seeedpool* currseeedpool;
   int varindex;

   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = userseeed->isFromUnpresolved()
         ? SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip)
         : SCIPconshdlrDecompGetSeeedpool(scip);

   varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname) );

   // if the block id is higher than expected, set the block to master
   if( blockid >= userseeed->getNBlocks() )
      userseeed->setNBlocks(blockid+1);
   userseeed->bookAsBlockVar(varindex, blockid);

   return SCIP_OKAY;
}


/**
 * @brief sets a constraint by name to master in the current user seeed
 * @param consname of the constraint that should be set to master
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPtoolboxUserSeeedSetConsToMaster(
   SCIP*                 scip,      /**< SCIP data structure */
   Seeed*                userseeed, /**< user seeed to edit */
   const char*           consname   /**< name of cons to book as master cons */
   )
{
   Seeedpool* currseeedpool;
   int consindex;
   SCIP_CONS* cons;

   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = userseeed->isFromUnpresolved()
      ? SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip)
      : SCIPconshdlrDecompGetSeeedpool(scip);

   cons = userseeed->isFromUnpresolved() ? SCIPfindOrigCons(scip, consname) : SCIPfindCons(scip, consname);
   consindex = currseeedpool->getIndexForCons( cons );

   userseeed->bookAsMasterCons(consindex);
   return SCIP_OKAY;
}


/**
 * @brief rejects and deletes the current user seeed
 * @returns SCIP return code
 */
static
SCIP_RETCODE GCGtoolboxUserSeeedReject(
  SCIP* scip, /**< SCIP data structure */
  Seeed* userseeed
  )
{
   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one  before you can reject it\n");
      return SCIP_OKAY;
   }
   delete userseeed;

   return SCIP_OKAY;
}


/**
 * finalizes and flushes the current user seeed, i.e. consider implicits, calc hashvalue, construct decdecomp if
 * complete etc
 * @returns SCIP return code
 */
static
SCIP_RETCODE GCGtoolboxUserSeeedFlush(
   SCIP* scip,       /**< SCIP data structure */
   Seeed* userseeed  /**< user seeed to be flushed */
   )
{
   char const* usergiveninfo;
   char const* presolvedinfo;

   if( userseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   Seeedpool* currseeedpool = userseeed->isFromUnpresolved()
      ? SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip)
      : SCIPconshdlrDecompGetSeeedpool(scip);

   userseeed->setSeeedpool(currseeedpool);
   userseeed->flushBooked();

   if( userseeed->shouldCompletedByConsToMaster() )
   {
      for(int opencons = 0; opencons < userseeed->getNOpenconss(); ++opencons)
         userseeed->bookAsMasterCons( userseeed->getOpenconss()[opencons] );
      userseeed->flushBooked();
   }

   userseeed->considerImplicits();
   currseeedpool->prepareSeeed(userseeed);

   if( !userseeed->checkConsistency() )
   {
      GCGtoolboxUserSeeedReject(scip, userseeed);
      SCIPwarningMessage(scip, "Your seeed was rejected because of inconsistencies! \n");
      return SCIP_OKAY;
   }
   userseeed->buildDecChainString();
   if( userseeed->isComplete() )
   {
      if( !userseeed->shouldCompletedByConsToMaster() )
         userseeed->setUsergiven( USERGIVEN::COMPLETE );

      if( !userseeed->isFromUnpresolved() )
      {
         SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, userseeed));
      }
      /* stems from unpresolved problem */
      else
      {

         SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForUnpresolved(scip, userseeed) );

         /* if seeedpool for presolved problem already exist try to translate seeed */
         if ( SCIPconshdlrDecompGetSeeedpool(scip) != NULL )          {
            std::vector<Seeed*> seeedtotranslate(0);
            std::vector<Seeed*> newseeeds(0);
            seeedtotranslate.push_back(userseeed);
            SCIPconshdlrDecompGetSeeedpool(scip)->translateSeeeds(SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip), seeedtotranslate, newseeeds);
            if( newseeeds.size() != 0 )
            {
               SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, newseeeds[0]) );
            }
         }
      }
   }
   else
   {
      assert( !userseeed->shouldCompletedByConsToMaster() );
      userseeed->setUsergiven( USERGIVEN::PARTIAL );

      if ( !userseeed->isFromUnpresolved() )
         SCIP_CALL(SCIPconshdlrDecompAddPartialSeeedForPresolved(scip, userseeed) );
      else
         SCIP_CALL(SCIPconshdlrDecompAddPartialSeeedForUnpresolved(scip, userseeed) );
   }

   /* set statistics */
   {
      int nvarstoblock = 0;
      int nconsstoblock = 0;

      for ( int b = 0; b < userseeed->getNBlocks(); ++b )
      {
         nvarstoblock += userseeed->getNVarsForBlock(b);
         nconsstoblock += userseeed->getNConssForBlock(b);
      }
      userseeed->setDetectorPropagated(NULL);

      userseeed->addClockTime(0.);
      userseeed->addPctVarsFromFree( (nvarstoblock + userseeed->getNMastervars() +userseeed->getNLinkingvars())/(SCIP_Real) userseeed->getNVars()  );
      userseeed->addPctVarsToBlock((nvarstoblock )/(SCIP_Real) userseeed->getNVars() );
      userseeed->addPctVarsToBorder( (userseeed->getNMastervars() +userseeed->getNLinkingvars())/(SCIP_Real) userseeed->getNVars() ) ;
      userseeed->addPctConssToBorder( (userseeed->getNMasterconss() ) / (SCIP_Real) userseeed->getNConss() ) ;
      userseeed->addPctConssFromFree( (userseeed->getNMasterconss() + nconsstoblock ) / (SCIP_Real) userseeed->getNConss() ) ;
      userseeed->addPctConssToBlock( (nconsstoblock ) / (SCIP_Real) userseeed->getNConss() );
      userseeed->addNNewBlocks(userseeed->getNBlocks());
   }

   userseeed->findVarsLinkingToMaster();
   userseeed->findVarsLinkingToStairlinking();


   if( userseeed->getUsergiven() == USERGIVEN::PARTIAL )
      usergiveninfo = "partial";
   if( userseeed->getUsergiven() == USERGIVEN::COMPLETE )
      usergiveninfo = "complete";
   if( userseeed->getUsergiven() == USERGIVEN::COMPLETED_CONSTOMASTER )
         usergiveninfo = "complete";
   if( userseeed->isFromUnpresolved() )
         presolvedinfo = "unpresolved";
   else presolvedinfo = "presolved";

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " added %s decomp for %s problem with %d blocks and %d masterconss, %d linkingvars, "
      "%d mastervars, and max white score of %s %f \n", usergiveninfo, presolvedinfo,
      userseeed->getNBlocks(), userseeed->getNMasterconss(),
      userseeed->getNLinkingvars(), userseeed->getNMastervars(), (userseeed->isComplete() ? " " : " at best "),
      userseeed->getScore(SCORETYPE::MAX_WHITE) );

   return SCIP_OKAY;
}


/**
 * @brief creates a user seeed for the problem
 *
 * creates a new seeed as specified by the user
 * @note this function assumes that the corresponding seeedpool exists
 * @returns new user seeed
 */
static
Seeed* GCGtoolboxUserSeeedCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             presolved,          /**< should the user seeed be created for the presolved problem */
   SCIP_Bool             markedincomplete    /**< should the user seeed be a partial one */
   )
{
   Seeed_Wrapper swpool = presolved ? SCIPconshdlrDecompGetSeeedpool(scip)
      : SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip);
   Seeedpool* currseeedpool = swpool.seeedpool;

   assert( currseeedpool != NULL );

   Seeed* userseeed = new Seeed(scip, currseeedpool->getNewIdForSeeed(), currseeedpool);
   userseeed->setIsFromUnpresolved( !presolved );

   if( markedincomplete )
      userseeed->setUsergiven(USERGIVEN::PARTIAL);
   else
      userseeed->setUsergiven(USERGIVEN::COMPLETED_CONSTOMASTER);

   return SCIP_OKAY;
}


/** Lets user modify conss during modification of seeed in toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxModifyConss(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   Seeed*                  userseeed   /**< seeed to be modified */
   )
{
    SCIP_Bool matching;
    char* consregex;
    char* command;
    char* command2;
    SCIP_Bool endoffile;
    int commandlen;

    assert(scip != NULL);
    matching = FALSE;

    Seeedpool* seeedpool;
    std::vector<int> matchingconss  = std::vector<int>(0);

    seeedpool = userseeed->getSeeedpool();
    /* Does user want to modify existing or create a new partial decomposition ?*/
    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
       "Please specify a regular expression (modified ECMAScript regular expression grammar) matching the names of unassigned constraints you want to assign : \nGCG/toolbox> ",
       &consregex, &endoffile) );

    /* case distinction: */

    std::regex expr;
    try  {
       expr = std::regex(consregex);
    }
    catch (const std::regex_error& e) {
       std::cout << "regex_error caught: " << e.what() << '\n';
       if (e.code() == std::regex_constants::error_brack) {
          std::cout << "The code was error_brack\n";
       }
    }

    for( int oc = 0; oc < userseeed->getNOpenconss(); ++oc )
    {
       const char* consname;

       consname = SCIPconsGetName( seeedpool->getConsForIndex(userseeed->getOpenconss()[oc] ) );

       if( std::regex_match(consname, expr) )
       {
          matching = TRUE;
          matchingconss.push_back(userseeed->getOpenconss()[oc]);
          SCIPdebugMessage(" consname %s matches regex %s \n", consname, consregex );
       } else
          SCIPdebugMessage(" consname %s does not match regex %s \n", consname, consregex);
    }

    if( !matching )
    {
       SCIPdialogMessage(scip, NULL, " There are no unassigned constraints with names matching given regular expression. Return to toolbox main menu.\n");
       return SCIP_OKAY;
    }

    //@todo note: here was the block updating last user seeed

    if( matchingconss.size() > 10 )
       SCIPdebugMessage(" There are %d unassigned constraints with names matching given regular expression. Showing the first 10:\n", (int) matchingconss.size());
    else
       SCIPdebugMessage(" There are %d unassigned constraints with names matching given regular expression: \n", (int) matchingconss.size());

    for( size_t mc = 0 ; mc < 10 && mc < matchingconss.size(); ++mc )
       SCIPdialogMessage(scip, NULL, " %s \n", SCIPconsGetName( seeedpool->getConsForIndex( matchingconss[mc] ) ));

    SCIPdialogMessage(scip, NULL, "\n Should these constraints be added to: \n");
    SCIPdialogMessage(scip, NULL, " master \n");
    SCIPdialogMessage(scip, NULL, " block (to be specified) \n");
    SCIPdialogMessage(scip, NULL, " nothing (return to toolbox main menu)? \n");


    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please specify how to proceed: \nGCG/toolbox> ", &command, &endoffile) );

    commandlen = strlen(command);

    /* case distinction: */
    if( strncmp( command, "master", commandlen) == 0 )
    {
       for( size_t mc = 0 ;  mc < matchingconss.size(); ++mc )
       {
          userseeed->bookAsMasterCons( matchingconss[mc] );
       }
    }
    else if( strncmp( command, "block", commandlen) == 0 )
    {
       SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please specify the block number these constraints should be assigned to: \nGCG/toolbox> ", &command2, &endoffile) );
       char* tail;
       int blockid = strtol(command2, &tail, 10);
       for( size_t mc = 0 ;  mc < matchingconss.size(); ++mc )
       {
          userseeed->bookAsBlockCons( matchingconss[mc], blockid );
       }
    }
    else
       return SCIP_OKAY;

    userseeed->flushBooked();

   return SCIP_OKAY;
}

/** Lets user specify how to finish the modified seeed while using the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPdialogToolboxModifyFinish(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   Seeed*                  userseeed
   )
{
   SCIP_Bool         choosenfinisher;

   char* command;
   SCIP_Bool endoffile;
   char* tail;
   int finisherid;
   SEEED_PROPAGATION_DATA* seeedPropData;
   DEC_DETECTOR* finisher;
   SCIP_Result result;

   assert(scip != NULL);

   Seeedpool* seeedpool;
   std::vector<int> matchingvars  = std::vector<int>(0);

   seeedpool = userseeed->getSeeedpool();
   choosenfinisher = FALSE;
   while ( !choosenfinisher )
   {
       SCIPdialogMessage(scip, NULL, " Available finisher: \n");
       /* 1) print out available finisher */
       SCIPdialogMessage(scip, NULL, "%d :  %s \n", -1, "abort" );
       for( int fi = 0; fi < seeedpool->getNFinishingDetectors(); ++fi )
       {
          SCIPdialogMessage(scip, NULL, "%d :  %s \n", fi, DECdetectorGetName(seeedpool->getFinishingDetectorForIndex(fi) ) );
       }

       /* Does user want to modify existing or create a new partial decomposition ?*/
       SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
          "Please specify the index of the finisher to use : \nGCG/toolbox> ", &command, &endoffile) );

       finisherid = strtol(command, &tail, 10);

       if( finisherid >= seeedpool->getNFinishingDetectors() || finisherid < -1 )
       {
            SCIPdialogMessage(scip, NULL, "The specified id is invalid \n"  );
            continue;
       }
       choosenfinisher = TRUE;
   }

   seeedPropData = new SEEED_PROPAGATION_DATA();
   seeedPropData->seeedpool = seeedpool;
   seeedPropData->nNewSeeeds = 0;
   seeedPropData->seeedToPropagate = userseeed;

   //@todo note: here was update of last user seeed

   finisher = seeedpool->getFinishingDetectorForIndex(finisherid);
   finisher->finishSeeed(scip, finisher, seeedPropData, &result);

   delete userseeed;

   for( int i = 0; i <  seeedPropData->nNewSeeeds; ++i)
   {
      delete seeedPropData->newSeeeds[i];
   }

   delete seeedPropData->seeedToPropagate;
   delete seeedPropData;

   return SCIP_OKAY;
}

/** Lets the user select a seeed to modify in toolbox
 *
 * @returns pointer to a copy of the chosen seeed */
static
Seeed* SCIPdialogToolboxChoose(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   char* ntochoose;
   SCIP_Bool endoffile;
   int idtochoose;

   Seeed* userseeed;

   int commandlen;

   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the (partial) decomposition to be chosen for modification:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntochoose, &endoffile) );
   commandlen = strlen(ntochoose);

   idtochoose = DEFAULT_MENULENGTH;
   if( commandlen != 0)
      idtochoose = atoi(ntochoose);

   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &idtochoose, &sw) );

   if( sw.seeed == NULL )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_PARAMETERWRONGVAL;
   }

   return new Seeed( sw.seeed );

   return SCIP_OKAY;
}

/** Lets user modify vars during use of the toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxModifyVars(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   Seeed*                  userseeed   /**< user seeed to edit */
   )
{
    SCIP_Bool matching;
    char* varregex;
    char* command;
    char* command2;
    SCIP_Bool endoffile;
    int commandlen;

    assert(scip != NULL);
    matching = FALSE;

    Seeedpool* seeedpool;
    std::vector<int> matchingvars  = std::vector<int>(0);

    seeedpool = userseeed->getSeeedpool();
    /* Does user want to modify existing or create a new partial decomposition ?*/
    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
       "Please specify a regular expression (modified ECMAScript regular expression grammar) matching the names of unassigned variables you want to assign : \nGCG/toolbox> ",
       &varregex, &endoffile) );

    std::regex expr;
    try  {
       expr = std::regex(varregex);
    }
    catch (const std::regex_error& e) {
       std::cout << "regex_error caught: " << e.what() << '\n';
       if (e.code() == std::regex_constants::error_brack) {
          SCIPdebugMessage("The code was error_brack\n");
       }
    }

    for( int oc = 0; oc < userseeed->getNOpenvars(); ++oc )
    {
       const char* varname;

       varname = SCIPvarGetName(seeedpool->getVarForIndex(userseeed->getOpenvars()[oc] ) );

       SCIPdebugMessage("check var %s for regex %s \n", varname, varregex);

       if( std::regex_match(varname, expr) )
       {
          matching = TRUE;
          matchingvars.push_back(userseeed->getOpenvars()[oc]);
          SCIPdebugMessage( " varname %s matches regex %s \n", varname, varregex );
       } else
          SCIPdebugMessage(" varname %s does not match regex %s \n", varname, varregex);
    }

    if( !matching )
    {
       SCIPdialogMessage(scip, NULL,
          " There are no unassigned variables with names matching given regular expression. Return to toolbox main menu.\n");
       return SCIP_OKAY;
    }

    //@todo here was the update of last user seeed

    if( matchingvars.size() > 10 )
       SCIPdialogMessage(scip, NULL,
          " There are %d unassigned variables with names matching given regular expression. Showing the first 10:\n",
          matchingvars.size());
    else
       SCIPdialogMessage(scip, NULL, " There are %d unassigned variables with names matching given regular expression: \n",
          matchingvars.size());

    for( size_t mc = 0 ; mc < 10 && mc < matchingvars.size(); ++mc )
       SCIPdialogMessage(scip, NULL, " %s \n", SCIPvarGetName( seeedpool->getVarForIndex( matchingvars[mc] ) ));

    SCIPdialogMessage(scip, NULL, "\n Should these variables be added to: \n");
    SCIPdialogMessage(scip, NULL, " master-only (static) \n");
    SCIPdialogMessage(scip, NULL, " linking \n");
    SCIPdialogMessage(scip, NULL, " block (to be specified) \n");
    SCIPdialogMessage(scip, NULL, " nothing (return to toolbox main menu)? \n");

    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please specify how to proceed: \nGCG/toolbox> ", &command,
       &endoffile) );

    commandlen = strlen(command);

    if( strncmp( command, "master", commandlen) == 0 )
    {
       for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
       {
          userseeed->bookAsMasterVar( matchingvars[mc] );
       }
    } else
       if( strncmp( command, "linking", commandlen) == 0 )
           {
              for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
              {
                 userseeed->bookAsLinkingVar( matchingvars[mc] );
              }
           }
    else if( strncmp( command, "block", commandlen) == 0 )
    {
       SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
          "Please specify the block number these variables should be assigned to: \nGCG/toolbox> ", &command2, &endoffile) );
       char* tail;
       int blockid = strtol(command2, &tail, 10);
       for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
       {
          userseeed->bookAsBlockVar( matchingvars[mc], blockid );
       }
    }
    else
       return SCIP_OKAY;

    userseeed->flushBooked();
    userseeed->deleteEmptyBlocks(true);

   return SCIP_OKAY;
}

/** Apply propagation, finishing or postprocessing to the current user seeed via dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxActOnSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   Seeed*                  userseeed,  /**< seeed to edit */
   toolboxtype             action      /**< what to do: can be set to PROPAGATE, FINISH or POSTPROCESS */
   )
{
   char* command;
   int commandlen;
   SCIP_Bool endoffile;
   SCIP_Result result;
   DEC_Detector** detectors;
   int ndetectors;
   int i, j;
   SEEED_PROPAGATION_DATA* seeedPropData;
   Seeedpool* seeedpool;
   SCIP_Bool finished, displayinfo;
   char stri[SCIP_MAXSTRLEN];
   const char* actiontype;

   /* set string for dialog */
   if( action == PROPAGATE )
     actiontype = "propagated";
   else if( action == FINISH )
      actiontype = "finished";
   else if( action == POSTPROCESS )
      actiontype = "postprocessed";
   else
      actiontype = "UNDEFINED_ACTION";

   if( action == POSTPROCESS && userseeed->isComplete() == FALSE )
   {
      SCIPinfoMessage(scip, NULL, "The currently selected seeed is not finished, postprocessing not possible.\n");
      return SCIP_OKAY;
   }

   if( SCIPconshdlrDecompGetNDetectors(scip) == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No detector available!\n\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &detectors, SCIPconshdlrDecompGetNDetectors(scip)) );

   /* determine the detectors that implement the specified callback */
   ndetectors = 0;
   for( i = 0; i < SCIPconshdlrDecompGetNDetectors(scip); ++i )
   {
      if( (action == PROPAGATE && SCIPconshdlrDecompGetDetectors(scip)[i]->propagateFromToolbox)
       || (action == FINISH && SCIPconshdlrDecompGetDetectors(scip)[i]->finishFromToolbox)
       || (action == POSTPROCESS && SCIPconshdlrDecompGetDetectors(scip)[i]->postprocessSeeed) )
      {
         detectors[ndetectors] = SCIPconshdlrDecompGetNDetectors(scip)[i];
         ++ndetectors;
      }
   }

   if( ndetectors == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No detector implements this callback, returning!\n\n");
      return SCIP_OKAY;
   }

   /* build seeed propagation data needed in callbacks */
   seeedpool = userseeed->getSeeedpool();

   seeedPropData = new SEEED_PROPAGATION_DATA();
   seeedPropData->seeedpool = seeedpool;
   seeedPropData->nNewSeeeds = 0;
   seeedPropData->seeedToPropagate = new Seeed(userseeed);
   seeedPropData->seeedToPropagate->setSeeedpool(seeedpool);
   if( action != POSTPROCESS )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropData->newSeeeds), 1) );
      seeedPropData->newSeeeds[0] = NULL;
   }

   /* user dialog to select wanted detector, apply it and handle the returned seeeds, if any */
   finished = FALSE;
   while( !finished )
   {
      result = SCIP_DIDNOTFIND;
      /* list the detectors implementing the specified callback by name with a leading number */
      j = 1;
      SCIPinfoMessage(scip, NULL, "Available detectors:\n");
      for( i = 0; i < ndetectors; ++i )
      {
         SCIPinfoMessage(scip, NULL, "%d)", j);
         SCIPinfoMessage(scip, NULL, "%s\n", detectors[i]->name);
         ++j;
      }
      commandlen = 0;
      while( commandlen == 0 )
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
            "Type in the name or number of the detector that you want to use (or \"none\"): \nGCG/toolbox> ", &command,
            &endoffile) );
         commandlen = strlen(command);
      }

      if( !strncmp( command, "none", commandlen) == 0 && !strncmp( command, "quit", commandlen) == 0 )
      {
         for( i = 0; i < ndetectors; ++i )
         {
            sprintf(stri, "%d", i+1); //used for matching numberings in the list, off-by-one since detectors start with 0
            if( strncmp( command, detectors[i]->name, commandlen) == 0 || strncmp( command, stri, commandlen ) == 0 )
            {
               if( action == PROPAGATE )
                  SCIP_CALL( detectors[i]->propagateFromToolbox(scip, detectors[i], seeedPropData, &result, dialoghdlr, dialog) );
               else if( action == FINISH )
                  SCIP_CALL( detectors[i]->finishFromToolbox(scip, detectors[i], seeedPropData, &result, dialoghdlr, dialog) );
               else if( action == POSTPROCESS )
                  SCIP_CALL( detectors[i]->postprocessSeeed(scip, detectors[i], seeedPropData, &result) );
               break;
            }
         }
      }
      else
      {
         finished = TRUE;
         continue;
      }
      if( result == SCIP_SUCCESS )
      {
         if( action != POSTPROCESS )
         {
            SCIPinfoMessage(scip, NULL, "Considering implicits of newly found seeed(s)...\n");
            for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
            {
               assert(seeedPropData->newSeeeds[i] != NULL);
               seeedPropData->newSeeeds[i]->considerImplicits( ); //There may be open vars/cons left that were not matched
            }

            SCIPinfoMessage(scip, NULL, "\nSeeed was successfully %s, %d potentially new seeed(s) found.\n", actiontype,
               seeedPropData->nNewSeeeds);

            displayinfo = TRUE;
            if( seeedPropData->nNewSeeeds > 1 )
            {
               commandlen = 0;
               while( commandlen == 0 )
               {
                  SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
                     "More than one seeed found. Do you want to display information about all found seeeds anyway? (\"yes\"/\"no\")?\nGCG/toolbox> ",
                     &command, &endoffile) );
                  commandlen = strlen(command);
               }
               if( strncmp( command, "no", commandlen) == 0 )
               {
                  displayinfo = FALSE;
               }
               else if( strncmp( command, "quit", commandlen) == 0 )
               {
                  finished = TRUE;
                  continue;
               }
            }

            if( displayinfo )
            {
               for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
               {
                  seeedPropData->newSeeeds[i]->displayInfo( 0 );
               }
            }

            if( seeedPropData->nNewSeeeds == 1 )
            {
               commandlen = 0;
               while( commandlen == 0 )
               {
                  SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
                     "Do you want to visualize the new seeed (\"yes\"/\"no\")?\nGCG/toolbox> ", &command, &endoffile) );
                  commandlen = strlen(command);
               }
               if( strncmp( command, "yes", commandlen) == 0 )
               {
                  SCIP_CALL( SCIPdialogSelectVisualize(scip, dialoghdlr, dialog ) );
               }
               else if( strncmp( command, "quit", commandlen) == 0 )
               {
                  finished = TRUE;
                  continue;
               }
            }

            SCIPinfoMessage(scip, NULL, "\nSaving newly found seeeds...\n\n");
            for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
            {
               userseeed = new Seeed( seeedPropData->newSeeeds[i] ); //@todo delete previous seeed beforehand?
               SCIP_CALL( GCGtoolboxUserSeeedFlush(scip, userseeed) );
               assert(userseeed == NULL);
            }

            if( seeedPropData->nNewSeeeds == 1 )
            {
               commandlen = 0;
               while( commandlen == 0 )
               {
                  SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
                     "\nDo you want to continue the decomposition with the new Seeed (\"continue\"), or continue with the previous Seeed (\"previous\")?\nGCG/toolbox> ",
                     &command, &endoffile) );
                  commandlen = strlen(command);
               }
               if( strncmp( command, "continue", commandlen) == 0 )
               {
                  userseeed = new Seeed(seeedPropData->newSeeeds[0]); //@todo delete previous seeed beforehand?
               }
               else
               {
                  userseeed = new Seeed(seeedPropData->seeedToPropagate); //@todo delete previous seeed beforehand?
               }
            }
            else
            {
               userseeed = new Seeed(seeedPropData->seeedToPropagate); //@todo delete previous seeed beforehand?
            }
            finished = TRUE;
            continue;
         }
         else if( action == POSTPROCESS )
         {
            SCIPinfoMessage(scip, NULL, "\nSeeed successfully %s. %d seeed(s) found in the process.\n", actiontype,
               seeedPropData->nNewSeeeds);

            commandlen = 0;
            while( commandlen == 0 )
            {
               SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
                  "Do you want to save all found seeeds (\"all\") or none (\"none\")?\nGCG/toolbox> ", &command, &endoffile) );
               commandlen = strlen(command);
            }
            if( strncmp(command, "all", commandlen) == 0 )
            {
               SCIPinfoMessage(scip, NULL, "Storing seeeds...\n");
               for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
               {
                  userseeed =  new Seeed(seeedPropData->newSeeeds[i]); //@todo delete previous seeed beforehand?
                  SCIP_CALL( GCGtoolboxUserSeeedFlush(scip, userseeed) );
               }
               userseeed = new Seeed(seeedPropData->seeedToPropagate); //@todo delete previous seeed beforehand?
               SCIPinfoMessage(scip, NULL, "\nAll seeeds stored successfully!\n");
            }
            finished = TRUE;
            continue;
         }
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "Seeed could not be %s.\n", actiontype);

         commandlen = 0;
         while( commandlen == 0 )
         {
            SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
               "Do you want to select another detector (\"detector\") or return to the previous menu (\"previous\")?\nGCG/toolbox> ",
               &command, &endoffile) );
            commandlen = strlen(command);
         }
         if( strncmp( command, "detector", commandlen) == 0 )
         {
            continue;
         }
         else
         {
            finished = TRUE;
            continue;
         }
      }
   }

   SCIPfreeMemoryArrayNull( scip, &(seeedPropData->newSeeeds) );
   delete seeedPropData->seeedToPropagate;
   seeedPropData->newSeeeds = NULL;
   seeedPropData->nNewSeeeds = 0;
   delete seeedPropData;

   SCIPfreeBufferArray(scip, &detectors);
   return SCIP_OKAY;
}


/** Finishes a seeed created/modified in the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPdialogToolboxFinishSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   Seeed*                  seeed       /**< seeed to finish */
   )
{
   return SCIPdialogToolboxActOnSeeed(scip, dialoghdlr, dialog, seeed, FINISH);
}


/** Propagates a seeed created/modified in the toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxPropagateSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   Seeed*                  seeed       /**< seeed to propagate */
   )
{
   return SCIPdialogToolboxActOnSeeed(scip, dialoghdlr, dialog, seeed, PROPAGATE);
}


/** Postprocesses a seeed created/modified in the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPdialogToolboxPostprocessSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   Seeed*                  seeed       /**< seeed to postprocess */
   )
{
   return SCIPdialogToolboxActOnSeeed(scip, dialoghdlr, dialog, seeed, POSTPROCESS);
}


/**
 * @brief method to handle and moderate user input for modifying decompositions
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPdialogExecToolboxModify(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   int                     startindex  /**< index in seeed list to start list extract at */
   )
{
   SCIP_Bool finished;
   char* command;
   SCIP_Bool endoffile;
   int commandlen;
   SCIP_Bool selectedsomeseeed;
   int nseeeds;

   selectedsomeseeed = TRUE;

   assert(scip != NULL);
   finished = FALSE;

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( SCIPconshdlrDecompGetNSeeedLeafs(scip) == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( SCIPtransformProb(scip) );
      SCIPinfoMessage(scip, NULL, "Applied tranformation to problem.\n");
   }
   /* 1) update list of interesting seeeds */
   SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );

   /* 2) while user has not aborted: show current list extract */
   while ( !finished )
   {
      SCIP_CALL( SCIPdialogShowListExtractHeader(scip) );

      SCIP_CALL( SCIPdialogShowListExtract(scip, startindex) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please choose an existing partial decomposition for modification (type \"choose <id>\" or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen = strlen(command);

      nseeeds = SCIPconshdlrDecompGetNSeeedLeafs(scip);
      if( strncmp( command, "back", commandlen) == 0 )
      {
         GCGsetSelectFirstIdToVisu(scip, GCGgetSelectFirstIdToVisu(scip) - DEFAULT_MENULENGTH);
         if(GCGgetSelectFirstIdToVisu(scip) < 0 )
            GCGsetSelectFirstIdToVisu(scip, 0);
         continue;
      }
      if( strncmp( command, "next", commandlen) == 0 )
      {
         GCGsetSelectFirstIdToVisu(scip, (GCGgetSelectFirstIdToVisu(scip) + DEFAULT_MENULENGTH));
         if( GCGgetSelectFirstIdToVisu(scip) > nseeeds - DEFAULT_MENULENGTH )
            GCGsetSelectFirstIdToVisu(scip, nseeeds - DEFAULT_MENULENGTH);
         continue;
      }
      if( strncmp( command, "top", commandlen) == 0 )
      {
         GCGsetSelectFirstIdToVisu(scip, 0);
         continue;
      }
      if( strncmp( command, "end", commandlen) == 0 )
      {
         GCGsetSelectFirstIdToVisu(scip, nseeeds - DEFAULT_MENULENGTH);
         continue;
      }

      if( strncmp( command, "quit", commandlen) == 0 )
      {
         finished = TRUE;
         selectedsomeseeed = FALSE;
         continue;
      }

      if( strncmp( command, "choose", commandlen) == 0 )
      {
         SCIP_RETCODE retcode = SCIPdialogToolboxChoose(scip, dialoghdlr, dialog );
    if (retcode != SCIP_OKAY)
    {
       selectedsomeseeed = FALSE;
       continue;
    }
    else
    {
       selectedsomeseeed = TRUE;
       finished = TRUE;
       break;
    }
      }

      if( strncmp( command, "abort", commandlen) == 0 )
      {
         finished = TRUE;
         selectedsomeseeed = FALSE;
         continue;
      }

      if( strncmp( command, "help", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogShowHelp(scip) );
         continue;
      }

      if( strncmp( command, "visualize", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogSelectVisualize(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "propagate", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "finishseeed", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxFinishSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "postprocess", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPostprocessSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
   }
   finished = FALSE;
   while ( !finished && selectedsomeseeed )
   {
      int commandlen2;
      SCIP_Bool success;

      userseeed.displaySeeed();

      SCIP_CALL( SCIPdialogShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen2 = strlen(command);

      /* case distinction: */
      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPdialogToolboxModifyConss(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPdialogToolboxModifyVars(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "finish", commandlen2) == 0 )
      {
         SCIPdialogToolboxModifyFinish(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         if(GCGgetSelectLastUserSeeed(scip) != NULL)
            delete GCGgetSelectLastUserSeeed(scip);
         GCGsetSelectLastUserSeeed(scip, new Seeed(GCGgetSelectCurrUserSeeed(scip)));
         GCGgetSelectCurrUserSeeed(scip)->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         /* determine whether there is a seeedpool for the presolved problem */
         SEEED_WRAPPER* wrapper = SCIPconshdlrDecompGetSeeedpool(scip);
         Seeedpool* internalseeedpool = wrapper->seeedpool;

         Seeedpool* seeedpool;
         Seeed* curruserseeed = GCGgetSelectCurrUserSeeed(scip);

         if( !curruserseeed->isFromUnpresolved() && internalseeedpool == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = curruserseeed->getSeeedpool();
         if( seeedpool == NULL )

         curruserseeed->sort();
         curruserseeed->considerImplicits();
         curruserseeed->calcHashvalue();
         assert( curruserseeed->checkConsistency() );

         if( curruserseeed->isComplete() )
         {
            seeedpool->addSeeedToFinished(curruserseeed, &success);
            if( !success )
            {
               delete curruserseeed;
            }
         } else
         {
            seeedpool->addSeeedToIncomplete(curruserseeed, &success);
            if( !success )
            {
               delete curruserseeed;
            }
         }
         curruserseeed = NULL;
         GCGsetSelectCurrUserSeeed(scip, curruserseeed);
         finished = TRUE;

         continue;
      }

      if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( GCGgetSelectLastUserSeeed(scip) == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete GCGgetSelectCurrUserSeeed(scip);
            GCGsetSelectCurrUserSeeed(scip, GCGgetSelectLastUserSeeed(scip));
            GCGsetSelectLastUserSeeed(scip, NULL);
         }
         continue;
      }

      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         userseeed.showVisualisation();
         continue;
      }

      if( strncmp( command, "propagate", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
      if( strncmp( command, "finishseeed", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxFinishSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "postprocess", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPostprocessSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
   }
   return SCIP_OKAY;
}

/**
 * @brief method to handle and moderate user input for creating new decompositions by the user
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return data structure
 */
static
SCIP_RETCODE SCIPdialogExecToolboxCreate(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   char* command;
   SCIP_Bool endoffile;
   SCIP_Bool         finished;
   int commandlen;

   assert(scip != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( SCIPconshdlrDecompGetNSeeedLeafs(scip) == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( SCIPtransformProb(scip) );
      SCIPinfoMessage(scip, NULL, "Applied tranformation to problem.\n");
   }

   /* create new decomposition */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ", &command, &endoffile) );
   commandlen = strlen(command);

   Seeedpool* seeedpool;
   SCIP_Bool isfromunpresolved;

   while( (strncmp( command, "presolved", commandlen) != 0 && strncmp( command, "unpresolved", commandlen) != 0) || commandlen == 0)
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Invalid input. Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ", &command, &endoffile) );
      commandlen = strlen(command);
   }

   if( strncmp( command, "presolved", commandlen) == 0 )
   {
      isfromunpresolved = FALSE;
      if (SCIPconshdlrDecompGetSeeedpool(scip) != NULL )
         seeedpool = SCIPconshdlrDecompGetSeeedpool(scip);
      else
      {
         if( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
         {
            SCIPinfoMessage(scip, NULL, "Problem is not presolved yet. Please presolve it first!\n");
            return SCIP_OKAY;
         }

         SCIPconshdlrDecompCreateSeeedpool(scip);
         seeedpool = SCIPconshdlrDecompGetSeeedpool(scip);
      }
   }
   else
   {
      isfromunpresolved = TRUE;
      if ( SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip) == NULL )
         SCIPconshdlrDecompCreateSeeedpoolUnpresolved(scip);
      seeedpool = SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip);

   }
   if( seeedpool == NULL )
   {
      if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED )

      {
         if(SCIPconshdlrDecompGetSeeedpool(scip) == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);
         seeedpool = SCIPconshdlrDecompGetSeeedpool(scip);
      }
      else
      {
         if(SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip) == NULL)
            SCIPconshdlrDecompCreateSeeedpoolUnpresolved(scip);
         seeedpool = SCIPconshdlrDecompGetSeeedpoolUnpresolved(scip);
      }

   }

   Seeed* newseeed = new Seeed( scip, SCIPconshdlrDecompGetNextSeeedID(scip), seeedpool );
   newseeed->setIsFromUnpresolved(isfromunpresolved);

   finished = FALSE;
   while ( !finished )
   {
      int commandlen2;
      SCIP_Bool success;

      newseeed->displaySeeed();
      SCIP_CALL( SCIPdialogShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen2 = strlen(command);

      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPdialogToolboxModifyConss(scip, dialoghdlr, dialog, newseeed);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPdialogToolboxModifyVars(scip, dialoghdlr, dialog, newseeed);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         //@todo here was update of last user seeed
         newseeed->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         if( !newseeed->isFromUnpresolved() && SCIPconshdlrDecompGetSeeedpool(scip) == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = newseeed->getSeeedpool();
         assert( seeedpool != NULL );

         newseeed->sort();
         newseeed->considerImplicits();
         newseeed->calcHashvalue();
         assert( newseeed->checkConsistency() );

         if(newseeed->isComplete())
         {
            seeedpool->addSeeedToFinished(newseeed, &success);
            if( !success )
            {
               delete newseeed;
            }
         }
         else
         {
            seeedpool->addSeeedToIncomplete(newseeed, &success);
            if( !success )
            {
               delete newseeed;
            }
         }

         finished = TRUE;
         continue;
      }

      //@todo reintroduce, maybe hand down pointer to last one?
      /*if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( GCGgetSelectLastUserSeeed(scip) == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete GCGgetSelectCurrUserSeeed(scip);
            GCGsetSelectCurrUserSeeed(scip, GCGgetSelectLastUserSeeed(scip));
            GCGsetSelectLastUserSeeed(scip, NULL);
         }
         continue;
      }*/

      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         newseeed->showVisualisation();
         continue;
      }

      if( strncmp( command, "propagate", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPropagateSeeed(scip, dialoghdlr, dialog, newseeed) );
         continue;
      }

      if( strncmp( command, "finish", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxFinishSeeed(scip, dialoghdlr, dialog, newseeed) );
         continue;
      }

      if( strncmp( command, "postprocess", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPostprocessSeeed(scip, dialoghdlr, dialog, newseeed) );
         continue;
      }
   }
   return SCIP_OKAY;
}


/** Lets the user select decompositions from the explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogExploreSelect(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;
   std::vector<int>* selected;

   int commandlen;

   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be selected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   Seeed_Wrapper sw;
   SCIP_CALL( GCGgetSeeedFromID(scip, &idtovisu, &sw) );

   if( sw.seeed == NULL )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_PARAMETERWRONGVAL;
   }

   sw.seeed->setSelected(!sw.seeed->isSelected() );

   if( sw.seeed->isSelected() )
      std::cout << "is selected!" << sw.seeed->isSelected() <<std::endl;

   return SCIP_OKAY;
}


/*
 * @brief method too handle user input for "explore" command
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPdialogExecSelect(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_Bool         finished;
   char* command;
   SCIP_Bool endoffile;
   int nseeeds;
   int startindex = 0;

   SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
   /* while user has not aborted: show current list extract */

   while ( !finished )
   {
      int commandlen;

      SCIP_CALL( SCIPdialogShowListExtractHeader(scip) );

      SCIP_CALL( SCIPdialogShowListExtract(scip, startindex) );


      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Please enter command or decomposition id to select (or \"h\" for help) : \nGCG/explore> ", &command, &endoffile) );

      commandlen = strlen(command);

      nseeeds = SCIPconshdlrDecompGetNSeeedLeafs(scip);
      if( strncmp( command, "back", commandlen) == 0 )
      {
         startindex = startindex - DEFAULT_MENULENGTH;
         if(startindex < 0 )
            startindex = 0;
         continue;
      }
      if( strncmp( command, "next", commandlen) == 0 )
      {
         startindex = startindex + DEFAULT_MENULENGTH;
         if( startindex > nseeeds - DEFAULT_MENULENGTH )
            startindex = nseeeds - DEFAULT_MENULENGTH;
         continue;
      }
      if( strncmp( command, "top", commandlen) == 0 )
      {
         startindex = 0;
         continue;
      }
      if( strncmp( command, "end", commandlen) == 0 )
      {
         startindex = nseeeds - DEFAULT_MENULENGTH;
         continue;
      }

      if( strncmp( command, "quit", commandlen) == 0 )
      {
         finished = TRUE;
         SCIP_CALL(SCIPconshdlrDecompChooseCandidatesFromSelected(scip, FALSE) );
         continue;
      }

      if( strncmp( command, "legend", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogShowLegend(scip) );
         continue;
      }

      if( strncmp( command, "help", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogShowHelp(scip) );
         continue;
      }

      if( strncmp( command, "visualize", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogSelectVisualize(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "inspect", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogSelectInspect( scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "calc_strong", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogSelectCalcStrongDecompositionScore( scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "select", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogExploreSelect(scip, dialoghdlr, dialog ) );
         continue;
      }
      if( strncmp( command, "modify", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogExecToolboxModify(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
      if( strncmp( command, "create", commandlen) == 0 )
      {
         SCIP_CALL( SCIPdialogExecToolboxCreate(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
   }

   return SCIP_OKAY;
}

} // namespace gcg
