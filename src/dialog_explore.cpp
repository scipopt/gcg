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

typedef gcg::Seeed* SeeedPtr;

struct SCIP_MenuData
{
};

/** Shows header for seeed information in explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogShowListExtractHeader(
   SCIP*                   scip  /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   int ndetectedpresolved;
   int ndetectedunpresolved;
   int nuserpresolvedfull;
   int nuserpresolvedpartial;
   int nuserunpresolvedfull;
   int nuserunpresolvedpartial;

   char* scorename;

   size_t i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   scorename = SCIPconshdlrDecompGetScoretypeShortName(scip, SCIPconshdlrdataGetScoretype(conshdlrdata) );

   ndetectedpresolved = 0;
   ndetectedunpresolved = 0;
   nuserpresolvedfull = 0;
   nuserpresolvedpartial = 0;
   nuserunpresolvedfull = 0;
   nuserunpresolvedpartial = 0;

   /* count corresponding seeeds */
   for ( i = 0; i < conshdlrdata->listall->size(); ++i )
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->listall->at(i);
      if( seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::NOT && !seeed->isFromUnpresolved() )
         ++ndetectedpresolved;
      if( seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::NOT && seeed->isFromUnpresolved() )
         ++ndetectedunpresolved;
      if( seeed->isComplete() && ( seeed->getUsergiven() == gcg::USERGIVEN::COMPLETE || seeed->getUsergiven() == gcg::USERGIVEN::COMPLETED_CONSTOMASTER) && !seeed->isFromUnpresolved() )
         ++nuserpresolvedfull;
      if( !seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::PARTIAL && !seeed->isFromUnpresolved() )
         ++nuserpresolvedpartial;
      if( seeed->isComplete() && ( seeed->getUsergiven() == gcg::USERGIVEN::COMPLETE || seeed->getUsergiven() == gcg::USERGIVEN::COMPLETED_CONSTOMASTER) && seeed->isFromUnpresolved() )
         ++nuserunpresolvedfull;
      if( !seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::PARTIAL && seeed->isFromUnpresolved() )
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

/** Shows information about the current user seeed in toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogShowCurrUserSeeedInfo(
   SCIP*                   scip     /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if ( conshdlrdata->curruserseeed->isFromUnpresolved() )
      conshdlrdata->curruserseeed->displaySeeed();
   else
      conshdlrdata->curruserseeed->displaySeeed();


   return SCIP_OKAY;
}

/** Shows detailed information about seeeds in explore menu
 *
 *@returns SCIP status
 * */
static
SCIP_RETCODE SCIPdialogShowListExtract(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   size_t i;

   for( i = conshdlrdata->startidvisu; i < (size_t) conshdlrdata->startidvisu + (size_t) conshdlrdata->selectvisulength && i < conshdlrdata->listall->size(); ++i)
   {
      SeeedPtr seeed;

      seeed = conshdlrdata->listall->at(i);

      assert( seeed->checkConsistency( ) );

      SCIPdialogMessage(scip, NULL, " %4d   ", i );
      SCIPdialogMessage(scip, NULL, "%5d  ", seeed->getNBlocks() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMasterconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNLinkingvars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMastervars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNTotalStairlinkingvars() );
      if( seeed->isComplete() )
         SCIPdialogMessage(scip, NULL, "%.4f  ",  seeed->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) );
      else
         SCIPdialogMessage(scip, NULL, "<=%.2f  ", seeed->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) );
      SCIPdialogMessage(scip, NULL, "%7s  ", seeed->getDetectorChainString() );
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->isFromUnpresolved() ? "no" : "yes")  );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenvars() );
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->getUsergiven() == gcg::USERGIVEN::NOT ? "no" : "yes")   );
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   char * scorename;
   char * scoredescr;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   scorename = SCIPconshdlrDecompGetScoretypeShortName(scip, SCIPconshdlrdataGetScoretype(conshdlrdata) );
   scoredescr = SCIPconshdlrDecompGetScoretypeDescription(scip, SCIPconshdlrdataGetScoretype(conshdlrdata) );


   SCIPdialogMessage(scip, NULL, "List of included detectors for decompositions histories: \n" );

   SCIPdialogMessage(scip, NULL, "\n%30s    %4s\n", "detector" , "char"  );
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "--------" , "----"  );

   for( int det = 0; det < conshdlrdata->ndetectors; ++det )
   {
      DEC_DETECTOR* detector;

      detector = conshdlrdata->detectors[det];

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

/** Modifies the number of presented seeeds in the explore menu via dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogModifyNVisualized(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int newval;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the maximum number of decompositions displayed at once in the table [%d]:\n",
      conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   newval = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      newval = atoi(ntovisualize);

   if (newval != 0)
      conshdlrdata->selectvisulength = newval;

   return SCIP_OKAY;
}

/** shows a visualization of current user seeed
 *
 * @returns SCip status*/
static
SCIP_RETCODE SCIPdialogSelectVisualizeCurrentUserSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->curruserseeed->showVisualisation();

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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be visualized:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0 )
      idtovisu = atoi(ntovisualize);

   /* check whether ID is in valid range */
   if( (int)conshdlrdata->listall->size() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( commandlen == 0 || idtovisu < 0 || idtovisu >= (int)conshdlrdata->listall->size() )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_OKAY;
   }

   conshdlrdata->listall->at(idtovisu)->showVisualisation();

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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntocalcstrong;
   SCIP_Bool endoffile;
   int idtocalcstrong;
   int commandlen;

   assert( scip != NULL );
   conshdlr = SCIPfindConshdlr( scip, CONSHDLR_NAME );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData( conshdlr );
   assert( conshdlrdata != NULL );

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
   if( 0 <= idtocalcstrong && idtocalcstrong < (int)conshdlrdata->listall->size() )
   {
      SCIP_Real score;
      gcg::Seeedpool* seeedpool = ( conshdlrdata->listall->at( idtocalcstrong )->isFromUnpresolved() ?
         conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool );
      seeedpool->calcStrongDecompositionScore(conshdlrdata->listall->at( idtocalcstrong ), &score);
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntoinspect;
   char* ndetaillevel;
   SCIP_Bool endoffile;
   int idtoinspect;
   int detaillevel;

   int commandlen;

   assert( scip != NULL );
   conshdlr = SCIPfindConshdlr( scip, CONSHDLR_NAME );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData( conshdlr );
   assert( conshdlrdata != NULL );

   /* read the id of the decomposition to be inspected */
   SCIPdialogMessage( scip, NULL, "Please specify the id of the decomposition to be inspected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntoinspect, &endoffile ) );
   commandlen = strlen( ntoinspect );

   idtoinspect = -1;
   if( commandlen != 0 )
      idtoinspect = atoi( ntoinspect );

   /* check whether ID is in valid range */
   if( idtoinspect < 0 || idtoinspect >= (int)conshdlrdata->listall->size() )
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
   assert( 0 <= idtoinspect && idtoinspect < (int)conshdlrdata->listall->size() );
   conshdlrdata->listall->at( idtoinspect )->displayInfo( detaillevel );

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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool         finished;
   char* command;
   SCIP_Bool endoffile;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   finished = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
   /* while user has not aborted: show current list extract */

   while ( !finished )
   {
      int commandlen;

      SCIP_CALL( SCIPdialogShowListExtractHeader(scip) );

      SCIP_CALL( SCIPdialogShowListExtract(scip) );


      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Please enter command or decomposition id to select (or \"h\" for help) : \nGCG/explore> ", &command, &endoffile) );

      commandlen = strlen(command);

      if( strncmp( command, "back", commandlen) == 0 )
      {
         conshdlrdata->startidvisu -= conshdlrdata->selectvisulength;
         if(conshdlrdata->startidvisu < 0 )
            conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "next", commandlen) == 0 )
      {
         conshdlrdata->startidvisu += conshdlrdata->selectvisulength;
         if( conshdlrdata->startidvisu > (int) conshdlrdata->listall->size() - conshdlrdata->selectvisulength )
            conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
         continue;
      }
      if( strncmp( command, "top", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "end", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
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
         SCIP_CALL(SCIPconshdlrDecompShowLegend(scip) );
         continue;
      }

      if( strncmp( command, "dispNEntries", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogModifyNVisualized(scip, dialoghdlr, dialog) );
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
         SCIP_CALL( SCIPconshdlrDecompSelectInspect( scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "calc_strong", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompSelectCalcStrongDecompositionScore( scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "select", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompExploreSelect(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "toolbox", commandlen) == 0 )
      {
    //@todo deprecated, use create/modify instead
         SCIP_CALL( SCIPconshdlrDecompExecToolbox(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
      if( strncmp( command, "modify", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompExecToolboxModify(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
      if( strncmp( command, "create", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompExecToolboxCreate(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
   }

   return SCIP_OKAY;
}


/** Lets user modify conss during modification of seeed in toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxModifyConss(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;
    SCIP_Bool         matching;
    char* consregex;
    char* command;
    char* command2;
    SCIP_Bool endoffile;
    int commandlen;

    assert(scip != NULL);
    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    assert( conshdlr != NULL );
    matching = FALSE;


    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != NULL);

    SeeedPtr seeed  = conshdlrdata->curruserseeed;
    gcg::Seeedpool* seeedpool;
    std::vector<int> matchingconss  = std::vector<int>(0);

    seeedpool = seeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
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

    for( int oc = 0; oc < seeed->getNOpenconss(); ++oc )
    {
       const char* consname;

       consname = SCIPconsGetName(  seeedpool->getConsForIndex(seeed->getOpenconss()[oc] ) );


       if( std::regex_match(consname, expr) )
       {
          matching = TRUE;
          matchingconss.push_back(seeed->getOpenconss()[oc]);
          SCIPdebugMessage(" consname %s matches regex %s \n", consname, consregex );
       } else
          SCIPdebugMessage(" consname %s does not match regex %s \n", consname, consregex);
    }

    if( !matching )
    {
       SCIPdialogMessage(scip, NULL, " There are no unassigned constraints with names matching given regular expression. Return to toolbox main menu.\n");
       return SCIP_OKAY;
    }

    if( conshdlrdata->lastuserseeed != NULL)
       delete conshdlrdata->lastuserseeed;
    conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;


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
          seeed->bookAsMasterCons( matchingconss[mc] );
       }
    }
    else if( strncmp( command, "block", commandlen) == 0 )
    {
       SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please specify the block number these constraints should be assigned to: \nGCG/toolbox> ", &command2, &endoffile) );
       char* tail;
       int blockid = strtol(command2, &tail, 10);
       for( size_t mc = 0 ;  mc < matchingconss.size(); ++mc )
       {
          seeed->bookAsBlockCons( matchingconss[mc], blockid );
       }
    }
    else
       return SCIP_OKAY;

    seeed->flushBooked();

   return SCIP_OKAY;
}

/** Lets user specify how to finish the modified seeed while using the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPdialogToolboxModifyFinish(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool         choosenfinisher;

   char* command;
   SCIP_Bool endoffile;
   char* tail;
   int finisherid;
   SEEED_PROPAGATION_DATA* seeedPropData;
   DEC_DETECTOR* finisher;
   SCIP_Result result;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SeeedPtr seeed  = conshdlrdata->curruserseeed;
   gcg::Seeedpool* seeedpool;
   std::vector<int> matchingvars  = std::vector<int>(0);

   seeedpool = seeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
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
   seeedPropData->seeedToPropagate = new gcg::Seeed(conshdlrdata->curruserseeed);

   if( conshdlrdata->lastuserseeed != NULL)
      delete conshdlrdata->lastuserseeed;
   conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;

   finisher = seeedpool->getFinishingDetectorForIndex(finisherid);
   finisher->finishSeeed(scip, finisher, seeedPropData, &result);

   delete conshdlrdata->curruserseeed;

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
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxChoose(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntochoose;
   SCIP_Bool endoffile;
   int idtochoose;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the (partial) decomposition to be chosen for modification:\n",
      conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntochoose, &endoffile) );
   commandlen = strlen(ntochoose);

   idtochoose = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      idtochoose = atoi(ntochoose);

   if ( commandlen == 0 || idtochoose < 0 || idtochoose >= (int)conshdlrdata->listall->size() )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_PARAMETERWRONGVAL;
   }

   if( conshdlrdata->curruserseeed != NULL )
      delete conshdlrdata->curruserseeed;

   conshdlrdata->curruserseeed = new gcg::Seeed( conshdlrdata->listall->at(idtochoose) );

   return SCIP_OKAY;
}

/** Lets user modify vars during use of the toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxModifyVars(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{

   SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;
    SCIP_Bool         matching;
    char* varregex;
    char* command;
    char* command2;
    SCIP_Bool endoffile;
    int commandlen;

    assert(scip != NULL);
    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    assert( conshdlr != NULL );
    matching = FALSE;


    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != NULL);

    SeeedPtr seeed  = conshdlrdata->curruserseeed;
    gcg::Seeedpool* seeedpool;
    std::vector<int> matchingvars  = std::vector<int>(0);

    seeedpool = seeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
    /* Does user want to modify existing or create a new partial decomposition ?*/
    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
       "Please specify a regular expression (modified ECMAScript regular expression grammar) matching the names of unassigned variables you want to assign : \nGCG/toolbox> ",
       &varregex, &endoffile) );

    /* case distinction: */

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

    for( int oc = 0; oc < seeed->getNOpenvars(); ++oc )
    {
       const char* varname;

       varname = SCIPvarGetName(  seeedpool->getVarForIndex(seeed->getOpenvars()[oc] ) );

       SCIPdebugMessage("check var %s for regex %s \n", varname, varregex);

       if( std::regex_match(varname, expr) )
       {
          matching = TRUE;
          matchingvars.push_back(seeed->getOpenvars()[oc]);
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

    if( conshdlrdata->lastuserseeed != NULL)
       delete conshdlrdata->lastuserseeed;
    conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;


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

    /* case distinction: */
    if( strncmp( command, "master", commandlen) == 0 )
    {
       for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
       {
          seeed->bookAsMasterVar( matchingvars[mc] );
       }
    } else
       if( strncmp( command, "linking", commandlen) == 0 )
           {
              for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
              {
                 seeed->bookAsLinkingVar( matchingvars[mc] );
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
          seeed->bookAsBlockVar( matchingvars[mc], blockid );
       }
    }
    else
       return SCIP_OKAY;

    seeed->flushBooked();
    seeed->deleteEmptyBlocks(true);

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
   toolboxtype             action      /**< what to do: can be set to PROPAGATE, FINISH or POSTPROCESS */
   )
{
   char* command;
   int commandlen;
   SCIP_Bool endoffile;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Result result;
   DEC_Detector** detectors;
   int ndetectors;
   int i, j;
   SEEED_PROPAGATION_DATA* seeedPropData;
   gcg::Seeedpool* seeedpool;
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

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   if( action == POSTPROCESS && conshdlrdata->curruserseeed->isComplete() == FALSE )
   {
      SCIPinfoMessage(scip, NULL, "The currently selected seeed is not finished, postprocessing not possible.\n");
      return SCIP_OKAY;
   }

   if( conshdlrdata->ndetectors == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No detector available!\n\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &detectors, conshdlrdata->ndetectors) );

   /* determine the detectors that implement the specified callback */
   ndetectors = 0;
   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      if( (action == PROPAGATE && conshdlrdata->detectors[i]->propagateFromToolbox)
       || (action == FINISH && conshdlrdata->detectors[i]->finishFromToolbox)
       || (action == POSTPROCESS && conshdlrdata->detectors[i]->postprocessSeeed) )
      {
         detectors[ndetectors] = conshdlrdata->detectors[i];
         ++ndetectors;
      }
   }

   if( ndetectors == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No detector implements this callback, returning!\n\n");
      return SCIP_OKAY;
   }

   /* build seeed propagation data needed in callbacks */
   seeedpool = conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;

   seeedPropData = new SEEED_PROPAGATION_DATA();
   seeedPropData->seeedpool = seeedpool;
   seeedPropData->nNewSeeeds = 0;
   seeedPropData->seeedToPropagate = new gcg::Seeed(conshdlrdata->curruserseeed);
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
               conshdlrdata->curruserseeed = new gcg::Seeed( seeedPropData->newSeeeds[i] );
               SCIP_CALL( SCIPconshdlrDecompUserSeeedFlush(scip) );
               assert(conshdlrdata->curruserseeed == NULL);
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
                  conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->newSeeeds[0]);
               }
               else
               {
                  conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->seeedToPropagate);
               }
            }
            else
            {
               conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->seeedToPropagate);
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
                  conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->newSeeeds[i]);
                  SCIP_CALL( SCIPconshdlrDecompUserSeeedFlush(scip) );
               }
               conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->seeedToPropagate);
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
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   return SCIPconshdlrDecompToolboxActOnSeeed(scip, dialoghdlr, dialog, FINISH);
}


/** Propagates a seeed created/modified in the toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPdialogToolboxPropagateSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   return SCIPconshdlrDecompToolboxActOnSeeed(scip, dialoghdlr, dialog, PROPAGATE);
}


/** Postprocesses a seeed created/modified in the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPdialogToolboxPostprocessSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   return SCIPconshdlrDecompToolboxActOnSeeed(scip, dialoghdlr, dialog, POSTPROCESS);
}


/*
 * @brief method to handle and moderate user input for modifying decompositions
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPdialogExecToolboxModify(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool         finished;
   char* command;
   SCIP_Bool endoffile;
   int commandlen;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool selectedsomeseeed;

   selectedsomeseeed = TRUE;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   finished = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( (int)conshdlrdata->listall->size() == 0 )
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

      SCIP_CALL( SCIPdialogShowListExtract(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please choose an existing partial decomposition for modification (type \"choose <id>\" or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen = strlen(command);

      /* case distinction: */
      if( strncmp( command, "back", commandlen) == 0 )
      {
         conshdlrdata->startidvisu -= conshdlrdata->selectvisulength;
         if(conshdlrdata->startidvisu < 0 )
            conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "next", commandlen) == 0 )
      {
         conshdlrdata->startidvisu += conshdlrdata->selectvisulength;
         if( conshdlrdata->startidvisu > (int) conshdlrdata->listall->size() - conshdlrdata->selectvisulength )
            conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
         continue;
      }
      if( strncmp( command, "top", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "end", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
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

      if( strncmp( command, "change number displayed", commandlen) == 0 )
      {
         SCIP_CALL(SCIPdialogModifyNVisualized(scip, dialoghdlr, dialog) );
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

      SCIP_CALL( SCIPconshdlrDecompShowCurrUserSeeedInfo(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen2 = strlen(command);

      /* case distinction: */
      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyConss(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyVars(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "finish", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyFinish(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         if( conshdlrdata->lastuserseeed != NULL)
            delete conshdlrdata->lastuserseeed;
         conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;
         conshdlrdata->curruserseeed->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         gcg::Seeedpool* seeedpool;
         if( !conshdlrdata->curruserseeed->isFromUnpresolved() && conshdlrdata->seeedpool == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = ( conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool);
         if( seeedpool == NULL )

         conshdlrdata->curruserseeed->sort();
         conshdlrdata->curruserseeed->considerImplicits();
         conshdlrdata->curruserseeed->calcHashvalue();
         assert( conshdlrdata->curruserseeed->checkConsistency() );



         if( conshdlrdata->curruserseeed->isComplete() )
         {
            seeedpool->addSeeedToFinished(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         } else
         {
            seeedpool->addSeeedToIncomplete(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         }
         conshdlrdata->curruserseeed = NULL;
         finished = TRUE;


         continue;
      }

      if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( conshdlrdata->lastuserseeed == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete conshdlrdata->curruserseeed;
            conshdlrdata->curruserseeed = conshdlrdata->lastuserseeed;
            conshdlrdata->lastuserseeed = NULL;
         }
         continue;
      }


      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         SCIP_CALL(SCIPdialogSelectVisualizeCurrentUserSeeed(scip, dialoghdlr, dialog ) );
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

/*
 * @brief method to handle and moderate user input for creating new decompositions by the user
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return data structure
 */
SCIP_RETCODE SCIPdialogExecToolboxCreate(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   char* command;
   SCIP_Bool endoffile;
   SCIP_Bool         finished;
   int commandlen;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( (int)conshdlrdata->listall->size() == 0 )
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

   if( conshdlrdata->curruserseeed != NULL )
      delete conshdlrdata->curruserseeed;

   gcg::Seeedpool* seeedpool;
   SCIP_Bool isfromunpresolved;

   while( (strncmp( command, "presolved", commandlen) != 0 && strncmp( command, "unpresolved", commandlen) != 0) || commandlen == 0)
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Invalid input. Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ", &command, &endoffile) );
      commandlen = strlen(command);
   }

   /* case distinction: */
   if( strncmp( command, "presolved", commandlen) == 0 )
   {
      isfromunpresolved = FALSE;
      if (conshdlrdata->seeedpool != NULL )
         seeedpool = conshdlrdata->seeedpool;
      else
      {
         if( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
         {
            SCIPinfoMessage(scip, NULL, "Problem is not presolved yet. Please presolve it first!\n");
            return SCIP_OKAY;
         }

         conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
         seeedpool = conshdlrdata->seeedpool;
      }
   }
   else
   {
      isfromunpresolved = TRUE;
      if ( conshdlrdata->seeedpoolunpresolved == NULL )
         conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));
       seeedpool = conshdlrdata->seeedpoolunpresolved;

   }
   if( seeedpool == NULL )
   {
      if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED )

      {
         if (conshdlrdata->seeedpool == NULL )
            conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
         seeedpool = conshdlrdata->seeedpool;
      }
      else
      {
         if ( conshdlrdata->seeedpoolunpresolved == NULL)
            conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));
         seeedpool = conshdlrdata->seeedpoolunpresolved;
      }

   }
   conshdlrdata->curruserseeed = new gcg::Seeed( scip, SCIPconshdlrDecompGetNextSeeedID(scip), seeedpool );
   conshdlrdata->curruserseeed->setIsFromUnpresolved(isfromunpresolved);
   finished = FALSE;
   while ( !finished )
   {
      int commandlen2;
      SCIP_Bool success;

      SCIP_CALL( SCIPconshdlrDecompShowCurrUserSeeedInfo(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen2 = strlen(command);

      /* case distinction: */
      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyConss(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyVars(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         if( conshdlrdata->curruserseeed->isFromUnpresolved() )
            seeedpool = conshdlrdata->seeedpoolunpresolved;
         else
            seeedpool = conshdlrdata->seeedpool;
         if( conshdlrdata->lastuserseeed != NULL)
            delete conshdlrdata->lastuserseeed;
         conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;
         conshdlrdata->curruserseeed->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         if( !conshdlrdata->curruserseeed->isFromUnpresolved() && conshdlrdata->seeedpool == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = ( conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool);
         if( seeedpool == NULL )

         conshdlrdata->curruserseeed->sort();
         conshdlrdata->curruserseeed->considerImplicits();
         conshdlrdata->curruserseeed->calcHashvalue();
         assert( conshdlrdata->curruserseeed->checkConsistency() );



         if( conshdlrdata->curruserseeed->isComplete() )
         {
            seeedpool->addSeeedToFinished(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         } else
         {
            seeedpool->addSeeedToIncomplete(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         }
         conshdlrdata->curruserseeed = NULL;
         finished = TRUE;


         continue;
      }

      if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( conshdlrdata->lastuserseeed == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete conshdlrdata->curruserseeed;
            conshdlrdata->curruserseeed = conshdlrdata->lastuserseeed;
            conshdlrdata->lastuserseeed = NULL;
         }
         continue;
      }


      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         SCIP_CALL(SCIPdialogSelectVisualizeCurrentUserSeeed(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "propagate", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "finish", commandlen2) == 0 )
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

/*
 * method to handle and moderate user input for creating new decompositions
 * and modifying existing decompositions by the user
 *
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPdialogExecToolbox(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool         finished;
   char* command;
   SCIP_Bool endoffile;
   int commandlen;
   SCIP_Bool selectedsomeseeed;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   finished = FALSE;

   selectedsomeseeed = TRUE;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( (int)conshdlrdata->listall->size() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( SCIPtransformProb(scip) );
      SCIPinfoMessage(scip, NULL, "Applied tranformation to problem.\n");
   }

   commandlen = 0;

   /* Does user want to modify existing or create a new partial decomposition ?*/
   while( (strncmp( command, "modify", commandlen) != 0 && strncmp( command, "create", commandlen) != 0) || commandlen == 0)
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Do you want to modify an existing (\"modify\") or create a new partial decomposition (\"create\")? : \nGCG/toolbox> ",
         &command, &endoffile) );
      commandlen = strlen(command);
   }


   /* case distinction: */
   if( strncmp( command, "modify", commandlen) == 0 )
   {
      /* 1) update list of interesting seeeds */

         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );


         /* 2) while user has not aborted: show current list extract */

         while ( !finished )
         {
            int commandlen2;

            SCIP_CALL( SCIPdialogShowListExtractHeader(scip) );

            SCIP_CALL( SCIPdialogShowListExtract(scip) );

            SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
               "Please choose an existing partial decomposition for modification (type \"choose <id>\" or \"h\" for help) : \nGCG/toolbox> ",
               &command, &endoffile) );

            commandlen2 = strlen(command);

            /* case distinction: */
            if( strncmp( command, "back", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu -= conshdlrdata->selectvisulength;
               if(conshdlrdata->startidvisu < 0 )
                  conshdlrdata->startidvisu = 0;
               continue;
            }
            if( strncmp( command, "next", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu += conshdlrdata->selectvisulength;
               if( conshdlrdata->startidvisu > (int) conshdlrdata->listall->size() - conshdlrdata->selectvisulength )
                  conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
               continue;
            }
            if( strncmp( command, "top", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu = 0;
               continue;
            }
            if( strncmp( command, "end", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
               continue;
            }

            if( strncmp( command, "quit", commandlen2) == 0 )
            {
               finished = TRUE;
               selectedsomeseeed = FALSE;
               continue;
            }


            if( strncmp( command, "choose", commandlen2) == 0 )
            {
               SCIP_RETCODE retcode = SCIPdialogToolboxChoose(scip, dialoghdlr, dialog );
          if (retcode != SCIP_OKAY)
          {
        selectedsomeseeed = FALSE;
        continue;
          }
          else
          {
        finished = TRUE;
        break;
          }
            }


            if( strncmp( command, "abort", commandlen2) == 0 )
            {
               finished = TRUE;
               selectedsomeseeed = FALSE;
               continue;
            }

            if( strncmp( command, "change number displayed", commandlen2) == 0 )
            {
               SCIP_CALL(SCIPdialogModifyNVisualized(scip, dialoghdlr, dialog) );
               continue;
            }

            if( strncmp( command, "help", commandlen2) == 0 )
            {
               SCIP_CALL(SCIPdialogShowHelp(scip) );
               continue;
            }

            if( strncmp( command, "visualize", commandlen2) == 0 )
            {
               SCIP_CALL(SCIPdialogSelectVisualize(scip, dialoghdlr, dialog ) );
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
   } /* finished yes == modify */
   else
   {
      /* create new decomposition */
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ",
         &command, &endoffile) );
      commandlen = strlen(command);

      if( conshdlrdata->curruserseeed != NULL )
         delete conshdlrdata->curruserseeed;

      gcg::Seeedpool* seeedpool;
      SCIP_Bool isfromunpresolved;

      while( (strncmp( command, "presolved", commandlen) != 0 && strncmp( command, "unpresolved", commandlen) != 0) || commandlen == 0)
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
            "Invalid input. Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ",
            &command, &endoffile) );
         commandlen = strlen(command);
      }

      /* case distinction: */
      if( strncmp( command, "presolved", commandlen) == 0 )
      {
         isfromunpresolved = FALSE;
         if (conshdlrdata->seeedpool != NULL )
            seeedpool = conshdlrdata->seeedpool;
         else
         {
            if( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
            {
               SCIPinfoMessage(scip, NULL, "Problem is not presolved yet. Please presolve it first!\n");
               return SCIP_OKAY;
            }

            conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
            seeedpool = conshdlrdata->seeedpool;
         }
      }
      else
      {
         isfromunpresolved = TRUE;
         if ( conshdlrdata->seeedpoolunpresolved == NULL )
            conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE,
               SCIPconshdlrDecompDetectBenders(scip));
          seeedpool = conshdlrdata->seeedpoolunpresolved;

      }
      if( seeedpool == NULL )
      {
         if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED )

         {
            if (conshdlrdata->seeedpool == NULL )
               conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
            seeedpool = conshdlrdata->seeedpool;
         }
         else
         {
            if ( conshdlrdata->seeedpoolunpresolved == NULL)
               conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE,
                  SCIPconshdlrDecompDetectBenders(scip));
            seeedpool = conshdlrdata->seeedpoolunpresolved;
         }

      }
      conshdlrdata->curruserseeed = new gcg::Seeed( scip, SCIPconshdlrDecompGetNextSeeedID(scip), seeedpool );
      conshdlrdata->curruserseeed->setIsFromUnpresolved(isfromunpresolved);
   }

   /* curruserseeed is ready to modify */

   finished = FALSE;
   while ( !finished && selectedsomeseeed )
   {
      int commandlen2;
      SCIP_Bool success;

      SCIP_CALL( SCIPconshdlrDecompShowCurrUserSeeedInfo(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ",
         &command, &endoffile) );

      commandlen2 = strlen(command);

      /* case distinction: */
      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyConss(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyVars(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         if( conshdlrdata->lastuserseeed != NULL)
            delete conshdlrdata->lastuserseeed;
         conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;
         conshdlrdata->curruserseeed->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         gcg::Seeedpool* seeedpool;
         if( !conshdlrdata->curruserseeed->isFromUnpresolved() && conshdlrdata->seeedpool == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = ( conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool);
         assert( seeedpool != NULL );

         conshdlrdata->curruserseeed->sort();
         conshdlrdata->curruserseeed->considerImplicits();
         conshdlrdata->curruserseeed->calcHashvalue();
         assert( conshdlrdata->curruserseeed->checkConsistency() );



         if( conshdlrdata->curruserseeed->isComplete() )
         {
            seeedpool->addSeeedToFinished(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         } else
         {
            seeedpool->addSeeedToIncomplete(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         }
         conshdlrdata->curruserseeed = NULL;
         finished = TRUE;


         continue;
      }

      if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( conshdlrdata->lastuserseeed == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete conshdlrdata->curruserseeed;
            conshdlrdata->curruserseeed = conshdlrdata->lastuserseeed;
            conshdlrdata->lastuserseeed = NULL;
         }
         continue;
      }


      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         SCIP_CALL(SCIPdialogSelectVisualizeCurrentUserSeeed(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "propagate", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPdialogToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
      if( strncmp( command, "finish", commandlen2) == 0 )
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
   SeeedPtr toselect;

   int commandlen;

   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be selected:\n", conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   toselect = conshdlrdata->listall->at(idtovisu);

   toselect->setSelected(!toselect->isSelected() );

   if( !toselect->isSelected() )
   {
      conshdlrdata->selected->erase(  find( conshdlrdata->selected->begin(), conshdlrdata->selected->end(), idtovisu) );
   }
   else
   {
      std::cout << "is selected!" << toselect->isSelected() <<std::endl;
      conshdlrdata->selected->push_back(idtovisu);
      assert(toselect->isSelected());
   }

   conshdlrdata->selectedexists = (conshdlrdata->selected->size() > 0);

   return SCIP_OKAY;
}



