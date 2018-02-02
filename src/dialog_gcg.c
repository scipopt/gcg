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

/**@file   dialog_gcg.c
 * @brief  GCG user interface dialog
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Gerald Gamrath
 * @author Martin Bergner
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/pub_dialog.h"
#include "scip/type_dialog.h"
#include "scip/dialog_default.h"

#include "gcg.h"

#include "dialog_gcg.h"
#include "relax_gcg.h"
#include "pricer_gcg.h"
#include "cons_decomp.h"
#include "pub_gcgheur.h"
#include "stat.h"
#include "reader_dec.h"

/* display the reader information */
static
void displayReaders(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             reader,             /**< display reader which can read */
   SCIP_Bool             writer              /**< display reader which can write */
   )
{
   SCIP_READER** readers;
   int nreaders;
   int r;

   assert( scip != NULL );

   readers = SCIPgetReaders(scip);
   nreaders = SCIPgetNReaders(scip);

   /* display list of readers */
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, " file reader          extension  description\n");
   SCIPdialogMessage(scip, NULL, " -----------          ---------  -----------\n");
   for( r = 0; r < nreaders; ++r )
   {
      if( (reader && SCIPreaderCanRead(readers[r])) || (writer && SCIPreaderCanWrite(readers[r])) )
      {
         SCIPdialogMessage(scip, NULL, " %-20s ", SCIPreaderGetName(readers[r]));
         if( strlen(SCIPreaderGetName(readers[r])) > 20 )
            SCIPdialogMessage(scip, NULL, "\n %20s ", "-->");
         SCIPdialogMessage(scip, NULL, "%9s  ", SCIPreaderGetExtension(readers[r]));
         SCIPdialogMessage(scip, NULL, "%s", SCIPreaderGetDesc(readers[r]));
         SCIPdialogMessage(scip, NULL, "\n");
      }
   }
   SCIPdialogMessage(scip, NULL, "\n");
}


/** writes out all decompositions currently known to cons_decomp */
static
SCIP_RETCODE writeAllDecompositions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog menu */
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG**         nextdialog          /**< pointer to store next dialog to execute */
   )
{

   char* filename;
   char* dirname;
   SCIP_Bool endoffile;

   if( SCIPconshdlrDecompGetNDecdecomps(scip) == 0 )
   {
      SCIPdialogMessage(scip, NULL, "No decomposition to write, please read or detect one first.\n");
      SCIPdialoghdlrClearBuffer(dialoghdlr);
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter directory and/or extension: ", &dirname, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if( SCIPdialoghdlrIsBufferEmpty(dialoghdlr) )
   {
      filename = dirname;
      dirname = NULL;
   }
   else
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter extension: ", &filename, &endoffile) );
   }

   if( filename[0] != '\0' )
   {
      char* extension;
      extension = filename;
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, extension, TRUE) );

      do
      {
         SCIP_RETCODE retcode = DECwriteAllDecomps(scip, dirname, extension);

         if( retcode == SCIP_FILECREATEERROR )
         {
            SCIPdialogMessage(scip, NULL, "error creating files\n");
            SCIPdialoghdlrClearBuffer(dialoghdlr);
            break;
         }
         else if( retcode == SCIP_WRITEERROR )
         {
            SCIPdialogMessage(scip, NULL, "error writing files\n");
            SCIPdialoghdlrClearBuffer(dialoghdlr);
            break;
         }
         else if( retcode == SCIP_PLUGINNOTFOUND )
         {
            /* ask user once for a suitable reader */
            if( extension == NULL )
            {
               SCIPdialogMessage(scip, NULL, "no reader for requested output format\n");

               SCIPdialogMessage(scip, NULL, "following readers are avaliable for writing:\n");
               displayReaders(scip, FALSE, TRUE);

               SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
                     "select a suitable reader by extension (or return): ", &extension, &endoffile) );

               if( extension[0] == '\0' )
                  break;
            }
            else
            {
               SCIPdialogMessage(scip, NULL, "no reader for output in <%s> format\n", extension);
               extension = NULL;
            }
         }
         else
         {
            /* check for unexpected errors */
            SCIP_CALL( retcode );

            /* print result message if writing was successful */
            SCIPdialogMessage(scip, NULL, "written all decompositions %s\n", extension);
            break;
         }
      }
      while (extension != NULL );
   }

   return SCIP_OKAY;
}

/** dialog execution method for the display statistics command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayStatistics)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIP_CALL( GCGprintStatistics(scip, NULL) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display decomposition command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayDecomposition)
{  /*lint --e{715}*/
   DEC_DECOMP* decomp;
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   decomp = DECgetBestDecomp(scip);
   if( decomp != NULL )
   {
      SCIP_CALL( GCGwriteDecomp(scip, NULL, decomp) );
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}


/** dialog execution method for the display additionalstatistics command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayAdditionalStatistics)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
   {
      if( SCIPgetStage(GCGgetMasterprob(scip)) < SCIP_STAGE_PRESOLVED )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), NULL, "No Dantzig-Wolfe reformulation applied. No decomposition statistics available.\n");
         *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
         return SCIP_OKAY;
      }

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), NULL, "\nAdditional statistics:\n");
      if( DECdecompGetType(DECgetBestDecomp(scip)) == DEC_DECTYPE_DIAGONAL )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(GCGgetMasterprob(scip)), NULL, "\n");
         SCIP_CALL( GCGwriteDecompositionData(scip) );

      }
      else
      {
         GCGpricerPrintStatistics(GCGgetMasterprob(scip), NULL);
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(GCGgetMasterprob(scip)), NULL, "\n");
         SCIP_CALL( GCGwriteDecompositionData(scip) );
         SCIP_CALL( GCGwriteVarCreationDetails(GCGgetMasterprob(scip)) );
      }
   }
   else
   {
      SCIPdialogMessage(scip, NULL, "Problem needs to solved first for additional statistics");
   }
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}


/** dialog execution method for the display detectors command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplayDetectors)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   DECprintListOfDetectors(scip);
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display solvers command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecDisplaySolvers)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   GCGpricerPrintListOfSolvers(GCGgetMasterprob(scip));
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the master command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetMaster)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(GCGgetMasterprob(scip)) != SCIP_STAGE_INIT )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "switching to the master problem shell is only possible before the solving process is started\n");

      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

      return SCIP_OKAY;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "switching to the master problem...\n");
   SCIP_CALL( SCIPstartInteraction(GCGgetMasterprob(scip)) );
   SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "back in the original problem...\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set loadmaster command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetLoadmaster)
{  /*lint --e{715}*/
   SCIP* masterprob;
   char* filename;
   SCIP_Bool endoffile;

   masterprob = GCGgetMasterprob(scip);
   assert(masterprob != NULL);

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }

   if( filename[0] != '\0' )
   {
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      if( SCIPfileExists(filename) )
      {
         SCIP_CALL( SCIPreadParams(masterprob, filename) );
         SCIPdialogMessage(scip, NULL, "loaded master parameter file <%s>\n", filename);
      }
      else
      {
         SCIPdialogMessage(scip, NULL, "file <%s> not found\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the detect command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecDetect)
{  /*lint --e{715}*/
   SCIP_RESULT result;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Starting detection\n");
   if( SCIPgetStage(scip) > SCIP_STAGE_INIT )
   {
      SCIP_CALL( DECdetectStructure(scip, &result) );
      if( result == SCIP_SUCCESS )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Detection was successful.\n");
      else
            SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Detection was not successful.\n");
   }
   else
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the optimize command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecOptimize)
{  /*lint --e{715}*/
   SCIP_RESULT result;
   int presolrounds;
   presolrounds = -1;

   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   SCIPdialogMessage(scip, NULL, "\n");
   switch( SCIPgetStage(scip) )
   {
   case SCIP_STAGE_INIT:
      SCIPdialogMessage(scip, NULL, "No problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      if( DEChasDetectionRun(scip) || (DECgetBestDecomp(scip) != NULL) )
      {
         SCIP_CALL( SCIPgetIntParam(scip, "presolving/maxrounds", &presolrounds) );
         SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
      }

      SCIP_CALL( SCIPpresolve(scip) ); /*lint -fallthrough*/

   case SCIP_STAGE_PRESOLVED:
      if( !DEChasDetectionRun(scip) )
      {
         SCIP_CALL( DECdetectStructure(scip, &result) );
         if( result == SCIP_DIDNOTFIND )
         {
            assert(DECgetBestDecomp(scip) == NULL && DEChasDetectionRun(scip));
            SCIPdialogMessage(scip, NULL, "No decomposition exists or could be detected. You need to specify one.\n");
            break;
         }
      }
      else if( DECgetBestDecomp(scip) == NULL )
      {
         assert(DECgetBestDecomp(scip) == NULL && DEChasDetectionRun(scip));
         SCIPdialogMessage(scip, NULL, "No decomposition exists or could be detected. You need to specify one.\n");
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPsolve(scip) );
      break;

   case SCIP_STAGE_SOLVED:
      SCIPdialogMessage(scip, NULL, "Problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_EXITSOLVE:
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_FREE:
   default:
      SCIPerrorMessage("Invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   if( presolrounds != -1 )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", presolrounds) );
   }

   return SCIP_OKAY;
}

/** dialog execution method for writing all known decompositions */
static
SCIP_DECL_DIALOGEXEC(GCGdialogExecWriteAllDecompositions)
{
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM )
   {
      SCIP_CALL( writeAllDecompositions(scip, dialog, dialoghdlr, nextdialog) );
   }
   else
      SCIPdialogMessage(scip, NULL, "no problem available\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for writing problem statistics */
static
SCIP_DECL_DIALOGEXEC(GCGdialogExecWriteStatistics)
{
   char* filename;
   SCIP_Bool endoffile;

   SCIPdialogMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      FILE* file;

      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename, TRUE) );

      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPdialogMessage(scip, NULL, "error creating file <%s>\n", filename);
         SCIPprintSysError(filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
      else
      {
         SCIP_RETCODE retcode;
         retcode = GCGprintStatistics(scip, file);
         if( retcode != SCIP_OKAY )
         {
            fclose(file);
            SCIP_CALL( retcode );
         }
         else
         {
            SCIPdialogMessage(scip, NULL, "written statistics to file <%s>\n", filename);
            fclose(file);
         }
      }
   }

   SCIPdialogMessage(scip, NULL, "\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}
/** dialog execution method for the set detectors aggressive command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsAggressive)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( GCGsetDetection(scip, SCIP_PARAMSETTING_AGGRESSIVE, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set detectors default command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsDefault)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( GCGsetDetection(scip, SCIP_PARAMSETTING_DEFAULT, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set detectors off command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsOff)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( GCGsetDetection(scip, SCIP_PARAMSETTING_OFF, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set detectors fast command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDetectorsFast)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( GCGsetDetection(scip, SCIP_PARAMSETTING_FAST, FALSE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set heuristics aggressive command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetHeuristicsAggressive)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, FALSE) );
   SCIP_CALL( GCGsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE) );

   return SCIP_OKAY;
}

/** dialog execution method for the set heuristics off command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetHeuristicsOff)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, FALSE) );
   SCIP_CALL( GCGsetHeuristics(scip, SCIP_PARAMSETTING_OFF) );

   return SCIP_OKAY;
}

/** dialog execution method for the set heuristics fast command */
SCIP_DECL_DIALOGEXEC(GCGdialogExecSetHeuristicsFast)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL, FALSE) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_FAST, FALSE) );
   SCIP_CALL( GCGsetHeuristics(scip, SCIP_PARAMSETTING_FAST) );

   return SCIP_OKAY;
}

/** creates a root dialog */
SCIP_RETCODE GCGcreateRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   )
{
   SCIP_CALL( SCIPincludeDialog(scip, root, NULL, SCIPdialogExecMenuLazy, NULL, NULL,
         "GCG", "GCG's main menu", TRUE, NULL) );

   SCIP_CALL( SCIPsetRootDialog(scip, *root) );
   SCIP_CALL( SCIPreleaseDialog(scip, root) );
   *root = SCIPgetRootDialog(scip);

   return SCIP_OKAY;
}

/** create a "emphasis" sub menu */
static
SCIP_RETCODE createEmphasisSubmenu(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          root,               /**< the menu to add the empty sub menu */
   SCIP_DIALOG**         submenu             /**< pointer to store the created emphasis sub menu */
   )
{
   if( !SCIPdialogHasEntry(root, "emphasis") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "emphasis", "predefined parameter settings", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, *submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, submenu) );
   }
   else if( SCIPdialogFindEntry(root, "emphasis", submenu) != 1 )
   {
      SCIPerrorMessage("emphasis sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   assert(*submenu != NULL);

   return SCIP_OKAY;
}

/** includes or updates the GCG dialog menus in SCIP */
SCIP_RETCODE SCIPincludeDialogGcg(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DIALOG* root;
   SCIP_DIALOG* submenu;
   SCIP_DIALOG* dialog;
   SCIP_DIALOG* setmenu;
   SCIP_DIALOG* emphasismenu;

   /* root menu */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIP_CALL( GCGcreateRootDialog(scip, &root) );
   }

   /* display */
   if( !SCIPdialogHasEntry(root, "display") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu, NULL, SCIPdialogExecMenu, NULL, NULL,
            "display", "display information", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "display", &submenu) != 1 )
   {
      SCIPerrorMessage("display sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* display statistics */
   if( !SCIPdialogHasEntry(submenu, "statistics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecDisplayStatistics, NULL, NULL,
            "statistics", "display problem and optimization statistics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }
   /* display decomposition */
   if( !SCIPdialogHasEntry(submenu, "decomposition") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecDisplayDecomposition, NULL, NULL,
            "decomposition", "display decomposition", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display additionalstatistics */
   if( !SCIPdialogHasEntry(submenu, "additionalstatistics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecDisplayAdditionalStatistics, NULL, NULL,
            "additionalstatistics", "display additional solving statistics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display detectors */
   if( !SCIPdialogHasEntry(submenu, "detectors") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecDisplayDetectors, NULL, NULL,
            "detectors", "display available detectors", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* display solvers */
   if( !SCIPdialogHasEntry(submenu, "solvers") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecDisplaySolvers, NULL, NULL,
            "solvers", "display available pricing problem solvers", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* master */
   if( !SCIPdialogHasEntry(root, "master") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecSetMaster, NULL, NULL,
            "master", "switch to the interactive shell of the master problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* optimize */
   if( !SCIPdialogHasEntry(root, "optimize") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            GCGdialogExecOptimize, NULL, NULL,
            "optimize", "solve the problem", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }


   /* detect */
   if( !SCIPdialogHasEntry(root, "detect") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecDetect, NULL, NULL,
            "detect", "Detect structure", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* quit */
   if( !SCIPdialogHasEntry(root, "quit") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecQuit, NULL, NULL,
            "quit", "leave GCG", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set */
   if( !SCIPdialogHasEntry(root, "set") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL, SCIPdialogExecMenu, NULL, NULL,
            "set", "load/save/change parameters", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "set", &setmenu) != 1 )
   {
      SCIPerrorMessage("set sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* set loadmaster */
   if( !SCIPdialogHasEntry(setmenu, "loadmaster") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            GCGdialogExecSetLoadmaster, NULL, NULL,
            "loadmaster", "load parameter settings for master problem from a file", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set detectors */
   if( !SCIPdialogHasEntry(setmenu, "detectors") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
            NULL,
            SCIPdialogExecMenu, NULL, NULL,
            "detectors", "change parameters for detectors", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "detectors", &submenu) != 1 )
   {
      SCIPerrorMessage("detectors sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create set detectors emphasis */
   SCIP_CALL( createEmphasisSubmenu(scip, submenu, &emphasismenu) );
   assert(emphasismenu != NULL);

   /* set detectors emphasis aggressive */
   if( !SCIPdialogHasEntry(emphasismenu, "aggressive") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetDetectorsAggressive, NULL, NULL,
            "aggressive", "sets detectors <aggressive>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set detectors emphasis default */
   if( !SCIPdialogHasEntry(emphasismenu, "default") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
         NULL, SCIPdialogExecSetDetectorsDefault, NULL, NULL,
         "default", "sets detectors <default>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set detectors emphasis fast */
   if( !SCIPdialogHasEntry(emphasismenu, "fast") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetDetectorsFast, NULL, NULL,
            "fast", "sets detectors <fast>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set detectors emphasis off */
   if( !SCIPdialogHasEntry(emphasismenu, "off") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL, SCIPdialogExecSetDetectorsOff, NULL, NULL,
            "off", "turns <off> all detectors", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set heuristics */
   if( !SCIPdialogHasEntry(setmenu, "heuristics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu,
         NULL,
         SCIPdialogExecMenu, NULL, NULL,
         "heuristics", "change parameters for primal heuristics", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, setmenu, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "heuristics", &submenu) != 1 )
   {
      SCIPerrorMessage("heuristics sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create set heuristics emphasis */
   SCIP_CALL( createEmphasisSubmenu(scip, submenu, &emphasismenu) );
   assert(emphasismenu != NULL);

   /* set heuristics emphasis aggressive */
   if( !SCIPdialogHasEntry(emphasismenu, "aggressive") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
         NULL, GCGdialogExecSetHeuristicsAggressive, NULL, NULL,
         "aggressive", "sets heuristics <aggressive>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set heuristics emphasis fast */
   if( !SCIPdialogHasEntry(emphasismenu, "fast") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
         NULL, GCGdialogExecSetHeuristicsFast, NULL, NULL,
         "fast", "sets heuristics <fast>", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set heuristics emphasis off */
   if( !SCIPdialogHasEntry(emphasismenu, "off") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
         NULL, GCGdialogExecSetHeuristicsOff, NULL, NULL,
         "off", "turns <off> all heuristics", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, emphasismenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write */
   if( !SCIPdialogHasEntry(root, "write") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &submenu, NULL, SCIPdialogExecMenu, NULL, NULL,
            "write", "write information to file", TRUE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, root, submenu) );
      SCIP_CALL( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "write", &submenu) != 1 )
   {
      SCIPerrorMessage("write sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* write alldecompositions */
   if( !SCIPdialogHasEntry(submenu, "alldecompositions") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecWriteAllDecompositions, NULL, NULL,
            "alldecompositions",
            "write all known decompostions to file (format is given by file extension, e.g., {dec,blk,ref})",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   /* write statistics */
   if( !SCIPdialogHasEntry(submenu, "statistics") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog, NULL, GCGdialogExecWriteStatistics, NULL, NULL,
            "statistics",
            "write statistics to file",
            FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }


   return SCIP_OKAY;
}
