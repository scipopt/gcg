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
#include <stdlib.h>

#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <sstream>

#include "class_seeed.h"
#include "cons_decomp.h"
#include "wrapper_seeed.h"

/* column headers */
#define DEFAULT_COLUMN_MIN_WIDTH  4 /**< min width of a column in the menu table */
#define DEFAULT_COLUMN_MAX_WIDTH 10 /**< max width of a column (also determines max width of column header abbreviation) */
#define DEFAULT_COLUMNS "nr id nbloc nmacon nlivar nmavar nstlva score history pre nopcon nopvar sel" /**< default column headers */
#define DEFAULT_SORT_HEADER "score"

#define DEFAULT_MENULENGTH 10 /**< initial number of entries in menu */

namespace gcg
{

/** rettype is used to store the return type of a callback function */
enum RETTYPE{
   UNKNOWN,    /**< dummy default to catch errors */
   INTEGER,    /**< integer */
   FLOAT,      /**< float */
   BOOLEAN,    /**< Boolean */
   STRING      /**< char* */
};

/** callback for getters of column infos */
typedef void (*callback)(SCIP*, int);

/** storage for column information */
struct Columninfo{
   std::string header;  /**< table header for the column */
   std::string desc;    /**< description of the column entries used in the menu help */
   callback getter;     /**< callback to get the infos displayed in the column for a given scip & seeed id */
   RETTYPE type;        /**< return type of the getter */

    Columninfo(std::string nheader, std::string ndesc, callback ngetter, RETTYPE ntype) :
    header(nheader), desc(ndesc), getter(ngetter), type(ntype) {}
    /**< constructor for easier initialization */
};

/** gets the seeed structure from a given id (local help function)
 *
 * @todo remove this help function once the seeed structure is depricated
 * @returns seeed for given id
*/
static
Seeed* getSeeed(
   SCIP* scip,    /**< SCIP data structure */
   int id         /**< id of seeed */
   )
{
   Seeed_Wrapper sw;
   GCGgetSeeedFromID(scip, &id, &sw);
   assert( sw.seeed != NULL );
   return sw.seeed;
}


/** @brief local sorting function for seeed id vectors
 *
 * avoids redundant sorting calls,
 * sorts by score in given order
 */
static
void sortSeeedList(
   SCIP* scip,                      /**< SCIP data structure */
   std::vector<int>* idlist,        /**< current list of seeed ids */
   std::string header,              /**< header of column to sort by */
   std::vector<Columninfo*> columns,/**< column infos */
   bool asc                         /**< whether to sort ascending or descending */
   )
{
   /* find the column infos for the given header */
   for(auto column : columns)
   {
      if( column->header.find( header ) == 0 )
      {
         /* sort the id list according to given order using the callback getter of the column */
         if(column->type == INTEGER)
         {
            /* the callback has to be parsed to expect an int output */
            if(asc)
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ((*( (int(*)(SCIP*, int)) column->getter))(scip, a) < (*( (int(*)(SCIP*, int)) column->getter))(scip, b)); });
            else
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ((*( (int(*)(SCIP*, int)) column->getter))(scip, a) > (*( (int(*)(SCIP*, int)) column->getter))(scip, b)); });
         }
         else if(column->type == FLOAT)
         {
            /* the callback has to be parsed to expect a float output */
            if(asc)
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ((*( (float(*)(SCIP*, int)) column->getter))(scip, a) < (*( (float(*)(SCIP*, int)) column->getter))(scip, b)); });
            else
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ((*( (float(*)(SCIP*, int)) column->getter))(scip, a) > (*( (float(*)(SCIP*, int)) column->getter))(scip, b)); });
         }
         else if(column->type == BOOLEAN)
         {
            /* the callback has to be parsed to expect a SCIP_Bool output */
            if(asc)
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ((*( (SCIP_Bool(*)(SCIP*, int)) column->getter))(scip, a) < (*( (SCIP_Bool(*)(SCIP*, int)) column->getter))(scip, b)); });
            else
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ((*( (SCIP_Bool(*)(SCIP*, int)) column->getter))(scip, a) > (*( (SCIP_Bool(*)(SCIP*, int)) column->getter))(scip, b)); });
         }
         else if(column->type == STRING)
         {
            /* the callback has to be parsed to expect a char* output, the comparison requires another cast into string */
            if(asc)
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ( (std::string) ((*( (char*(*)(SCIP*, int)) column->getter))(scip, a)) < (std::string) ((*( (char*(*)(SCIP*, int)) column->getter))(scip, b))); });
            else
               std::sort(idlist->begin(), idlist->end(), [&](const int a, const int b) {return ( (std::string) ((*( (char*(*)(SCIP*, int)) column->getter))(scip, a)) > (std::string) ((*( (char*(*)(SCIP*, int)) column->getter))(scip, b))); });
         }
         break;
     }
  }
}


/** modifies menulength according to input and updates menu accordingly
 * @returns SCIP return code */
static
SCIP_RETCODE GCGdialogSetNEntires(
   SCIP* scip,                   /**< SCIP data structure */
   SCIP_DIALOGHDLR* dialoghdlr,  /**< dialog handler for user input management */
   SCIP_DIALOG* dialog,          /**< dialog for user input management */
   int listlength,               /**< length of seeed id list */
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

   /* check for invalid input */
   if( commandlen == 0 || newlength < 1 )
   {
      SCIPdialogMessage( scip, NULL, "The input was not a valid number." );
      return SCIP_OKAY;
   }

   /* set new length (max listlength) */
   if( newlength < listlength )
      *menulength = newlength;
   else
      *menulength = listlength;

   return SCIP_OKAY;
}


/** sets the used score according to user input
 *
 * @returns SCIP return code
 */
static
SCIP_RETCODE GCGdialogChangeScore(
   SCIP* scip,                         /**< SCIP data structure */
   SCIP_DIALOGHDLR* dialoghdlr,        /**< dialog handler for user input management */
   SCIP_DIALOG* dialog                 /**< dialog for user input management */
   )
{
   char* getscore;
   SCIP_Bool endoffile;
   int commandlen;

   SCIPdialogMessage(scip, NULL, "\nPlease specify the new score:\n");
   SCIPdialogMessage(scip, NULL, "0: max white, \n1: border area, \n2: classic, \n3: max foreseeing white, \n4: ppc-max-white, \n");
   SCIPdialogMessage(scip, NULL, "5: max foreseeing white with aggregation info, \n6: ppc-max-white with aggregation info, \n7: experimental benders score\n");
   SCIPdialogMessage(scip, NULL, "8: strong decomposition score\n");
   SCIPdialogMessage(scip, NULL, "Note: Sets the detection/scoretype parameter to the given score.\n");

   /* get input */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &getscore, &endoffile) );
   commandlen = strlen(getscore);
   if( commandlen != 0 )
   {
      /* convert to int (if value is invalid, this results in 0) */
      int scorenr = atoi(getscore);

      /* check if the value is in valid range */
      if(scorenr >= 0 && scorenr <= 8)
      {
         /* set score */
         SCIPsetIntParam(scip, "detection/scoretype", scorenr);
         SCIPconshdlrDecompSetScoretype(scip, static_cast<SCORETYPE>(scorenr));
         SCIPdialogMessage(scip, NULL, "Score set to %d.\n", scorenr);
      }
   }

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
 * Update seeed list in case it changed since the last command
 * and show the table of seeeds.
 * @returns SCIP status
 */
static
SCIP_RETCODE GCGdialogShowMenu(
   SCIP* scip,                            /**< SCIP data structure */
   std::vector<Columninfo*> columns,      /**< list of column headers/ info sources */
   int* nseeeds,                          /**< max number of seeeds */
   const int startindex,                  /**< index (in seeed list) of uppermost seeed in extract */
   int menulength,                        /**< number of menu entries */
   std::vector<int>* idlist,              /**< current list of seeed ids */
   bool sortasc,                          /**< true iff sorting should be ascending */
   std::string sortby                     /**< table header of column to sort by */
   )
{
   assert(scip != NULL);

   /* update seeed list in case it changed (in which case the amount of seeeds should have changed)*/
   if(*nseeeds < SCIPconshdlrDecompGetNSeeeds(scip))
   {
      *nseeeds = SCIPconshdlrDecompGetNSeeeds(scip);
      int* idarray;
      int listlength;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idarray, *nseeeds) );
      SCIPconshdlrDecompGetSeeedLeafList(scip, &idarray, &listlength);

      /* reset idlist to the new idarray */
      idlist->clear();
      for(int i = 0; i < listlength; i++)
      {
         idlist->push_back(idarray[i]);
      }

      /* free idarray */
      SCIPfreeBlockMemoryArray(scip, &idarray, *nseeeds);
   }

   /* sort seeed ids by score, descending (in case score was changed or id list was updated)*/
   sortSeeedList(scip, idlist, sortby, columns, sortasc);

   /* count corresponding seeeds for overview statistics */
   int ndetectedpresolved = 0;
   int ndetectedunpresolved = 0;

   for(int i = 0; i < (int) idlist->size(); ++i)
   {
      Seeed* seeed = getSeeed(scip, idlist->at(i));
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

   /* add each column header */
   for(auto column : columns)
   {
      /* "score" is a wildcard for the current score, relace it with actual scoretype */
      std::string header = column->header;
      std::string newheader;
      if(header != "score")
         newheader = header;
      else
         newheader = SCIPconshdlrDecompGetScoretypeShortName(scip, SCIPconshdlrDecompGetScoretype(scip));

      /* make sure the header name is unique and add a length for header */
      assert(columnlength.find(header) == columnlength.end());
      columnlength.insert(std::pair<std::string,int>(header, 0));
      /* if header is smaller than min column width, add spaces to header first */
      if(newheader.size() < DEFAULT_COLUMN_MIN_WIDTH)
      {
         for(int i = 0; i < (DEFAULT_COLUMN_MIN_WIDTH - (int) newheader.size()); i++)
         {
            headerline += " ";
            borderline += "-";
            columnlength.at(header)++;
         }
      }
      /* add header to headerline and add #chars of header as '-' to borderline*/
      headerline += newheader;
      for(int i = 0; i < (int) newheader.size(); i++)
         borderline += "-";
      columnlength.at(header) += (int) newheader.size();
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

   /* go through all seeeds that should currently be displayed,
    * so from startindex on menulength many entries if there are that much left in the list */
   for(int i = startindex; i < startindex + menulength && i < (int) idlist->size(); ++i)
   {
      /* get current seeed id */
      int seeedid = idlist->at(i);

      /* each line starts with a space */
      SCIPdialogMessage(scip, NULL, " ");

      /* go through the columns and write the entry for each one */
      for(auto column : columns)
      {
         std::string towrite;
         /* get the type of the callback function that is the getter for this column */
         RETTYPE type = column->type;
         /* get the header of this column */
         std::string header = column->header;

         /* call the getter of the column depending on the given return type */
         if(type == UNKNOWN)
         {
            /* "nr" and "id" are special cases and should be the only ones where the type is UNKNOWN */
            if(header == "nr")
               /* "nr" is the current position of the seeed in the menu list */
               towrite = std::to_string(i);
            else if(header == "id")
               /* "id" is the seeed's id */
               towrite = std::to_string(seeedid);
         }
         else if(type == INTEGER)
            /* convert the callback function to int rettype, call it on (scip, seeedid) and convert the result to a string */
            towrite = std::to_string( (*( (int(*)(SCIP*, int)) column->getter ))(scip, seeedid) );
         else if(type == FLOAT)
         {
            /* convert the callback function to float rettype, call it on (scip, seeedid) */
            float number = (*( (float(*)(SCIP*, int)) column->getter ))(scip, seeedid);
            /* convert the result to a string and set the number of digits (i.e.letters) to the width of the column */
            towrite = std::to_string(number).substr(0, columnlength.at(header));
         }
         else if(type == BOOLEAN)
            /* convert the callback function to SCIP_Bool rettype, call it on (scip, seeedid) and check the result */
            towrite = ( (*( (SCIP_Bool(*)(SCIP*, int)) column->getter ))(scip, seeedid) ) ? "yes" : "no";
         else if(type == STRING)
            /* convert the callback function to char* rettype and call it on (scip, seeedid) */
            towrite = (*( (char* (*)(SCIP*, int)) column->getter ))(scip, seeedid);
         else
            towrite = " ";

         /* write spaces to fill out the columnwidth until towrite */
         outputCharXTimes(scip, ' ', (columnlength.at(header) - (int) towrite.size()));
         /* write actual value of the column +1 space for border */
         SCIPdialogMessage(scip, NULL, "%s ", towrite.c_str());
      }

      /* continue to next line */
      SCIPdialogMessage(scip, NULL, "\n");
   }

   /* at the end of the table add a line */
   outputCharXTimes(scip, '=', linelength);

   return SCIP_OKAY;
}


/** Shows information about the explore screen and its abbreviations
 *
 * @returns SCIP status */
static
SCIP_RETCODE GCGdialogShowLegend(
   SCIP* scip,                         /**< SCIP data structure */
   std::vector<Columninfo*> columns    /**< list of column headers/ info sources */
   )
{
   assert(scip != NULL);
   DEC_DETECTOR** detectors;

   /* print header for detector list */
   SCIPdialogMessage(scip, NULL, "List of included detectors for decompositions histories: \n");

   SCIPdialogMessage(scip, NULL, "\n%30s    %4s\n", "detector" , "char");
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "--------" , "----");

   /* get and print char of each detector */
   detectors = SCIPconshdlrDecompGetDetectors(scip);

   for( int det = 0; det < SCIPconshdlrDecompGetNDetectors(scip); ++det )
   {
      DEC_DETECTOR* detector;
      detector = detectors[det];

      SCIPdialogMessage(scip, NULL, "%30s    %4c\n", DECdetectorGetName(detector), DECdetectorGetChar(detector));
   }
   /* print usergiven as part of detector chars */
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "given by user" , "U");
   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "=================================================================================================== \n");

   SCIPdialogMessage(scip, NULL, "\n" );

   /* print header of abbreviation table */
   SCIPdialogMessage(scip, NULL, "List of abbreviations of decomposition table \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "abbreviation", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "------------", "-----------");

   /* add legend entry for each header abbreviation */
   for(auto column : columns)
   {
      /* get table header */
      std::string header = column->header;

      /* print the header with the description */
      if(header != "score")
      {
         SCIPdialogMessage(scip, NULL, "%30s     %s\n", header.c_str(), column->desc.c_str());
      }
      /* if the header is "score" replace with shortname of the current score */
      else
      {
         SCIPdialogMessage(scip, NULL, "%30s     %s\n", SCIPconshdlrDecompGetScoretypeShortName(scip, SCIPconshdlrDecompGetScoretype(scip)),
            SCIPconshdlrDecompGetScoretypeDescription(scip, SCIPconshdlrDecompGetScoretype(scip)) );
      }

   }
   SCIPdialogMessage(scip, NULL, "\n=================================================================================================== \n");

   return SCIP_OKAY;
}

/** @brief Shows help section of explore menu
 *
 * Outputs al ist of commands and a description of their function
 * @returns SCIP status */
static
SCIP_RETCODE GCGdialogShowHelp(
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
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "help", "displays this help");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "legend", "displays the legend for table header and history abbreviations");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "select", "selects/unselects decomposition with given nr");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "previous", "displays the preceding decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "next", "displays the subsequent decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "top", "displays the first decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "end", "displays the last decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "entries", "modifies the number of decompositions to display per page");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "visualizes the specified decomposition (requires gnuplot)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "inspect", "displays detailed information for the specified decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "score", "sets the score by which the quality of decompositions is evaluated");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "sort", "sets the column by which the decompositions are sorted (default: by score)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "ascending", "sort decompositions in ascending (true) or descending (false) order");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "return to main menu");

   SCIPdialogMessage(scip, NULL, "\n=================================================================================================== \n");

   return SCIP_OKAY;
}


/** Shows a visualization of the seeed specified by the user via the dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE GCGdialogSelectVisualize(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   std::vector<int>        idlist      /**< current list of seeed ids */
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
   if( commandlen == 0 || idtovisu < 0 || idtovisu >= (int) idlist.size() )
   {
      SCIPdialogMessage( scip, NULL, "This nr is out of range." );
      return SCIP_OKAY;
   }

   /* get and show seeed */
   Seeed* seeed = getSeeed(scip, idlist.at(idtovisu));
   seeed->showVisualisation();

   return SCIP_OKAY;
}


/**
 * Displays information about a seeed that is chosen by the user in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE GCGdialogInspectSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   std::vector<int>        idlist      /**< current list of seeed ids */
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

   if(idtoinspect < 0 || idtoinspect >= (int) idlist.size()){
      SCIPdialogMessage( scip, NULL, "This nr is out of range." );
      return SCIP_OKAY;
   }

   /* check whether ID is in valid range */
   Seeed* seeed = getSeeed(scip, idlist.at(idtoinspect));

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
   seeed->displayInfo( detaillevel );

   return SCIP_OKAY;
}


/** Lets the user select decompositions from the explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE GCGdialogSelect(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   std::vector<int>        idlist      /**< current list of seeed ids */
   )
{
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;

   int commandlen;

   assert(scip != NULL);

   /* get input */
   SCIPdialogMessage(scip, NULL, "\nPlease specify the nr of the decomposition to be selected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   /* check if the input is a valid number */
   if( commandlen == 0 || idtovisu < 0 || idtovisu >= (int) idlist.size() )
   {
      SCIPdialogMessage( scip, NULL, "This nr is out of range, nothing was selected." );
      return SCIP_OKAY;
   }

   /* get seeed from id*/
   Seeed* seeed = getSeeed(scip, idlist.at(idtovisu));

   /* reverse selection (deselects if seeed was previously selected) */
   seeed->setSelected(!seeed->isSelected() );

   return SCIP_OKAY;
}

/** Set whether order in menu should be ascending/descending
 *
 * @returns SCIP return code */
static
SCIP_RETCODE GCGdialogSortAsc(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   bool*                   asc         /**< true iff sorting should be ascending */
   )
{
   char* ascen;
   SCIP_Bool endoffile;
   int commandlen;

   assert(scip != NULL);

   /* get input */
   SCIPdialogMessage(scip, NULL, "\nPlease enter \"true\"/\"1\" for ascending or \"false\"/\"0\" for descending order:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ascen, &endoffile) );
   commandlen = strlen(ascen);

   std::string input = ascen;
   if( commandlen != 0)
   {
      /* check if input if is true, set ascending */
      if(input == "true" || input == "1")
         *asc = true;
      /* check if input if is false, set descending */
      else if(input == "false" || input == "0")
         *asc = false;
      /* all other inputs are considered invalid and do nothing */
   }

   return SCIP_OKAY;
}


/** Checks whether the given header is valid
 * @returns true iff header is valid
 */
static
bool isHeader(
   std::string header,                 /**< header to check */
   std::vector<Columninfo*> columns    /**< list of column headers/ info sources */  
   )
{
   /* check if the given header is a (prefix of a) registered table header */
   for(auto column : columns)
   {
      if(column->header.find( header ) == 0 )
         return true;
   }
   /* else return false */
   return false;
}


/** Set whether order in menu should be ascending/descending
 *
 * @returns SCIP return code */
static
SCIP_RETCODE GCGdialogSortBy(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   std::vector<Columninfo*> columns,   /**< list of column headers/ info sources */
   std::string*            sortby      /**< table header, identifies by which column to sort by */
   )
{
   char* newsort;
   SCIP_Bool endoffile;
   int commandlen;

   assert(scip != NULL);

   /* get input */
   SCIPdialogMessage(scip, NULL, "\nPlease enter the table header of the column by which you would like to sort:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &newsort, &endoffile) );
   commandlen = strlen(newsort);

   /* if the input is a valid table header, change sortby */
   std::string input = newsort;
   if( commandlen != 0)
   {
      /* all headers (including the "score" wildcard) are valid */
      if(isHeader(input, columns))
         *sortby = input;
      /* if the score abbreviation is entered, the header would not be in the column info */
      else if( input == SCIPconshdlrDecompGetScoretypeShortName(scip, SCIPconshdlrDecompGetScoretype(scip)))
         *sortby = "score";
   }
   return SCIP_OKAY;
}


static
SCIP_RETCODE GCGdialogExecCommand(
   SCIP*                   scip,          /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr,    /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,        /**< dialog for user input management */
   std::vector<Columninfo*> columns,      /**< list of column headers/ info sources */
   char*                   command,       /**< the command that was entered */
   int*                    startindex,    /**< number of seeed there the menu extract starts */
   int*                    menulength,    /**< current menu length to be modified */
   SCIP_Bool*              finished,      /**< whether to quit the menu */
   std::vector<int>*       idlist,        /**< current list of seeed ids */
   bool*                   sortasc,       /**< true iff sorting should be ascending */
   std::string*            sortby         /**< name of table header to identify sorting column */
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
         if( *startindex > (int) idlist->size() - *menulength )
            *startindex = (int) idlist->size() - *menulength;
      }
      else if( strncmp( command, "top", commandlen) == 0 )
      {
         *startindex = 0;
      }
      else if( strncmp( command, "end", commandlen) == 0 )
      {
         *startindex = (int) idlist->size() - *menulength;
      }

      else if( strncmp( command, "quit", commandlen) == 0 || strncmp( command, "..", commandlen) == 0 )
      {
         *finished = TRUE;
         SCIP_CALL( SCIPconshdlrDecompChooseCandidatesFromSelected(scip) );
      }

      else if( strncmp( command, "legend", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogShowLegend(scip, columns) );
      }

      else if( strncmp( command, "help", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogShowHelp(scip) );
      }

      else if( strncmp( command, "entries", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSetNEntires(scip, dialoghdlr, dialog, (int) idlist->size(), menulength) );
      }

      else if( strncmp( command, "visualize", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSelectVisualize(scip, dialoghdlr, dialog, *idlist) );
      }

      else if( strncmp( command, "inspect", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogInspectSeeed( scip, dialoghdlr, dialog, *idlist) );
      }

      else if( strncmp( command, "select", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSelect(scip, dialoghdlr, dialog, *idlist) );
      }

      else if( strncmp( command, "score", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogChangeScore(scip, dialoghdlr, dialog) );
      }

      else if( strncmp( command, "ascending", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSortAsc(scip, dialoghdlr, dialog, sortasc) );
      }
      else if( strncmp( command, "sort", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSortBy(scip, dialoghdlr, dialog, columns, sortby) );
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
   int startindex = 0;                       /**< number of seeed there the menu extract starts */
   int menulength = DEFAULT_MENULENGTH;      /**< number of entries shown in menu */
   bool sortasc = false;                     /**< whether to show entries in ascending order (score) */
   std::string sortby = DEFAULT_SORT_HEADER; /**< table header, identifies by which column to sort by */

   /* check for available seeeds */
   int nseeeds;   /**< stores the last known number of seeeds, is handed down to check for changes in seeed number */
   nseeeds = SCIPconshdlrDecompGetNSeeeds(scip);
   if(nseeeds == 0)
   {
      SCIPdialogMessage( scip, NULL, "There are no decompositions to explore yet, please detect first.\n" );
      return SCIP_OKAY;
   }

   /* get initial seeed id list */
   int* idarray;
   int listlength;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idarray, nseeeds) );
   SCIPconshdlrDecompGetSeeedLeafList(scip, &idarray, &listlength);

   /* put ids into vector for easier handling */
   std::vector<int> idlist = std::vector<int>();
   for(int i = 0; i < listlength; i++)
   {
      idlist.push_back(idarray[i]);
   }

   /* free idarray */
   SCIPfreeBlockMemoryArray(scip, &idarray, nseeeds);

   /* set initial columns */
   /* the column information has the header, a getter for the column info and the return type of the getter */
   std::vector<Columninfo*> columns;
   char columnstr[] = DEFAULT_COLUMNS;
   char* tempcolumns = strtok(columnstr, " ");
   /* go through each column header and determine its getter */
   while(tempcolumns != NULL)
   {
      /* get each column header of default */
      char newchar[DEFAULT_COLUMN_MAX_WIDTH]; // cutting string at max column width if longer
      strcpy(newchar, tempcolumns);
      /**@note: 'score' is a wildcard! replace by score name later*/

      /* determine what callback the header should receive as its getter */
      if( strcmp(newchar, "nr") == 0)
         /* "nr" represents the position in the menu table and is determined by the menu */
         columns.push_back(new Columninfo(newchar, "number of the decomposition (use this number for selecting the decomposition)", NULL, UNKNOWN));
      else if(strcmp(newchar, "id") == 0)
         /* "id" is the seeed id, the list of ids is known to the menu */
         columns.push_back(new Columninfo(newchar, "id of the decomposition (identifies the decomposition in reports/statistics/visualizations/etc.)", NULL, UNKNOWN));
      else
      {
         /**@note devs: if you want to add new headers, please specify their getters here! */
         RETTYPE type = UNKNOWN;
         callback funct;
         std::string desc;

         if(strcmp(newchar, "nbloc") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetNBlocksBySeeedId;
            type = INTEGER;
            desc = "number of blocks";
         }
         else if(strcmp(newchar, "nmacon") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetNMasterConssBySeeedId;
            type = INTEGER;
            desc = "number of master constraints";
         }
         else if(strcmp(newchar, "nmavar") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetNMasterVarsBySeeedId;
            type = INTEGER;
            desc = "number of \"master only\" variables (also called \"static\", do not occur in blocks)";
         }
         else if(strcmp(newchar, "nlivar") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetNLinkingVarsBySeeedId;
            type = INTEGER;
            desc = "number of linking variables";
         }
         else if(strcmp(newchar, "nstlva") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetNStairlinkingVarsBySeeedId;
            type = INTEGER;
            desc = "number of stair linking variables";
         }
         else if(strcmp(newchar, "score") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetScoreBySeeedId;
            type = FLOAT;
            /**@note as "score" is a wildcard, its description is determined only when needed */
            desc = " ";
         }
         else if(strcmp(newchar, "history") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetDetectorHistoryBySeeedId;
            type = STRING;
            desc = "list of detectors (their chars) which  worked on this decomposition ";
         }
         else if(strcmp(newchar, "pre") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGisPresolvedBySeeedId;
            type = BOOLEAN;
            desc = "is this decomposition for the presolved problem?";
         }
         else if(strcmp(newchar, "nopcon") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetNOpenConssBySeeedId;
            type = INTEGER;
            desc = "number of open (=unassigned) constraints";
         }
         else if(strcmp(newchar, "nopvar") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGgetNOpenVarsBySeeedId;
            type = INTEGER;
            desc = "number of open (=unassigned) variables";
         }
         else if(strcmp(newchar, "sel") == 0)
         {
            funct = (void(*)(SCIP*, int)) &GCGisSelectedBySeeedId;
            type = BOOLEAN;
            desc = "is this decomposition selected?";
         }

         /* add the column if a corresponding callback was found */
         if(type != UNKNOWN)
            columns.push_back(new Columninfo(newchar, desc, funct, type));
      }

      /* get the next item in the list of headers */
      tempcolumns = strtok (NULL, " ");
   }

   /* check that the given default sorting header is valid */
   assert(isHeader(sortby, columns));

   /* sort by default, descending */
   sortSeeedList(scip, &idlist, sortby, columns, sortasc);

   /* while user has not aborted: show current list extract and catch commands */
   SCIP_Bool finished = false;
   char* command;
   SCIP_Bool endoffile;
   while( !finished )
   {
      GCGdialogShowMenu(scip, columns, &nseeeds, startindex, menulength, &idlist, sortasc, sortby);

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Please enter command   (or \"h\" for help) : \nGCG/explore> ", &command, &endoffile) );

      GCGdialogExecCommand(scip, dialoghdlr, dialog, columns, command, &startindex, &menulength, &finished, &idlist, &sortasc, &sortby);
   }

   return SCIP_OKAY;
}

} // extern "C"

} // namespace gcg
