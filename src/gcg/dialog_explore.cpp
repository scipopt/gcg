/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dialog_explore.cpp
 * @brief  dialog menu for exploring decompositions
 * @author Michael Bastubbe
 * @author Hanna Franzen
 * @author Erik Muehmer
 *
 * This file contains all dialog calls to build and use the explore menu.
 * The explore menu gives the user detailed information about all decompositions and a possibility to edit such.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <cstddef>
#include <cstdlib>

#include <string>
#include <iostream>
#include <regex>
#include <map>
#include <sstream>
#include <iomanip>
#include <functional>

#include "gcg/class_partialdecomp.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/score.h"

/* column headers */
#define DEFAULT_COLUMN_MIN_WIDTH  4 /**< min width of a column in the menu table */
#define DEFAULT_SORT_HEADER "score"

#define DEFAULT_MENULENGTH 10 /**< initial number of entries in menu */

/**< default column headers */
const char* DEFAULT_COLUMNS[] {"nr", "id", "nbloc", "nmacon", "nlivar", "nmavar", "nstlva", "score", "history", "pre", "nopcon", "nopvar", "sel"};

namespace gcg
{

/** RETTYPE is used to store the return type of a callback function */
enum RETTYPE
{
   UNKNOWN,    /**< dummy default */
   INTEGER,    /**< integer */
   REAL,       /**< SCIP_Real */
   BOOLEAN,    /**< Boolean */
   STRING      /**< char* */
};

/** storage for column information */
class AbstractColumn
{
public:
   AbstractColumn(const char* columnHeader, const char* columnDesc, RETTYPE columnType) :
         header(columnHeader), desc(columnDesc), type(columnType) {}

   virtual ~AbstractColumn() = default;

   virtual int compareValues(PARTIALDECOMP* firstdec, PARTIALDECOMP* seconddec) = 0;

   virtual std::string getValueAsString(PARTIALDECOMP* partialdec) = 0;

   RETTYPE getReturnType() const { return type; }

   std::string header;  /**< table header of the column */
   std::string desc;    /**< description of the column entries used in the menu help */

private:
   RETTYPE type;        /**< return type of the getter */
};

template <class T>
class Column: public AbstractColumn
{
public:
   Column(const char* columnHeader, const char* columnDesc, std::function<T(PARTIALDECOMP&)>&& columnCallback, RETTYPE columnType) :
         AbstractColumn(columnHeader, columnDesc, columnType), callback(columnCallback) {}

   ~Column() override = default;

   T getValue(PARTIALDECOMP* partialdec) { return callback(*partialdec); }

   std::string getValueAsString(PARTIALDECOMP* partialdec) override
   {
      std::ostringstream sstream;
      sstream << std::fixed << std::setprecision(2) << getValue(partialdec);
      return sstream.str();
   }

   int compareValues(PARTIALDECOMP* firstdec, PARTIALDECOMP* seconddec) override
   {
       if( getReturnType() == UNKNOWN )
       {
          return 0;
       }
       else
       {
          T val1 = getValue(firstdec);
          T val2 = getValue(seconddec);
          if( val1 == val2 )
             return 0;
          else if( val1 < val2 )
             return -1;
          else
             return 1;
       }
   }

private:
   std::function<T(PARTIALDECOMP&)> callback;
};

bool updatePartialdecList(
   GCG* gcg,
   std::vector<PARTIALDECOMP*>& partialdeclist,
   unsigned int& npartialdecs,
   bool includeopenpartialdecs
   )
{
   unsigned int newnpartialdecs = GCGconshdlrDecompGetNPartialdecs(gcg);
   if( newnpartialdecs != npartialdecs || (includeopenpartialdecs != (npartialdecs == partialdeclist.size())) )
   {
      std::vector<PARTIALDECOMP*>* partialdecs = GCGconshdlrDecompGetPartialdecs(gcg);

      npartialdecs = newnpartialdecs;
      partialdeclist.clear();
      partialdeclist.reserve(npartialdecs);

      for( PARTIALDECOMP* partialdec: *partialdecs )
      {
         if( includeopenpartialdecs || partialdec->isComplete() )
         {
            partialdeclist.push_back(partialdec);
         }
      }
      return true;
   }
   return false;
}

/** @brief local sorting function for partialdec id vectors
 *
 * avoids redundant sorting calls,
 * sorts by score in given order
 */
static
void sortPartialdecList(
   std::vector<PARTIALDECOMP*>& partialdeclist,    /**< current list of partialdecs */
   std::string& header,                            /**< header of column to sort by */
   std::vector<AbstractColumn*>& columns,          /**< vector of pointers to columns */
   bool asc                                        /**< whether to sort ascending or descending */
   )
{
   /* find the column infos for the given header */
   for( auto column : columns )
   {
      if( column->header.find(header) == 0 )
      {
         /* sort the id list according to given order using the callback getter of the column */
         if( column->getReturnType() != UNKNOWN )
         {
            /* the callback has to be parsed to expect an int output */
            if( asc )
               std::stable_sort(partialdeclist.begin(), partialdeclist.end(), [&](PARTIALDECOMP* const a, PARTIALDECOMP* const b) { return column->compareValues(a, b) < 0; });
            else
               std::stable_sort(partialdeclist.begin(), partialdeclist.end(), [&](PARTIALDECOMP* const a, PARTIALDECOMP* const b) { return column->compareValues(b, a) < 0; });
         }
         break;
      }
   }
}


/** modifies menulength according to input and updates menu accordingly
 * @returns SCIP return code */
static
SCIP_RETCODE GCGdialogSetNEntires(
   GCG* gcg,                     /**< GCG data structure */
   SCIP_DIALOGHDLR* dialoghdlr,  /**< dialog handler for user input management */
   SCIP_DIALOG* dialog,          /**< dialog for user input management */
   int listlength,               /**< length of partialdec id list */
   int& menulength               /**< current menu length to be modified */
   )
{
   SCIP* scip;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int newlength;
   int commandlen;

   scip = GCGgetOrigprob(gcg);

   SCIPdialogMessage(scip, NULL, "Please specify the amount of entries to be shown in this menu:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = (int) strlen(ntovisualize);

   newlength = -1;
   if( commandlen != 0 )
      newlength = atoi(ntovisualize);

   /* check whether there are decompositions,
    * (preventing "Why doesn't it show anything? Maybe the entry number is 0") */
   if( GCGconshdlrDecompGetNPartialdecs(gcg) == 0 )
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
      menulength = newlength;
   else
      menulength = listlength;

   return SCIP_OKAY;
}


/** sets the used score according to user input
 *
 * @returns SCIP return code
 */
static
SCIP_RETCODE GCGdialogChangeScore(
   GCG* gcg,                           /**< GCG data structure */
   SCIP_DIALOGHDLR* dialoghdlr,        /**< dialog handler for user input management */
   SCIP_DIALOG* dialog                 /**< dialog for user input management */
   )
{
   SCIP* scip;
   char* getscore;
   SCIP_Bool endoffile;
   int commandlen;

   scip = GCGgetOrigprob(gcg);

   SCIPdialogMessage(scip, NULL, "\nPlease specify the new score:\n");
   for( int i = 0; i < GCGconshdlrDecompGetNScores(gcg); i++)
   {
      GCG_SCORE* score = GCGconshdlrDecompGetScores(gcg)[i];

      SCIPdialogMessage(scip, NULL, "%d: %s\n", i, GCGscoreGetName(score));      
   }
   SCIPdialogMessage(scip, NULL, "Note: Sets the detection/scores/selected parameter to the score\'s shortname.\n");

   /* get input */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &getscore, &endoffile) );
   commandlen = (int) strlen(getscore);
   if( commandlen != 0 )
   {
      /* convert to int (if value is invalid, this results in 0) */
      int scorenr = atoi(getscore);

      /* check if the value is in valid range */
      if( scorenr >= 0 && scorenr <= GCGconshdlrDecompGetNScores(gcg) - 1 )
      {
         /* set score */
         GCG_SCORE* score = GCGconshdlrDecompGetScores(gcg)[scorenr];

         SCIP_CALL( SCIPsetStringParam(scip, "detection/scores/selected", GCGscoreGetShortname(score)) );
         SCIPdialogMessage(scip, NULL, "Score set to %s.\n", GCGscoreGetName(score));
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

/** @brief show current menu containing partialdec information
 *
 * Update partialdec list in case it changed since the last command
 * and show the table of partialdecs.
 * @returns SCIP status
 */
static
SCIP_RETCODE GCGdialogShowMenu(
   GCG* gcg,                                       /**< GCG data structure */
   std::vector<AbstractColumn*>& columns,          /**< vector of pointers to columns */
   unsigned int& npartialdecs,                     /**< max number of partialdecs */
   const int startindex,                           /**< index (in partialdec list) of uppermost partialdec in extract */
   int menulength,                                 /**< number of menu entries */
   std::vector<PARTIALDECOMP*>& partialdeclist,    /**< current list of partialdecs */
   bool sortasc,                                   /**< true iff sorting should be ascending */
   std::string sortby,                             /**< table header of column to sort by */
   bool listopenpartialdecs                        /**< open partialdecs will be listed iff set to true*/
   )
{
   SCIP* scip;
   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

   /* update partialdec list in case it changed (in which case the amount of partialdecs should have changed)*/
   updatePartialdecList(gcg, partialdeclist, npartialdecs, listopenpartialdecs);

   /* sort partialdec ids by score, descending (in case score was changed or id list was updated)*/
   sortPartialdecList(partialdeclist, sortby, columns, sortasc);

   /* count corresponding partialdecs for overview statistics */
   int ndetectedpresolved = 0;
   int ndetectedunpresolved = 0;

   for(PARTIALDECOMP* partialdec : partialdeclist)
   {
      if(partialdec->isComplete())
      {
         /* from presolved problem */
         if(partialdec->isAssignedToOrigProb())
         {
            ++ndetectedunpresolved;
         }
         /* from original problem */
         else
         {
            ++ndetectedpresolved;
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
   for( auto column : columns )
   {
      /* "score" is a wildcard for the current score, relace it with actual scoretype */
      std::string header = column->header;
      std::string newheader;
      if( header != "score" )
         newheader = header;
      else
         newheader = GCGscoreGetShortname(GCGgetCurrentScore(gcg));

      /* make sure the header name is unique and add a length for header */
      assert(columnlength.find(header) == columnlength.end());
      columnlength.insert(std::pair<std::string,int>(header, 0));
      /* if header is smaller than min column width, add spaces to header first */
      if( newheader.size() < DEFAULT_COLUMN_MIN_WIDTH )
      {
         for( int i = 0; i < (DEFAULT_COLUMN_MIN_WIDTH - (int) newheader.size()); i++ )
         {
            headerline += " ";
            borderline += "-";
            columnlength.at(header)++;
         }
      }
      /* add header to headerline and add #chars of header as '-' to borderline*/
      headerline += newheader;
      for( int i = 0; i < (int) newheader.size(); i++ )
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

   /* go through all partialdecs that should currently be displayed,
    * so from startindex on menulength many entries if there are that much left in the list */
   for( int i = startindex; i < startindex + menulength && i < (int) partialdeclist.size(); ++i )
   {
      /* get current partialdec id */
      PARTIALDECOMP* partialdec = partialdeclist.at(i);

      /* each line starts with a space */
      SCIPdialogMessage(scip, NULL, " ");

      /* go through the columns and write the entry for each one */
      for( auto column : columns )
      {
         std::ostringstream sstream;
         std::string header = column->header;

         /* write spaces to fill out the columnwidth */
         sstream << std::setw(columnlength.at(header));

         /* call the getter of the column depending on the given return type */
         if( header == "nr" )
         {
            sstream << i;
         }
         else
         {
            sstream << column->getValueAsString(partialdec);
         }
         SCIPdialogMessage(scip, NULL, "%s ", sstream.str().c_str());
      }

      /* continue to next line */
      SCIPdialogMessage(scip, NULL, "\n");
   }

   /* at the end of the table add a line */
   outputCharXTimes(scip, '=', linelength);
   SCIPdialogMessage(scip, NULL, "\n");

   return SCIP_OKAY;
}


/** Shows information about the explore screen and its abbreviations
 *
 * @returns SCIP status */
static
SCIP_RETCODE GCGdialogShowLegend(
   GCG* gcg,                              /**< GCG data structure */
   std::vector<AbstractColumn*>& columns  /**< vector of pointers to columns */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   GCG_DETECTOR** detectors;

   scip = GCGgetOrigprob(gcg);

   /* print header for detector list */
   SCIPdialogMessage(scip, NULL, "List of included detectors for decompositions histories: \n");

   SCIPdialogMessage(scip, NULL, "\n%30s    %4s\n", "detector" , "char");
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "--------" , "----");

   /* get and print char of each detector */
   detectors = GCGconshdlrDecompGetDetectors(gcg);

   for( int det = 0; det < GCGconshdlrDecompGetNDetectors(gcg); ++det )
   {
      GCG_DETECTOR* detector;
      detector = detectors[det];

      SCIPdialogMessage(scip, NULL, "%30s    %4c\n", GCGdetectorGetName(detector), GCGdetectorGetChar(detector));
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
   for( auto column : columns )
   {
      /* get table header */
      std::string& header = column->header;

      /* print the header with the description */
      if( header != "score" )
      {
         SCIPdialogMessage(scip, NULL, "%30s     %s\n", header.c_str(), column->desc.c_str());
      }
      /* if the header is "score" replace with shortname of the current score */
      else
      {
         GCG_SCORE* score = GCGgetCurrentScore(gcg);
         SCIPdialogMessage(scip, NULL, "%30s     %s\n", GCGscoreGetShortname(score), GCGscoreGetDesc(score));
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
   GCG* gcg  /**< GCG data structure */
   )
{
   SCIP* scip;
   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

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
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "export", "generates visualization of the specified decomposition in gnuplot format");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "generates visualization and opens it (requires gnuplot)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "inspect", "displays detailed information for the specified decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "score", "sets the score by which the quality of decompositions is evaluated");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "sort", "sets the column by which the decompositions are sorted (default: by score)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "ascending", "sort decompositions in ascending (true) or descending (false) order");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "list", "specify whether all decompositions should be listed");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "return to main menu");

   SCIPdialogMessage(scip, NULL, "\n=================================================================================================== \n");

   return SCIP_OKAY;
}


/** Shows a visualization of the partialdec specified by the user via the dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE GCGdialogSelectVisualize(
   GCG*                          gcg,            /**< GCG data structure */
   SCIP_DIALOGHDLR*              dialoghdlr,     /**< dialog handler for user input management */
   SCIP_DIALOG*                  dialog,         /**< dialog for user input management */
   std::vector<PARTIALDECOMP*>&  partialdeclist, /**< current list of partialdecs */
   SCIP_Bool                     open            /**< compile and open gnuplot file? */
   )
{
   SCIP* scip;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;
   int commandlen;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

   SCIPdialogMessage(scip, NULL, "Please specify the nr of the decomposition to be visualized:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = (int) strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0 )
      idtovisu = atoi(ntovisualize);

   /* check whether the partialdec exists */
   if( commandlen == 0 || idtovisu < 0 || idtovisu >= (int) partialdeclist.size() )
   {
      SCIPdialogMessage( scip, NULL, "This nr is out of range." );
      return SCIP_OKAY;
   }

   /* get and show partialdec */
   PARTIALDECOMP* partialdec = partialdeclist.at(idtovisu);
   if( open ) 
      partialdec->showVisualization();
   else
      partialdec->exportVisualization();

   return SCIP_OKAY;
}


/**
 * Displays information about a partialdec that is chosen by the user in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE GCGdialogInspectPartialdec(
   GCG*                          gcg,            /**< GCG data structure */
   SCIP_DIALOGHDLR*              dialoghdlr,     /**< dialog handler for user input management */
   SCIP_DIALOG*                  dialog,         /**< dialog for user input management */
   std::vector<PARTIALDECOMP*>&  partialdeclist  /**< current list of partialdecs */
   )
{
   SCIP* scip;
   char* ntoinspect;
   char* ndetaillevel;
   SCIP_Bool endoffile;
   int idtoinspect;
   int detaillevel;

   int commandlen;

   assert( gcg != NULL );

   scip = GCGgetOrigprob(gcg);

   /* read the id of the decomposition to be inspected */
   SCIPdialogMessage( scip, NULL, "Please specify the nr of the decomposition to be inspected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntoinspect, &endoffile ) );
   commandlen = (int) strlen( ntoinspect );

   idtoinspect = -1;
   if( commandlen != 0 )
      idtoinspect = atoi( ntoinspect );

   if(idtoinspect < 0 || idtoinspect >= (int) partialdeclist.size()){
      SCIPdialogMessage( scip, NULL, "This nr is out of range." );
      return SCIP_OKAY;
   }

   /* check whether ID is in valid range */
   PARTIALDECOMP* partialdec = partialdeclist.at(idtoinspect);

   /* read the desired detail level; for wrong input, it is set to 1 by default */
   SCIPdialogMessage( scip, NULL,
      "Please specify the detail level:\n  0 - brief overview\n  1 - block and detector info (default)\n  2 - cons and var assignments\n" );
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ndetaillevel, &endoffile ) );
   commandlen = (int) strlen( ndetaillevel );

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
   partialdec->displayInfo( detaillevel );

   return SCIP_OKAY;
}


/** Lets the user select decompositions from the explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE GCGdialogSelect(
   GCG*                          gcg,            /**< GCG data structure */
   SCIP_DIALOGHDLR*              dialoghdlr,     /**< dialog handler for user input management */
   SCIP_DIALOG*                  dialog,         /**< dialog for user input management */
   std::vector<PARTIALDECOMP*>&  partialdeclist  /**< current list of partialdecs */
   )
{
   SCIP* scip;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;

   int commandlen;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

   /* get input */
   SCIPdialogMessage(scip, NULL, "\nPlease specify the nr of the decomposition to be selected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = (int) strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   /* check if the input is a valid number */
   if( commandlen == 0 || idtovisu < 0 || idtovisu >= (int) partialdeclist.size() )
   {
      SCIPdialogMessage( scip, NULL, "This nr is out of range, nothing was selected." );
      return SCIP_OKAY;
   }

   /* get partialdec from id*/
   PARTIALDECOMP* partialdec = partialdeclist.at(idtovisu);

   /* reverse selection (deselects if partialdec was previously selected) */
   partialdec->setSelected(!partialdec->isSelected() );

   return SCIP_OKAY;
}

/** Set whether order in menu should be ascending/descending
 *
 * @returns SCIP return code */
static
SCIP_RETCODE GCGdialogSortAsc(
   GCG*                    gcg,        /**< GCG data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   bool&                   asc         /**< true iff sorting should be ascending */
   )
{
   SCIP* scip;
   char* ascen;
   SCIP_Bool endoffile;
   int commandlen;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

   /* get input */
   SCIPdialogMessage(scip, NULL, "\nPlease enter \"true\"/\"1\" for ascending or \"false\"/\"0\" for descending order:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ascen, &endoffile) );
   commandlen = (int) strlen(ascen);

   std::string input = ascen;
   if( commandlen != 0)
   {
      /* check if input if is true, set ascending */
      if(input == "true" || input == "1")
         asc = true;
      /* check if input if is false, set descending */
      else if(input == "false" || input == "0")
         asc = false;
      /* all other inputs are considered invalid and do nothing */
   }

   return SCIP_OKAY;
}


/** Checks whether the given header is valid
 * @returns true iff header is valid
 */
static
bool isHeader(
   std::string& header,                  /**< header to check */
   std::vector<AbstractColumn*>& columns /**< vector of pointers to columns */
   )
{
   /* check if the given header is a (prefix of a) registered table header */
   for( auto column : columns )
   {
      if(column->header.find(header) == 0 )
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
   GCG*                          gcg,        /**< GCG data structure */
   SCIP_DIALOGHDLR*              dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*                  dialog,     /**< dialog for user input management */
   std::vector<AbstractColumn*>& columns,    /**< vector of pointers to columns */
   std::string&                  sortby      /**< table header, identifies by which column to sort by */
   )
{
   SCIP* scip;
   char* newsort;
   SCIP_Bool endoffile;
   int commandlen;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

   /* get input */
   SCIPdialogMessage(scip, NULL, "\nPlease enter the table header of the column by which you would like to sort:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &newsort, &endoffile) );
   commandlen = (int) strlen(newsort);

   /* if the input is a valid table header, change sortby */
   std::string input = newsort;
   if( commandlen != 0 )
   {
      /* all headers (including the "score" wildcard) are valid */
      if( isHeader(input, columns) )
         sortby = input;
      /* if the score abbreviation is entered, the header would not be in the column info */
      else if( input == GCGscoreGetShortname(GCGgetCurrentScore(gcg)) )
         sortby = "score";
   }
   return SCIP_OKAY;
}


/** Set whether order in menu should be ascending/descending
 *
 * @returns SCIP return code */
static
SCIP_RETCODE GCGdialogChangeListMode(
   GCG*                    gcg,        /**< GCG data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   bool&                   listopenpartialdecs /** will be updated */
   )
{
   SCIP* scip;
   char* input;
   SCIP_Bool endoffile;
   int commandlen;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

   /* get input */
   SCIPdialogMessage(scip, NULL, "\nShould incomplete decompositions be listed? Please enter \"true\" or \"false\":\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &input, &endoffile) );
   commandlen = (int) strlen(input);

   if( commandlen != 0)
   {
      /* check if input if is true, set ascending */
      if( strncmp(input, "true", commandlen) == 0 || strncmp(input, "1", commandlen) == 0 )
         listopenpartialdecs = true;
         /* check if input if is false, set descending */
      else if( strncmp(input, "false", commandlen) == 0 || strncmp(input, "0", commandlen) == 0 )
         listopenpartialdecs = false;
      /* all other inputs are considered invalid and do nothing */
   }
   return SCIP_OKAY;
}


static
SCIP_RETCODE GCGdialogExecCommand(
   GCG*                          gcg,                    /**< GCG data structure */
   SCIP_DIALOGHDLR*              dialoghdlr,             /**< dialog handler for user input management */
   SCIP_DIALOG*                  dialog,                 /**< dialog for user input management */
   std::vector<AbstractColumn*>& columns,                /**< vector of pointers to columns */
   char*                         command,                /**< the command that was entered */
   int&                          startindex,             /**< number of partialdec there the menu extract starts */
   int&                          menulength,             /**< current menu length to be modified */
   bool&                         finished,               /**< whether to quit the menu */
   std::vector<PARTIALDECOMP*>&  partialdeclist,         /**< current list of partialdec ids */
   bool&                         sortasc,                /**< true iff sorting should be ascending */
   std::string&                  sortby,                 /**< name of table header to identify sorting column */
   bool&                         listopenpartialdecs     /** whether open patialdecs should be listed*/
   )
{
   int commandlen = (int) strlen(command);

      if( strncmp( command, "previous", commandlen) == 0 )
      {
         startindex = startindex - menulength;
         if(startindex < 0 )
            startindex = 0;
      }
      else if( strncmp( command, "next", commandlen) == 0 )
      {
         startindex = startindex + menulength;
         if( startindex > (int) partialdeclist.size() - menulength )
            startindex = (int) partialdeclist.size() - menulength;
      }
      else if( strncmp( command, "top", commandlen) == 0 )
      {
         startindex = 0;
      }
      else if( strncmp( command, "end", commandlen) == 0 )
      {
         startindex = (int) partialdeclist.size() - menulength;
         if (startindex < 0)
            startindex = 0;
      }

      else if( strncmp( command, "quit", commandlen) == 0 || strncmp( command, "..", commandlen) == 0 )
      {
         finished = TRUE;
      }

      else if( strncmp( command, "legend", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogShowLegend(gcg, columns) );
      }

      else if( strncmp( command, "help", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogShowHelp(gcg) );
      }

      else if( strncmp( command, "entries", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSetNEntires(gcg, dialoghdlr, dialog, (int) partialdeclist.size(), menulength) );
      }

      else if( strncmp( command, "visualize", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSelectVisualize(gcg, dialoghdlr, dialog, partialdeclist, TRUE) );
      }

      else if( strncmp( command, "export", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSelectVisualize(gcg, dialoghdlr, dialog, partialdeclist, FALSE) );
      }

      else if( strncmp( command, "inspect", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogInspectPartialdec(gcg, dialoghdlr, dialog, partialdeclist) );
      }

      else if( strncmp( command, "select", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSelect(gcg, dialoghdlr, dialog, partialdeclist) );
      }

      else if( strncmp( command, "score", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogChangeScore(gcg, dialoghdlr, dialog) );
      }

      else if( strncmp( command, "ascending", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSortAsc(gcg, dialoghdlr, dialog, sortasc) );
      }
      else if( strncmp( command, "sort", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogSortBy(gcg, dialoghdlr, dialog, columns, sortby) );
      }
      else if( strncmp( command, "list", commandlen) == 0 )
      {
         SCIP_CALL( GCGdialogChangeListMode(gcg, dialoghdlr, dialog, listopenpartialdecs) );
      }

   return SCIP_OKAY;
}

extern "C" {

SCIP_RETCODE GCGdialogExecExplore(
   GCG*                    gcg,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   )
{
   /* set navigation defaults */
   int startindex = 0;                       /**< number of partialdec there the menu extract starts */
   int menulength = DEFAULT_MENULENGTH;      /**< number of entries shown in menu */
   bool sortasc = false;                     /**< whether to show entries in ascending order (score) */
   std::string sortby = DEFAULT_SORT_HEADER; /**< table header, identifies by which column to sort by */
   bool listopenpartialdecs = false;

   /* check for available partialdecs */
   unsigned int npartialdecs = 0;   /**< stores the last known number of partialdecs, is handed down to check for changes in partialdec number */
   std::vector<PARTIALDECOMP*> partialdeclist;
   updatePartialdecList(gcg, partialdeclist, npartialdecs, listopenpartialdecs);

   if(npartialdecs == 0)
   {
      SCIPdialogMessage(GCGgetOrigprob(gcg), NULL, "There are no decompositions to explore yet, please detect first.\n");
      return SCIP_OKAY;
   }

   /* set initial columns */
   /* the column information has the header, a getter for the column info and the return type of the getter */
   std::vector<AbstractColumn*> columns;
   /* go through each column header and determine its getter */
   for( auto columnname : DEFAULT_COLUMNS )
   {
      /**@note devs: if you want to add new headers, please specify their getters here! */
      AbstractColumn* column;

      if( strcmp(columnname, "nr") == 0 )
      {
         /* special case: "nr" represents the position in the menu table and is determined by the menu */
         column = new Column<int>(
            columnname,
            "number of the decomposition (use this number for selecting the decomposition)",
            NULL,
            UNKNOWN);
      }
      else if( strcmp(columnname, "id") == 0 )
      {
         /* "id" is the partialdec id, the list of ids is known to the menu */
         column = new Column<int>(
            columnname,
            "id of the decomposition (identifies the decomposition in reports/statistics/visualizations/etc.)",
            &PARTIALDECOMP::getID,
            INTEGER);
      }
      else if( strcmp(columnname, "nbloc") == 0 )
      {
         column = new Column<int>(
            columnname,
            "number of blocks",
            &PARTIALDECOMP::getNBlocks,
            INTEGER);
      }
      else if( strcmp(columnname, "nmacon") == 0 )
      {
         column = new Column<int>(
            columnname,
            "number of master constraints",
            &PARTIALDECOMP::getNMasterconss,
            INTEGER);
      }
      else if( strcmp(columnname, "nmavar") == 0 )
      {
         column = new Column<int>(
            columnname,
            "number of \"master only\" variables (also called \"static\", do not occur in blocks)",
            &PARTIALDECOMP::getNMastervars,
            INTEGER);
      }
      else if( strcmp(columnname, "nlivar") == 0 )
      {
         column = new Column<int>(
            columnname,
            "number of linking variables",
            &PARTIALDECOMP::getNLinkingvars,
            INTEGER);
      }
      else if( strcmp(columnname, "nstlva") == 0 )
      {
         column = new Column<int>(
            columnname,
            "number of stair linking variables",
            &PARTIALDECOMP::getNTotalStairlinkingvars,
            INTEGER);
      }
      else if( strcmp(columnname, "score") == 0 )
      {
         column = new Column<SCIP_Real>(
            columnname,
            " ",
            static_cast<SCIP_Real(PARTIALDECOMP::*)()>(&PARTIALDECOMP::getScore),
            REAL);
      }
      else if( strcmp(columnname, "history") == 0 )
      {
         column = new Column<std::string>(
            columnname,
            "list of detectors (their chars) which  worked on this decomposition",
            static_cast<std::string(PARTIALDECOMP::*)()>(&PARTIALDECOMP::buildDecChainString),
            STRING);
      }
      else if( strcmp(columnname, "pre") == 0 )
      {
         column = new Column<bool>(
            columnname,
            "is this decomposition for the presolved problem?",
            &PARTIALDECOMP::isAssignedToPresolvedProb,
            BOOLEAN);
      }
      else if( strcmp(columnname, "nopcon") == 0 )
      {
         column = new Column<int>(
            columnname,
            "number of open (=unassigned) constraints",
            &PARTIALDECOMP::getNOpenconss,
            INTEGER);
      }
      else if( strcmp(columnname, "nopvar") == 0 )
      {
         column = new Column<int>(
            columnname,
            "number of open (=unassigned) variables",
            &PARTIALDECOMP::getNOpenvars,
            INTEGER);
      }
      else if( strcmp(columnname, "sel") == 0 )
      {
         column = new Column<bool>(
            columnname,
            "is this decomposition selected?",
            &PARTIALDECOMP::isSelected,
            BOOLEAN);
      }

      columns.push_back(column);
   }

   /* check that the given default sorting header is valid */
   assert(isHeader(sortby, columns));

   /* sort by default, descending */
   sortPartialdecList(partialdeclist, sortby, columns, sortasc);

   /* while user has not aborted: show current list extract and catch commands */
   bool finished = false;
   char* command;
   SCIP_Bool endoffile;
   while( !finished )
   {
      GCGdialogShowMenu(gcg, columns, npartialdecs, startindex, menulength, partialdeclist, sortasc, sortby, listopenpartialdecs);

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Please enter command   (or \"h\" for help) : \nGCG/explore> ", &command, &endoffile) );

      GCGdialogExecCommand(gcg, dialoghdlr, dialog, columns, command, startindex, menulength, finished, partialdeclist, sortasc, sortby, listopenpartialdecs);
   }

   for( auto column : columns )
      delete column;

   return SCIP_OKAY;
}

} // extern "C"

} // namespace gcg
