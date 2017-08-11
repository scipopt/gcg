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

/**@file   reader_gp.cpp
 * @brief  GP file reader writing seeeds to gnuplot files
 * @author Martin Bergner
 * @author Hanna Franzen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include <fstream>

#include "scip/scip.h"

#include "reader_gp.h"
#include "scip_misc.h"
#include "struct_decomp.h"
#include "cons_decomp.h"
#include "pub_decomp.h"
#include "params_visu.h"

#include "class_seeed.h"
#include "class_seeedpool.h"
#include "class_miscvisualization.h"

#define READER_NAME             "gpreader"
#define READER_DESC             "gnuplot file writer for seeed visualization"
#define READER_EXTENSION        "gp"


namespace gcg{

/** write file header with terminal etc. */
static
SCIP_RETCODE writeGpHeader(
   char*                 filename,           /**< filename (including path) to write to */
   const char*           outputname          /**< the filename to which gnuplot should compile the visualization */
   )
{
   std::ofstream ofs;
   ofs.open( filename, std::ofstream::out );

   /* set output format and file */
   ofs << "set terminal pdf" << std::endl;
   ofs << "set output \"" << outputname << ".pdf\"" << std::endl;

   ofs.close();

   return SCIP_OKAY;
}


/* writes gp code to given file that contains a box with given coordinates and color */
static
SCIP_RETCODE drawGpBox(
   char* filename,   /**< filename (including path) to write to */
   int objectid,     /**< id number of box (>0), must be unique */
   int x1,           /**< x value of lower left vertex coordinate */
   int y1,           /**< y value of lower left vertex coordinate */
   int x2,           /**< x value of upper right vertex coordinate */
   int y2,           /**< y value of upper right vertex coordinate */
   char* color       /**< color hex code (e.g. #000000) for box filling */
   )
{
   std::ofstream ofs;
   ofs.open( filename, std::ofstream::out );

   ofs << "set object " << objectid << " rect from " << x1 << "," << y1 << " to " << x2 << "," << y2 << " fc rgb \""
      << color << "\"" << " lc rgb \"" << SCIPvisuGetColorLine() << "\"" << std::endl;

   ofs.close();
   return SCIP_OKAY;
}


/** writes gp code to given file that contains all nonzero points */
static
SCIP_RETCODE writeGpNonzeros(
   const char* filename,   /**< filename to write to (including path & extension) */
   Seeed* seeed,           /**< Seeed for which the nonzeros should be visualized */
   Seeedpool* seeedpool,   /**< current Seeedpool */
   float radius            /**< radius of the dots */
   )
{
   std::vector<int> orderToRows(seeed->getNConss(), -1);
   std::vector<int> rowToOrder(seeed->getNConss(), -1);
   std::vector<int> orderToCols(seeed->getNVars(), -1);
   std::vector<int> colsToOrder(seeed->getNVars(), -1);
   int counterrows = 0;
   int countercols = 0;
   std::ofstream ofs;

   /* --- order of constraints --- */

   /* master constraints */
   for ( int i = 0; i < seeed->getNMasterconss() ; ++i )
   {
      int rowidx = seeed->getMasterconss()[i];
      orderToRows[counterrows] = rowidx;
      rowToOrder[rowidx] = counterrows;
      ++counterrows;
   }

   /* block constraints */
   for ( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      for (int i = 0; i < seeed->getNConssForBlock(b) ; ++i )
      {
         int rowidx = seeed->getConssForBlock(b)[i];
         orderToRows[counterrows] = rowidx;
         rowToOrder[rowidx] = counterrows;
         ++counterrows;
      }
   }

   /** open constraints */
   for ( int i = 0; i < seeed->getNOpenconss() ; ++i )
   {
      int rowidx = seeed->getOpenconss()[i];
      orderToRows[counterrows] = rowidx;
      rowToOrder[rowidx] = counterrows;
      ++counterrows;
   }

   /* --- order of variables --- */

   /* linking variables */
   for ( int i = 0; i < seeed->getNLinkingvars() ; ++i )
   {
      int colidx = seeed->getLinkingvars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* master variables */
   for ( int i = 0; i < seeed->getNMastervars() ; ++i )
   {
      int colidx = seeed->getMastervars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* block variables */
   for ( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      for (int i = 0; i < seeed->getNVarsForBlock(b) ; ++i )
      {
         int colidx = seeed->getVarsForBlock(b)[i];
         orderToCols[countercols] = colidx;
         colsToOrder[colidx] = countercols;
         ++countercols;
      }
      for (int i = 0; i < seeed->getNStairlinkingvars(b) ; ++i )
      {
         int colidx = seeed->getStairlinkingvars(b)[i];
         orderToCols[countercols] = colidx;
         colsToOrder[colidx] = countercols;
         ++countercols;
      }
   }

   /** open vars */
   for ( int i = 0; i < seeed->getNOpenvars() ; ++i )
   {
      int colidx = seeed->getOpenvars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   ofs.open (filename, std::ofstream::out );

   /* start writing dots */
   ofs << "plot \"-\" using 1:2:(" << radius << ") notitle with circles fc rgb \"" << SCIPvisuGetColorNonzero()
      << "\" fill solid" << std::endl;

   /* write scatter plot */
   for( int row = 0; row < seeed->getNConss(); ++row )
      for ( int col = 0; col < seeed->getNVars(); ++col )
      {
         assert( orderToRows[row] != -1);
         assert( orderToCols[col] != -1);
         if( seeedpool->getVal( orderToRows[row], orderToCols[col]  ) != 0 )
            ofs << col+0.5 << " " << row+0.5 << std::endl;
      }

   /* end writing dots */
   ofs << "e" << std::endl;

   ofs.close();

   return SCIP_OKAY;
}


static
SCIP_RETCODE writeGpSeeed(
   char* filename,         /**< filename (including path) to write to */
   Seeed* seeed,           /**< Seeed for which the nonzeros should be visualized */
   Seeedpool* seeedpool    /**< current Seeedpool */
   )
{
   int rowboxcounter = 0;
   int colboxcounter = 0;

   std::ofstream ofs;
   ofs.open( filename, std::ofstream::out );

   /* set coordinate range */
   ofs << "set xrange [-1:" << seeed->getNVars() << "]" << std::endl;
   ofs << "set yrange[" << seeed->getNConss() << ":-1]" << std::endl;

   /* --- draw boxes ---*/

   /* linking vars */
   drawGpBox( filename,1, 0, 0, seeed->getNLinkingvars(), seeed->getNConss(), SCIPvisuGetColorLinking() );
   colboxcounter += seeed->getNLinkingvars();

   /* mastervars */
   drawGpBox( filename, 2, colboxcounter, 0, seeed->getNMastervars()+colboxcounter, seeed->getNConss(),
      SCIPvisuGetColorMastervars() );
   colboxcounter += seeed->getNMastervars();

   /* specific interface methodsmasterconss */
   drawGpBox( filename, 3, 0, 0, seeed->getNVars(), seeed->getNMasterconss(), SCIPvisuGetColorMasterconss() );
   rowboxcounter += seeed->getNMasterconss();

   /* blocks */
   for( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      drawGpBox( filename, 2 * b + 4, colboxcounter, rowboxcounter, colboxcounter + seeed->getNVarsForBlock(b),
         rowboxcounter + seeed->getNConssForBlock(b), SCIPvisuGetColorBlock() );
      colboxcounter += seeed->getNVarsForBlock(b);

      if( seeed->getNStairlinkingvars(b) != 0 )
      {
         drawGpBox( filename, 2 * b + 5, colboxcounter, rowboxcounter, colboxcounter + seeed->getNStairlinkingvars(b),
            rowboxcounter + seeed->getNConssForBlock(b) + seeed->getNConssForBlock(b+1), SCIPvisuGetColorStairlinking() );
      }
      colboxcounter += seeed->getNStairlinkingvars(b);
      rowboxcounter += seeed->getNConssForBlock(b);
   }

   /* open */
   drawGpBox( filename, 2 * seeed->getNBlocks() + 4, colboxcounter, rowboxcounter, colboxcounter + seeed->getNOpenvars(),
      rowboxcounter+seeed->getNOpenconss(), SCIPvisuGetColorOpen() );
   colboxcounter += seeed->getNOpenvars();
   rowboxcounter += seeed->getNOpenconss();

   /* --- draw nonzeros --- */
   /*@todo get radius */
   writeGpNonzeros( filename, seeed, seeedpool, SCIPvisuGetNonzeroRadius( seeed->getNVars(), seeed->getNConss(), 1 ) );

   ofs.close();

   return SCIP_OKAY;
}


/*
 * public methods
 */


/** writes a visualization for the given seeed */
SCIP_RETCODE GCGwriteGpVisualization(
   SCIP* scip,             /**< SCIP data structure */
   char* filename,         /**< filename (including path) to write to */
   char* outputname,       /**< filename for compiled output file */
   int seeedid             /**< id of seeed to visualize */
   )
{
   SeeedPtr seeed;
   Seeedpool* seeedpool;
   MiscVisualization* misc = new MiscVisualization();

   /* get seeed and seeedpool */
   seeed = misc->GCGgetSeeed(scip, seeedid, seeedpool);
   if( seeed == NULL )
   {
      SCIPerrorMessage("Could not find Seeed!\n");
      return SCIP_ERROR;
   }
   if( seeedpool == NULL )
   {
      SCIPerrorMessage("Could not find Seeedpool!\n");
      return SCIP_ERROR;
   }

   /* write file */
   writeGpHeader( filename, outputname );
   writeGpSeeed( filename, seeed, seeedpool );

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */


/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeGp)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIPfreeMemory(scip, &readerdata);

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGp)
{
   MiscVisualization misc = new MiscVisualization();
   SeeedPtr seeed;
   char* filename;
   char* outputname;

   assert(scip != NULL);
   assert(file != NULL);

   /*@todo get seeed */

   filename = misc.GCGgetFilePath( scip, file );

   /* get filename for compiled file */
   outputname = misc.GCGgetVisualizationFilename( scip, seeed, "pdf" );

   GCGwriteGpVisualization( scip, filename, outputname, seeed->getID() );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/*
 * reader include
 */


/** includes the gp file reader into SCIP */
SCIP_RETCODE SCIPincludeReaderGp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include gp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
      NULL, readerFreeGp, NULL, readerWriteGp, NULL) );

   return SCIP_OKAY;
}

} /* namespace gcg */