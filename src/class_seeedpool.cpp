/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2015 Operations Research, RWTH Aachen University       */
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

/**@file   class_seeedpool.cpp
 * @brief  class with functions for seeedpoolnTotalSeeeds
 * @author Michael Bastubbe
 * @author Julius Hense
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*@todo don't disable lint */
/*lint -e64 disable useless and wrong lint warning */

/*@todo this looks like a workaround, is disabling warnings a good idea? */
#ifdef __INTEL_COMPILER
#ifndef _OPENMP
#pragma warning disable 3180  /* disable wrong and useless omp warnings */
#endif
#endif

//#define WRITE_ORIG_CONSTYPES
//#define SCIP_DEBUG

#include "gcg.h"
#include "objscip/objscip.h"
#include "scip/scip.h"
#include "class_seeedpool.h"
#include "struct_detector.h"
#include "pub_decomp.h"
#include "struct_decomp.h"
#include "cons_decomp.h"
#include "decomp.h"
#include "scip_misc.h"
#include "scip/clock.h"
#include "scip/cons.h"
#include "scip/scip.h"
#include <algorithm>
#include <list>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <iomanip>
#include <queue>
#include <fstream>
#include <exception>



#ifdef _OPENMP
#include <omp.h>
#endif

#define SCIP_CALL_EXC( x ) do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( ( _restat_ = ( x ) ) != SCIP_OKAY )                                             \
                          {                                                                                   \
                             SCIPerrorMessage( "Error <%d> in function call\n", _restat_);                    \
                             throw std::exception();                                                          \
                          }                                                                                   \
                       }                                                                                      \
                       while( false )

#define ENUM_TO_STRING( x ) # x
#define DEFAULT_THREADS    0     /**< number of threads (0 is OpenMP default) */


//#ifdef WITH_PRINTORIGCONSTYPES
/** constraint type */
enum SCIP_Constype_orig
{
   SCIP_CONSTYPE_EMPTY         =  0,         /**<  */
   SCIP_CONSTYPE_FREE          =  1,         /**<  */
   SCIP_CONSTYPE_SINGLETON     =  2,         /**<  */
   SCIP_CONSTYPE_AGGREGATION   =  3,         /**<  */
   SCIP_CONSTYPE_VARBOUND      =  4,         /**<  */
   SCIP_CONSTYPE_SETPARTITION  =  5,         /**<  */
   SCIP_CONSTYPE_SETPACKING    =  6,         /**<  */
   SCIP_CONSTYPE_SETCOVERING   =  7,         /**<  */
   SCIP_CONSTYPE_CARDINALITY   =  8,         /**<  */
   SCIP_CONSTYPE_INVKNAPSACK   =  9,         /**<  */
   SCIP_CONSTYPE_EQKNAPSACK    = 10,         /**<  */
   SCIP_CONSTYPE_BINPACKING    = 11,         /**<  */
   SCIP_CONSTYPE_KNAPSACK      = 12,         /**<  */
   SCIP_CONSTYPE_INTKNAPSACK   = 13,         /**<  */
   SCIP_CONSTYPE_MIXEDBINARY   = 14,         /**<  */
   SCIP_CONSTYPE_GENERAL       = 15          /**<  */
};
typedef enum SCIP_Constype_orig SCIP_CONSTYPE_ORIG;
//#endif




namespace gcg{

/** local methods */

struct sort_decr
{
   bool operator()(
      const std::pair<int, int> &left,
      const std::pair<int, int> &right)
   {
      return left.second > right.second;
   }
};

struct sort_pred
{
   bool operator()(
      const std::pair<int, int> &left,
      const std::pair<int, int> &right)
   {
      return left.second < right.second;
   }
};

/** is constraint ranged row, i.e., -inf < lhs < rhs < inf? */
static
SCIP_Bool isRangedRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             lhs,
   SCIP_Real             rhs
   )
{
   assert(scip != NULL);

   return !(SCIPisEQ(scip, lhs, rhs)
      || SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs) );
}

/** is constraint ranged row, i.e., -inf < lhs < rhs < inf? */
static
SCIP_Bool isFiniteNonnegativeIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             x                   /**< value */
   )
{
   assert(scip != NULL);

   return (!SCIPisInfinity(scip, x) && !SCIPisNegative(scip, x) && SCIPisIntegral(scip, x));
}





/** returns a folder name for a seeed */
std::string getSeeedFolderLatex(
   SeeedPtr seeed
   )
{
   std::stringstream decompfilename;
   decompfilename << "dec" << seeed->getID() << ".pdf";

   return decompfilename.str();
}

/** returns true if there exists an unfinished child in childsfinished array */
SCIP_Bool unfinishedChildExists(
   std::vector<SCIP_Bool> const& childsfinished
   )
{
   for( size_t s = 0; s < childsfinished.size(); ++ s )
   {
      if( ! childsfinished[s] )
         return true;
   }
   return false;
}

/** returns first unfinished child in childfinished array (-1 if there is none) */
int getFirstUnfinishedChild(
   std::vector<SCIP_Bool> const& childsfinished,
   std::vector<int> const& childs
   )
{
   for( size_t s = 0; s < childsfinished.size(); ++ s )
   {
      if( ! childsfinished[s] )
         return childs[s];
   }
   return - 1;
}

/** returns index of first unfinished child in childfinished array (-1 if there is none) */
int getFirstUnfinishedChildId(
   std::vector<SCIP_Bool> const& childsfinished,
   std::vector<int> const& childs
   )
{
   for( size_t s = 0; s < childsfinished.size(); ++ s )
   {
      if( ! childsfinished[s] )
         return (int) s;
   }
   return - 1;
}

/** sets next possible child finished (should equal child)
 *  returns true if next child is the last unfinished child */
SCIP_Bool finishNextChild(
   std::vector<int>& childs,
   std::vector<SCIP_Bool>& childsfinished,
   int child
   )
{
   for( size_t s = 0; s < childsfinished.size(); ++ s )
   {
      if( ! childsfinished[s] )
      {
         assert( childs[s] == child );
         childsfinished[s] = true;
         return s == childsfinished.size() - 1;
      }
   }
   return false;
}

/** returns the detector chain info of a seeed in a latex format */
std::string writeSeeedDetectorChainInfoLatex(
   SeeedPtr seeed,
   int currheight,
   int visucounter
   )
{
   std::stringstream line;
   std::string relposition;
   int position = visucounter % 3;
   if( position == 0 )
      relposition = "above";
   else if( position == 1 )
      relposition = "";
   else if( position == 2 )
      relposition = "below";
   else
      relposition = "below left";

   if( currheight != 1 )
      relposition = "";

   if( currheight > seeed->getNDetectorchainInfo() )
      line << "edge from parent node [" << relposition << "] {no info" << seeed->getID() << "-" << currheight - 1 << " } ";
   else
   {
      std::string oldinfo = seeed->getDetectorchainInfo( currheight - 1 );
      /* take latexified detctorchaininfo */
      size_t index = 0;
      while( true )
      {
         /* Locate the substring to replace. */
         index = oldinfo.find( "_", index );
         if( index == std::string::npos )
            break;
         if( index > 0 && oldinfo.at( index - 1 ) == '\\' )
         {
            ++ index;
            continue;
         }

         /* Make the replacement. */
         oldinfo.replace( index, 1, "\\_" );

         /* Advance index forward so the next iteration doesn't pick it up as well. */
         index += 2;
      }

 //     std::cout << "oldinfo: " << oldinfo << std::endl;
      line << "edge from parent node [" << relposition << "] {" << oldinfo << "} ";
   }
   return line.str();
}

/** returns the seeed info of a seeed in a latex format */
std::string writeSeeedInfoLatex(
   SeeedPtr seeed
   )
{
   std::stringstream line;
   line << "\\node[below = \\belowcaptionskip of s" << seeed->getID() << "] (caps" << seeed->getID() << ") {\\scriptsize "
      << seeed->getShortCaption() << "}; " << std::endl;

   return line.str();
}

/** returns the include graphics line for a seeed in a latex format */
std::string writeSeeedIncludeLatex(
   SeeedPtr seeed,
   std::string workfolder
   )
{
   std::stringstream line;
   line << " (s" << seeed->getID() << ") { \\includegraphics[width=0.15\\textwidth]{" << getSeeedFolderLatex( seeed )
      << "} }" << std::endl;

   return line.str();
}

/** writes detector call round information to passed parameter */
SCIP_RETCODE getDetectorCallRoundInfo(
   SCIP* scip,
   const char* detectorname,
   SCIP_Bool transformed,
   int* maxcallround,
   int* mincallround,
   int* freqcallround
   )
{
   char setstr[SCIP_MAXSTRLEN];
   if( transformed )
   {
      (void) SCIPsnprintf( setstr, SCIP_MAXSTRLEN, "detectors/%s/maxcallround", detectorname );
      SCIP_CALL( SCIPgetIntParam( scip, setstr, maxcallround ) );
      (void) SCIPsnprintf( setstr, SCIP_MAXSTRLEN, "detectors/%s/mincallround", detectorname );
      SCIP_CALL( SCIPgetIntParam( scip, setstr, mincallround ) );
      (void) SCIPsnprintf( setstr, SCIP_MAXSTRLEN, "detectors/%s/freqcallround", detectorname );
      SCIP_CALL_ABORT( SCIPgetIntParam( scip, setstr, freqcallround ) );
   }
   else
   {
      (void) SCIPsnprintf( setstr, SCIP_MAXSTRLEN, "detectors/%s/origmaxcallround", detectorname );
      SCIP_CALL( SCIPgetIntParam( scip, setstr, maxcallround ) );
      (void) SCIPsnprintf( setstr, SCIP_MAXSTRLEN, "detectors/%s/origmincallround", detectorname );
      SCIP_CALL( SCIPgetIntParam( scip, setstr, mincallround ) );
      (void) SCIPsnprintf( setstr, SCIP_MAXSTRLEN, "detectors/%s/origfreqcallround", detectorname );
      SCIP_CALL( SCIPgetIntParam( scip, setstr, freqcallround ) );
   }

   return SCIP_OKAY;
}

/** returns TRUE if seeed i has a greater MaxWhiteScore than seeed j */
SCIP_Bool cmpSeeedsMaxWhite(
   SeeedPtr i,
   SeeedPtr j
   )
{
   return ( i->getMaxWhiteScore() > j->getMaxWhiteScore() );
}

/** returns TRUE if seeed i has a greater border area score than seeed j */
SCIP_Bool cmpSeeedsBorderArea(
   SeeedPtr i,
   SeeedPtr j
   )
{
   return ( i->getScore( BORDER_AREA ) > j->getScore( BORDER_AREA ) );
}

/** returns TRUE if seeed i has a greater score than seeed j */
SCIP_Bool cmpSeeedsClassic(
   SeeedPtr i,
   SeeedPtr j
   )
{
   return ( i->getScore( CLASSIC )  > j->getScore( CLASSIC ) );
}

/** returns TRUE if seeed i has a greater score than seeed j */
SCIP_Bool cmpSeeedsFWhite(
   SeeedPtr i,
   SeeedPtr j
   )
{
   return ( i->getScore( MAX_FORESSEEING_WHITE )  > j->getScore( MAX_FORESSEEING_WHITE ) );
}

/** returns TRUE if seeed i has a greater score than seeed j */
SCIP_Bool cmpSeeedsAggFWhite(
   SeeedPtr i,
   SeeedPtr j
   )
{
   return ( i->getScore( MAX_FORESSEEING_AGG_WHITE )  > j->getScore( MAX_FORESSEEING_AGG_WHITE ) );
}


/** returns TRUE if seeed i has a greater score than seeed j */
SCIP_Bool cmpSeeedsPPCfWhite(
   SeeedPtr i,
   SeeedPtr j
   )
{
   return ( i->getScore( SETPART_FWHITE )  > j->getScore( SETPART_FWHITE ) );
}

/** returns TRUE if seeed i has a greater score than seeed j */
SCIP_Bool cmpSeeedsPPCaggFWhite(
   SeeedPtr i,
   SeeedPtr j
   )
{
   return ( i->getScore( SETPART_AGG_FWHITE )  > j->getScore( SETPART_AGG_FWHITE ) );
}




/* method to thin out the vector of given seeeds */
std::vector<SeeedPtr> thinout(
   std::vector<SeeedPtr> finishedseeeds,
   size_t ndecomps,
   SCIP_Bool addtrivialdecomp
   )
{
   std::vector<SeeedPtr> justbest( 0 );
   for( size_t dec = 0; dec < ndecomps && dec < finishedseeeds.size(); ++ dec )
   {
      justbest.push_back( finishedseeeds[dec] );
   }

   if( addtrivialdecomp )
   {
      for( size_t dec = 0; dec < finishedseeeds.size(); ++ dec )
      {
         if( finishedseeeds[dec]->getNMasterconss() == 0 && finishedseeeds[dec]->getNLinkingvars() == 0
            && finishedseeeds[dec]->getNBlocks() == 1 )
         {
            justbest.push_back( finishedseeeds[dec] );
         }
      }
   }
   return justbest;
}

/** returns levenshtein distance between two strings */
int calcLevenshteinDistance(
   std::string s,
   std::string t
   )
{
   // trivial cases
   if( s.compare( t ) == 0 )
      return 0;
   if( s.length() == 0 )
      return t.length();
   if( t.length() == 0 )
      return s.length();

   // create two work vectors of integer distances
   std::vector<int> v0( t.length() + 1 );
   std::vector<int> v1( t.length() + 1 );

   /* initialize v0 (the previous row of distances)
    * this row is A[0][i]: edit distance for an empty s
    * the distance is just the number of characters to delete from t */
   for( size_t i = 0; i < v0.size(); i ++ )
   {
      v0[i] = i;
   }

   for( size_t i = 0; i < s.length(); i ++ )
   {
      // calculate v1 (current row distances) from the previous row v0

      /* first element of v1 is A[i+1][0]
       * edit distance is delete (i+1) chars from s to match empty t */
      v1[0] = i + 1;

      // use formula to fill in the rest of the row
      for( size_t j = 0; j < t.length(); j ++ )
      {
         int cost = ( s[i] == t[j] ) ? 0 : 1;
         v1[j + 1] = std::min( v1[j] + 1, std::min( v0[j + 1] + 1, v0[j] + cost ) );
      }

      // copy v1 (current row) to v0 (previous row) for next iteration
      for( size_t j = 0; j < v0.size(); j ++ )
         v0[j] = v1[j];
   }

   return v1[t.length()];
}

/** removes all digits from string str */
void removeDigits(
   char *str,
   int *nremoved
   )
{
   char digits[11] = "0123456789";
   * nremoved = 0;

   for( int i = 0; i < 10; ++ i )
   {
      char digit = digits[i];
      size_t j = 0;
      while( j < strlen( str ) )
      {
         if( str[j] == digit )
         {
            * nremoved = * nremoved + 1;
            for( size_t k = j; k < strlen( str ); ++ k )
            {
               str[k] = str[k + 1];
            }
         }
         else
            ++ j;
      }
   }
}

/** method to calculate the greatest common divisor */
int gcd(
   int a,
   int b
   )
{
   return b == 0 ? a : gcd( b, a % b );
}

/** returns the relevant representative of a cons */
SCIP_CONS* consGetRelevantRepr(
   SCIP* scip,
   SCIP_CONS* cons
   )
{
   return cons;
}

/** returns the relevant representative of a var */
SCIP_VAR* varGetRelevantRepr(
   SCIP* scip,
   SCIP_VAR* var
   )
{
   return SCIPvarGetProbvar( var );
}





/** returns FALSE if there exists a seeed in seeeds that is a duplicate of compseeed */
SCIP_Bool seeedIsNoDuplicateOfSeeeds(
   SeeedPtr compseeed,
   std::vector<SeeedPtr> const & seeeds,
   bool sort
   )
{
   assert( compseeed != NULL );
   SCIP_Bool isduplicate;

   for( size_t i = 0; i < seeeds.size(); ++ i )
   {
      assert( seeeds[i] != NULL );

      compseeed->isEqual( seeeds[i], & isduplicate, sort );
      if( isduplicate )
         return false;
   }
   return true;
}


/** returns FALSE if there exists a seed in currSeeeds or finishedSeeeds that is a duplicate of seeed */
SCIP_Bool seeedIsNoDuplicate(
   SeeedPtr seeed,
   std::vector<SeeedPtr> const & currseeeds,
   std::vector<SeeedPtr> const & finishedseeeds,
   bool sort
   )
{
   SCIP_Bool bool1 = seeedIsNoDuplicateOfSeeeds( seeed, currseeeds, sort );
   SCIP_Bool bool2 = seeedIsNoDuplicateOfSeeeds( seeed, finishedseeeds, sort );
   return ( bool1 && bool2 );
}

/** constructor */
Seeedpool::Seeedpool(
   SCIP* givenScip,
   const char* conshdlrName,
   SCIP_Bool _transformed
   ) :
   scip( givenScip ), incompleteSeeeds( 0 ), currSeeeds( 0 ), ancestorseeeds( 0 ),
   nVars( SCIPgetNVars( givenScip ) ), nConss( SCIPgetNConss( givenScip ) ), nDetectors( 0 ),
   nFinishingDetectors( 0 ), nPostprocessingDetectors(0), nnonzeros( 0 ), candidatesNBlocks( 0 ), transformed( _transformed )
{
   SCIP_CONS** conss;
   SCIP_VAR** vars;

   int ndetectors;
   DEC_Detector** detectors;
   SCIP_Bool createconssadj;
   SCIP_Bool useconnected;
   SCIP_Bool useconssadj;

   createconssadj = TRUE;

   if( transformed )
      SCIPgetBoolParam(scip, "detectors/connectedbase/enabled", &useconnected);
   else
      SCIPgetBoolParam(scip, "detectors/connectedbase/origenabled", &useconnected);

   SCIPgetBoolParam(scip, "detectors/connectedbase/useconssadj", &useconssadj);

  // createconssadj = useconnected && useconssadj;

   if( ! transformed )
   {
      nVars = SCIPgetNOrigVars( scip );
      nConss = SCIPgetNOrigConss( scip );
   }

   int relevantVarCounter = 0;
   int relevantConsCounter = 0;

   detectors = SCIPconshdlrDecompGetDetectors(scip);
   ndetectors = SCIPconshdlrDecompGetNDetectors(scip);

   /** set detection data */
   SCIP_CALL_ABORT( SCIPgetIntParam( givenScip, "detection/maxrounds", & maxndetectionrounds ) );

   assert( ndetectors > 0 );


   /** set up enabled detectors and store them */
   for( int d = 0; d < ndetectors; ++d )
   {
      DEC_DETECTOR* detector;
      detector = detectors[d];

      assert( detector != NULL );
      if( transformed )
      {
         if( ! detector->enabled || detector->propagateSeeed == NULL )
            continue;
      }
      else
      {
         if( ! detector->enabledOrig || detector->propagateSeeed == NULL )
            continue;
      }

      scipDetectorToIndex[detector] = nDetectors;
      detectorToScipDetector.push_back( detector );
      ++ nDetectors;
   }

   /** set up enabled finishing detectors */
   for( int d = 0; d < ndetectors; ++d )
   {
      DEC_DETECTOR* detector;

      detector = detectors[d];
      assert( detector != NULL );
      if( ! detector->enabledFinishing || detector->finishSeeed == NULL )
         continue;

      scipFinishingDetectorToIndex[detector] = nFinishingDetectors;
      detectorToFinishingScipDetector.push_back( detector );
      ++ nFinishingDetectors;
   }

   /** set up enabled postprocessing detectors */
   for( int d = 0; d < ndetectors; ++d )
   {
      DEC_DETECTOR* detector;

      detector = detectors[d];
      assert( detector != NULL );
      if( ! detector->enabledPostprocessing || detector->postprocessSeeed == NULL )
         continue;

      scipPostprocessingDetectorToIndex[detector] = nPostprocessingDetectors;
      detectorToPostprocessingScipDetector.push_back( detector );
      ++nPostprocessingDetectors;
   }


   /** initilize matrix datastructures */
   if( transformed )
   {
      conss = SCIPgetConss( scip );
      vars = SCIPgetVars( scip );
   }
   else
   {
      conss = SCIPgetOrigConss( scip );
      vars = SCIPgetOrigVars( scip );
   }

   /** assign an index to every cons and var
    * @TODO: are all constraints/variables relevant? (probvars etc)  */
   for( int i = 0; i < nConss; ++ i )
   {
      SCIP_CONS* relevantCons;

      relevantCons = transformed ? consGetRelevantRepr( scip, conss[i] ) : conss[i];

      if( SCIPconsIsDeleted( relevantCons ) || SCIPconsIsObsolete(relevantCons) )
         continue;

      if( relevantCons != NULL )
      {
         scipConsToIndex[relevantCons] = relevantConsCounter;
         consToScipCons.push_back( relevantCons );

        //SCIPcaptureCons(scip, relevantCons);

//         if( relevantConsCounter == 7712 && transformed )
//         {
//            ++ relevantConsCounter;
//            -- relevantConsCounter;
//         }
         ++relevantConsCounter;
      }
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "relevant cons is NULL\n");
      }
   }

   for( int i = 0; i < nVars; ++ i )
   {
      SCIP_VAR* relevantVar;

      if( transformed )
         relevantVar = varGetRelevantRepr( scip, vars[i] );
      else
         relevantVar = vars[i];

      if( relevantVar != NULL )
      {
         scipVarToIndex[relevantVar] = relevantVarCounter;
         varToScipVar.push_back( relevantVar );
         ++ relevantVarCounter;
      }
   }

   /** from here on nVars and nConss represents the relevant numbers */
   nVars = relevantVarCounter;
   nConss = relevantConsCounter;
   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, " nvars: %d / nconss: %d \n", nVars, nConss  );
   varsForConss = std::vector<std::vector<int>>( nConss );
   valsForConss = std::vector < std::vector < SCIP_Real >> ( nConss );
   conssForVars = std::vector<std::vector<int>>( nVars );

   assert( (int) varToScipVar.size() == nVars );
   assert( (int) consToScipCons.size() == nConss );

   /** assumption: now every relevant constraint and variable has its index
    * and is stored in the corresponding unordered_map */
   /** find constraint <-> variable relationships and store them in both directions */
   for( int i = 0; i < (int) consToScipCons.size(); ++ i )
   {
      SCIP_CONS* cons;
      SCIP_VAR** currVars;
      SCIP_Real* currVals;
      int nCurrVars;

      cons = consToScipCons[i];

      nCurrVars = GCGconsGetNVars( scip, cons );

      if( nCurrVars == 0 )
         continue;

      assert(SCIPconsGetName( cons) != NULL);

      SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & currVars, nCurrVars ) );
      SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & currVals, nCurrVars ) );
      SCIP_CALL_ABORT( GCGconsGetVars( scip, cons, currVars, nCurrVars ) );
      SCIP_CALL_ABORT( GCGconsGetVals( scip, cons, currVals, nCurrVars ) );

      for( int currVar = 0; currVar < nCurrVars; ++ currVar )
      {
         int varIndex;
         std::tr1::unordered_map<SCIP_VAR*, int>::const_iterator iterVar;

         /*@todo remove this after the bug is fixed */
         /* because of the bug of GCGconsGet*()-methods some variables have to be negated */
         if( ! SCIPvarIsNegated( currVars[currVar] ) )
            iterVar = scipVarToIndex.find( currVars[currVar] );
         else
            iterVar = scipVarToIndex.find( SCIPvarGetNegatedVar( currVars[currVar] ) );

         if( iterVar == scipVarToIndex.end() )
            continue;

         varIndex = iterVar->second;

         varsForConss[i].push_back( varIndex );
         conssForVars[varIndex].push_back( i );
         valsForConss[i].push_back( currVals[currVar] );
         valsMap[std::pair<int, int>( i, varIndex )] = currVals[currVar];
         ++ nnonzeros;
      }
      SCIPfreeBufferArrayNull( scip, & currVals );
      SCIPfreeBufferArrayNull( scip, & currVars );

   }

   if( createconssadj )
   {
      std::vector<std::list<int>> conssadjacenciestemp( consToScipCons.size(), std::list<int>(0) );

      /** find constraint <-> constraint relationships and store them in both directions */
      for( size_t i = 0; i < consToScipCons.size(); ++i )
      {
         for( size_t varid = 0; varid < varsForConss[i].size(); ++varid )
         {
            int var = varsForConss[i][varid];

            for( size_t otherconsid = 0; otherconsid < conssForVars[var].size(); ++otherconsid )
            {
               int othercons = conssForVars[var][otherconsid];
               if( othercons == (int) i )
                  continue;

               std::list<int>::iterator consiter = std::lower_bound( conssadjacenciestemp[i].begin(),conssadjacenciestemp[i].end(), othercons);

               if( consiter == conssadjacenciestemp[i].end() || *consiter != othercons )
                  conssadjacenciestemp[i].insert(consiter, othercons);
            }
         }
      }

      for( size_t i = 0; i < consToScipCons.size(); ++ i )
      {
         conssadjacencies.push_back(std::vector<int>(0));
         std::list<int>::iterator consiter = conssadjacenciestemp[i].begin();
         std::list<int>::iterator consiterend = conssadjacenciestemp[i].end();
         for( ; consiter != consiterend; ++consiter )
         {
            conssadjacencies[i].push_back(*consiter);
         }
      }
   }
   /*  init  seeedpool with empty seeed */
   SeeedPtr emptyseeed = new Seeed( scip, SCIPconshdlrDecompGetNextSeeedID( scip ), nConss, nVars );

   addSeeedToCurr( emptyseeed );
   addSeeedToAncestor(emptyseeed);

   for( int i = 0; i < SCIPconshdlrDecompGetNBlockNumberCandidates( scip ); ++i )
   {
      addUserCandidatesNBlocks(SCIPconshdlrDecompGetBlockNumberCandidate( scip, i ) );
   }

} //end constructor

/** destructor */
Seeedpool::~Seeedpool()
{


   for( size_t i = 0; i < ancestorseeeds.size(); ++i )
   {
      size_t help = ancestorseeeds.size() - i - 1;
      if( ancestorseeeds[help] != NULL && ancestorseeeds[help]->getID() >= 0 )
      {
         delete ancestorseeeds[help];
      }

   }


   for( size_t i = 0; i < finishedSeeeds.size(); ++i )
   {
      size_t help = finishedSeeeds.size() - i - 1;
      if( finishedSeeeds[help] != NULL && finishedSeeeds[help]->getID() >= 0 )
      {
         delete finishedSeeeds[help];
      }
   }

   for( size_t i = 0; i < incompleteSeeeds.size(); ++i )
   {
      size_t help = incompleteSeeeds.size() - i - 1;
      if( incompleteSeeeds[help] != NULL && incompleteSeeeds[help]->getID() >= 0 )
         delete incompleteSeeeds[help];
   }


   for( size_t i = 0; i < consclassescollection.size(); ++ i )
   {
      size_t help = consclassescollection.size() - i - 1;
      if( consclassescollection[help] != NULL )
         delete consclassescollection[help];
   }

   for( size_t i = 0; i < varclassescollection.size(); ++ i )
   {
      size_t help = varclassescollection.size() - i - 1;
      if( varclassescollection[help] != NULL )
         delete varclassescollection[help];
   }
}

/** creates constraint and variable classifiers, and deduces block number candidates */
SCIP_RETCODE Seeedpool::calcClassifierAndNBlockCandidates(
   SCIP* givenScip /**< SCIP data structure */
   )
{
   SCIP_Bool conssclassnnonzeros;
   SCIP_Bool conssclassscipconstypes;
   SCIP_Bool conssclassmiplibconstypes;
   SCIP_Bool conssclassconsnamenonumbers;
   SCIP_Bool conssclassconsnamelevenshtein;
   SCIP_Bool varclassscipvartypes;
   SCIP_Bool varclassobjvals;
   SCIP_Bool varclassobjvalsigns;

   if( transformed )
   {
      SCIPgetBoolParam( scip, "detection/consclassifier/nnonzeros/enabled", & conssclassnnonzeros );
      SCIPgetBoolParam( scip, "detection/consclassifier/scipconstype/enabled", & conssclassscipconstypes );
      SCIPgetBoolParam( scip, "detection/consclassifier/miplibconstype/enabled", & conssclassmiplibconstypes );
      SCIPgetBoolParam( scip, "detection/consclassifier/consnamenonumbers/enabled", & conssclassconsnamenonumbers );
      SCIPgetBoolParam( scip, "detection/consclassifier/consnamelevenshtein/enabled", & conssclassconsnamelevenshtein );
      SCIPgetBoolParam( scip, "detection/varclassifier/scipvartype/enabled", & varclassscipvartypes );
        SCIPgetBoolParam(scip, "detection/varclassifier/objectivevalues/enabled", &varclassobjvals);
        SCIPgetBoolParam(scip, "detection/varclassifier/objectivevaluesigns/enabled", &varclassobjvalsigns);
   }
   else
   {
      SCIPgetBoolParam( scip, "detection/consclassifier/nnonzeros/origenabled", & conssclassnnonzeros );
      SCIPgetBoolParam( scip, "detection/consclassifier/scipconstype/origenabled", & conssclassscipconstypes );
      SCIPgetBoolParam( scip, "detection/consclassifier/miplibconstype/origenabled", & conssclassmiplibconstypes );
      SCIPgetBoolParam( scip, "detection/consclassifier/consnamenonumbers/origenabled", & conssclassconsnamenonumbers );
      SCIPgetBoolParam( scip, "detection/consclassifier/consnamelevenshtein/origenabled", & conssclassconsnamelevenshtein );
      SCIPgetBoolParam( scip, "detection/varclassifier/scipvartype/origenabled", & varclassscipvartypes );
        SCIPgetBoolParam(scip, "detection/varclassifier/objectivevalues/origenabled", &varclassobjvals);
        SCIPgetBoolParam(scip, "detection/varclassifier/objectivevaluesigns/origenabled", &varclassobjvalsigns);
   }

   if( conssclassnnonzeros )
      addConsClassifier( createConsClassifierForNNonzeros() );
   if( conssclassscipconstypes )
      addConsClassifier( createConsClassifierForSCIPConstypes() );
   if( conssclassmiplibconstypes )
      addConsClassifier( createConsClassifierForMiplibConstypes() );

   if( conssclassconsnamenonumbers )
      addConsClassifier( createConsClassifierForConsnamesDigitFreeIdentical() );
   if( conssclassconsnamelevenshtein )
      addConsClassifier( createConsClassifierForConsnamesLevenshteinDistanceConnectivity( 1 ) );

   if( varclassscipvartypes )
      addVarClassifier( createVarClassifierForSCIPVartypes() );
     if ( varclassobjvals )
        addVarClassifier( createVarClassifierForObjValues() );
     if ( varclassobjvalsigns )
        addVarClassifier( createVarClassifierForObjValueSigns() );

   reduceConsclasses();
   reduceVarclasses();

   calcCandidatesNBlocks();

   return SCIP_OKAY;
}

/** constructs seeeds using the registered detectors
 *  @return user has to free seeeds */
std::vector<SeeedPtr> Seeedpool::findSeeeds()
{
   /** 1) read parameter, as there are: maxrounds
    *  2) loop rounds
    *  3) every seeed in seeeds
    *  4) every detector not registered yet propagates seeed
    *  5)  */

   bool displaySeeeds = false;
   int verboseLevel;
   std::vector<int> successDetectors;
   std::vector<SeeedPtr> delSeeeds;
   bool duplicate;

   successDetectors = std::vector<int>( nDetectors, 0 );

   delSeeeds = std::vector < SeeedPtr > ( 0 );

   verboseLevel = 1;

   /** set detection data */
   SCIP_CALL_ABORT( SCIPgetIntParam( scip, "detection/maxrounds", & maxndetectionrounds ) );


    /** @TODO this does not look well streamlined: currseeeds should be empty here, and seeedstopopulate should be the only seeeds to poopulate */
   for( size_t i = 0; i < currSeeeds.size(); ++ i )
   {
      SCIP_CALL_ABORT( prepareSeeed( currSeeeds[i] ) );
       if( currSeeeds[i]->getID() < 0 )
          currSeeeds[i]->setID(getNewIdForSeeed() );
   }

   /** add translated original seeeds (of unpresolved problem) */
   for( size_t i = 0; i < seeedstopopulate.size(); ++ i )
   {
      SCIP_CALL_ABORT( prepareSeeed( seeedstopopulate[i] ) );
      if( seeedIsNoDuplicateOfSeeeds( seeedstopopulate[i], currSeeeds, true ) )
         currSeeeds.push_back( seeedstopopulate[i] );
       else
          continue;
       if( seeedstopopulate[i]->getID() < 0 )
          seeedstopopulate[i]->setID(getNewIdForSeeed() );
   }

   for( int round = 0; round < maxndetectionrounds; ++ round )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Begin of detection round %d of %d total rounds \n", round, maxndetectionrounds);
      std::vector<SeeedPtr> nextSeeeds = std::vector < SeeedPtr > ( 0 );
      std::vector<SeeedPtr> currSeeedsToDelete = std::vector < SeeedPtr > ( 0 );

#pragma omp parallel for schedule( static, 1 )
      for( size_t s = 0; s < currSeeeds.size(); ++ s )
      {
         SeeedPtr seeedPtr;
         seeedPtr = currSeeeds[s];

#pragma omp critical ( ostream )
         {
            if( displaySeeeds || verboseLevel >= 1 )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Start to propagate seeed with id %d (%d of %d in round %d) \n", seeedPtr->getID(), s, currSeeeds.size(), round);
               if( displaySeeeds )
                  seeedPtr->displaySeeed();
            }
         }

         /** the current seeed is handled by all detectors */
         for( int d = 0; d < nDetectors; ++ d )
         {
            SEEED_PROPAGATION_DATA* seeedPropData;
            seeedPropData = new SEEED_PROPAGATION_DATA();
            seeedPropData->seeedpool = this;
            seeedPropData->nNewSeeeds = 0;
            DEC_DETECTOR* detector;
            std::vector<SeeedPtr>::const_iterator newSIter;
            std::vector<SeeedPtr>::const_iterator newSIterEnd;
            int maxcallround;
            int mincallround;
            int freqcallround;
            const char* detectorname;
            SCIP_CLOCK* detectorclock;

            detector = detectorToScipDetector[d];
            detectorname = DECdetectorGetName( detector );
            SCIP_RESULT result = SCIP_DIDNOTFIND;

            /** if the seeed is also propagated by the detector go on with the next detector */
            if( seeedPtr->isPropagatedBy( detector ) && ! detector->usefulRecall )
               continue;

            /** check if detector is callable in current detection round */
            SCIP_CALL_ABORT(
               getDetectorCallRoundInfo( scip, detectorname, transformed, & maxcallround, & mincallround,
                  & freqcallround ) );

            if( maxcallround < round || mincallround > round || ( ( round - mincallround ) % freqcallround != 0 ) )
               continue;

#pragma omp critical ( seeedcount )
            seeedPropData->seeedToPropagate = new gcg::Seeed( seeedPtr );

            /** new seeeds are created by the current detector */
#pragma omp critical ( clock )
            SCIPcreateClock( scip, & detectorclock );
            SCIP_CALL_ABORT( SCIPstartClock( scip, detectorclock ) );

            if( verboseLevel >= 1 )
            {
#pragma omp critical ( scipinfo )
               SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "Detector %s started to propagate seeed with id %d )  \n",
                  DECdetectorGetName( detectorToScipDetector[d] ), seeedPtr->getID() );
            }

            SCIP_CALL_ABORT(
               detectorToScipDetector[d]->propagateSeeed( scip, detectorToScipDetector[d], seeedPropData, & result ) );

            for( int j = 0; j < seeedPropData->nNewSeeeds; ++ j )
            {
#pragma omp critical ( seeedcount )

               seeedPropData->newSeeeds[j]->setID( getNewIdForSeeed() );
               prepareSeeed( seeedPropData->newSeeeds[j] );
               assert( seeedPropData->newSeeeds[j]->checkConsistency( this ) );
               seeedPropData->newSeeeds[j]->addDecChangesFromAncestor( seeedPtr );
            }

            SCIP_CALL_ABORT( SCIPstopClock( scip, detectorclock ) );

#pragma omp critical ( clockcount )
            detectorToScipDetector[d]->dectime += SCIPgetClockTime( scip, detectorclock );

#pragma omp critical ( clock )
            SCIPfreeClock( scip, & detectorclock );

            if( seeedPropData->nNewSeeeds != 0 && ( displaySeeeds ) )
            {
#pragma omp critical ( ostream )
               SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "Detector %s found %d new seeed%s: \n",
                                 DECdetectorGetName( detectorToScipDetector[d] ), seeedPropData->nNewSeeeds, (seeedPropData->nNewSeeeds == 1 ? "": "s") );
#pragma omp critical ( ostream )
               SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "%d", seeedPropData->newSeeeds[0]->getID() );
               for( int j = 1; j < seeedPropData->nNewSeeeds; ++ j )
               {
#pragma omp critical ( ostream )
                  SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, ", %d", seeedPropData->newSeeeds[0]->getID() );
               }
#pragma omp critical ( ostream )
               SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "\n" );

               if( displaySeeeds )
               {
                  for( int j = 0; j < seeedPropData->nNewSeeeds; ++ j )
                  {
#pragma omp critical ( ostream )
                     seeedPropData->newSeeeds[j]->displaySeeed();
                  }

               }
            }
            else if( displaySeeeds )
            {
#pragma omp critical ( ostream )
               SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "Detector %s found 0 new seeeds.\n", DECdetectorGetName( detectorToScipDetector[d] ) );
            }

            /** if a new seeed is no duplicate it is either added to the nextRoundSeeeds or the finishedSeeeds  */
            for( int seeed = 0; seeed < seeedPropData->nNewSeeeds; ++ seeed )
            {
               SCIP_Bool noduplicate;
#pragma omp critical ( seeedptrstore )
               noduplicate = seeedIsNoDuplicate( seeedPropData->newSeeeds[seeed], nextSeeeds, finishedSeeeds, false );
               if( ! seeedPropData->newSeeeds[seeed]->isTrivial() && noduplicate )
               {
                  if( seeedPropData->newSeeeds[seeed]->getNOpenconss() == 0
                     && seeedPropData->newSeeeds[seeed]->getNOpenvars() == 0 )
                  {
                     if( verboseLevel > 2 )
                     {
#pragma omp critical ( ostream )
                        {
                           SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "Seeed %d is added to finished seeeds.\n", seeedPropData->newSeeeds[seeed]->getID() );
                           seeedPropData->newSeeeds[seeed]->showVisualisation( this );
                        }
                     }
#pragma omp critical ( seeedptrstore )
                     {
                        assert( seeedPropData->newSeeeds[seeed]->getID() >= 0 );
                        addSeeedToFinished( seeedPropData->newSeeeds[seeed], &noduplicate );
                     }
                  }
                  else
                  {
                     if( verboseLevel > 2 )
                     {
#pragma omp critical ( ostream )
                        {
                           SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "Seeed %d is addded to next round seeeds!\n", seeedPropData->newSeeeds[seeed]->getID() );
                           seeedPropData->newSeeeds[seeed]->showVisualisation( this );
                        }
                     }
#pragma omp critical ( seeedptrstore )
                     {
                        nextSeeeds.push_back( seeedPropData->newSeeeds[seeed] );
                     }
                  }
               }
               else
               {
                  delete seeedPropData->newSeeeds[seeed];
                  seeedPropData->newSeeeds[seeed] = NULL;
               }
            }
            /** cleanup propagation data structure */
            delete seeedPropData->seeedToPropagate;
            SCIPfreeMemoryArrayNull( scip, & seeedPropData->newSeeeds );
            seeedPropData->newSeeeds = NULL;
            seeedPropData->nNewSeeeds = 0;
            delete seeedPropData;
         } // end for detectors

         SCIPverbMessage( scip, SCIP_VERBLEVEL_HIGH, NULL, "Start finishing of partial decomposition %d.\n", seeedPtr->getID() );
         for( int d = 0; d < nFinishingDetectors; ++ d )
         {
            DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
            SCIP_RESULT result = SCIP_DIDNOTFIND;
            SEEED_PROPAGATION_DATA* seeedPropData;
            seeedPropData = new SEEED_PROPAGATION_DATA();
            seeedPropData->seeedpool = this;
            seeedPropData->nNewSeeeds = 0;
#pragma omp critical ( seeedcount )
            seeedPropData->seeedToPropagate = new gcg::Seeed( seeedPtr );

            if( verboseLevel > 2 )
#pragma omp critical ( ostream )
            {
               SCIPverbMessage( scip, SCIP_VERBLEVEL_FULL, NULL, "Check if finisher of detector %s  is enabled. \n", DECdetectorGetName( detectorToFinishingScipDetector[d] ) );
            }

            /** if the finishing of the detector is not enabled go on with the next detector */
            if( ! detector->enabledFinishing )
               continue;

            if( verboseLevel > 2 )
#pragma omp critical ( ostream )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Call finisher for detector %s\n",
                  DECdetectorGetName( detectorToFinishingScipDetector[d] ) );
            }
            SCIP_CALL_ABORT(
               detectorToFinishingScipDetector[d]->finishSeeed( scip, detectorToFinishingScipDetector[d], seeedPropData,
                  & result ) );

            for( int finished = 0; finished < seeedPropData->nNewSeeeds; ++ finished )
            {
               SeeedPtr seeed = seeedPropData->newSeeeds[finished];
#pragma omp critical ( seeedcount )
               seeed->setID( getNewIdForSeeed() );
               seeed->sort();
               seeed->calcHashvalue();
               seeed->setSeeedpool(this);
               seeed->addDecChangesFromAncestor( seeedPtr );
               seeed->setFinishedByFinisher( true );
#pragma omp critical ( seeedptrstore )
               {
                  SCIP_Bool success;
                  addSeeedToFinished( seeed, &success );
                  if( !success )
                  {
                     bool isIdentical = false;
                     for( size_t h = 0; h < finishedSeeeds.size(); ++ h )
                     {
                        if( seeed == finishedSeeeds[h] )
                        {
                           isIdentical = true;
                           break;
                        }
                     }

                     if( ! isIdentical )
                     {
                        currSeeedsToDelete.push_back( seeed );
                     }
                  }
               }
            }
            SCIPfreeMemoryArrayNull( scip, & seeedPropData->newSeeeds );
            delete seeedPropData->seeedToPropagate;
            seeedPropData->newSeeeds = NULL;
            seeedPropData->nNewSeeeds = 0;
            delete seeedPropData;
         }
          #pragma omp critical (seeedptrstore)
          addSeeedToAncestor(seeedPtr);
      } // end for currseeeds

      for( size_t s = 0; s < currSeeedsToDelete.size(); ++ s )
      {
         delete currSeeedsToDelete[s];
         currSeeedsToDelete[s] = NULL;
      }

      currSeeeds = nextSeeeds;
   } // end for rounds

   /** complete the currseeeds with finishing detectors and add them to finished seeeds */
#pragma omp parallel for schedule( static, 1 )
   for( size_t i = 0; i < currSeeeds.size(); ++ i )
   {
      SeeedPtr seeedPtr = currSeeeds[i];
      for( int d = 0; d < nFinishingDetectors; ++ d )
      {
         DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
         SCIP_RESULT result = SCIP_DIDNOTFIND;
         SEEED_PROPAGATION_DATA* seeedPropData;
         seeedPropData = new SEEED_PROPAGATION_DATA();
         seeedPropData->seeedpool = this;
         seeedPropData->nNewSeeeds = 0;

#pragma omp critical ( seeedcount )
         seeedPropData->seeedToPropagate = new gcg::Seeed( seeedPtr );

         if( verboseLevel > 2 )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "check if finisher of detector %s is enabled\n",
               DECdetectorGetName( detectorToFinishingScipDetector[d] ) );

         /** if the finishing of the detector is not enabled go on with the next detector */
         if( ! detector->enabledFinishing )
            continue;

         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "call finisher for detector %s \n ", DECdetectorGetName( detectorToFinishingScipDetector[d] ) );

         SCIP_CALL_ABORT(
            detectorToFinishingScipDetector[d]->finishSeeed( scip, detectorToFinishingScipDetector[d], seeedPropData,
               & result ) );

         for( int finished = 0; finished < seeedPropData->nNewSeeeds; ++ finished )
         {
            SeeedPtr seeed = seeedPropData->newSeeeds[finished];
#pragma omp critical ( seeedcount )
            seeed->setID( getNewIdForSeeed() );

            seeed->calcHashvalue();
            seeed->addDecChangesFromAncestor( seeedPtr );
            seeed->setFinishedByFinisher( true );
            seeed->setSeeedpool(this);

            if( seeedIsNoDuplicateOfSeeeds( seeed, finishedSeeeds, false ) )
            {
               if( verboseLevel > 2 )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, " Seeed %d is finished from next round seeeds!\n", seeed->getID() );
                  seeed->showVisualisation( this );
               }
#pragma omp critical ( seeedptrstore )
               {
                  SCIP_Bool success;
                  assert( seeed->getID() >= 0 );
                  addSeeedToFinished( seeed, &success  );
               }
            }
            else
               delete seeed;

            SCIPfreeMemoryArrayNull( scip, & seeedPropData->newSeeeds );
            seeedPropData->newSeeeds = NULL;
            seeedPropData->nNewSeeeds = 0;
         }

         delete seeedPropData->seeedToPropagate;
         delete seeedPropData;
      }
#pragma omp critical ( seeedptrstore )
       addSeeedToAncestor(seeedPtr);
   } // end for finishing curr seeeds

   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%d  finished seeeds are found.\n", (int) finishedSeeeds.size() );

   if( displaySeeeds )
   {
      for( size_t i = 0; i < finishedSeeeds.size(); ++ i )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%d-th finished seeed:\n", (int) i );
         finishedSeeeds[i]->displaySeeed();
      }
   }

   /** count the successful refinement calls for each detector */
   for( size_t i = 0; i < finishedSeeeds.size(); ++ i )
   {
      assert( finishedSeeeds[i]->checkConsistency( this ) );
      assert( finishedSeeeds[i]->getNOpenconss() == 0 );
      assert( finishedSeeeds[i]->getNOpenvars() == 0 );

      SCIP_CALL_ABORT( finishedSeeeds[i]->buildDecChainString() );
      for( int d = 0; d < nDetectors; ++ d )
      {
         if( finishedSeeeds[i]->isPropagatedBy( detectorToScipDetector[d] ) )
            successDetectors[d] += 1;
      }
   }

   /** preliminary output detector stats */

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Measured running time per detector: \n" );

   for( int i = 0; i < nDetectors; ++ i )
   {
      if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_HIGH )
      {
         std::cout << "Detector " << std::setw( 25 ) << std::setiosflags( std::ios::left )
         << DECdetectorGetName( detectorToScipDetector[i] ) << " \t worked on \t " << successDetectors[i] << " of "
         << finishedSeeeds.size() << "\t and took a total time of \t" << detectorToScipDetector[i]->dectime << std::endl;
      }
   }

   if( (int) finishedSeeeds.size() != 0 )
   {
       SCIP_Real maxscore = finishedSeeeds[0]->getScore( SCIPconshdlrDecompGetCurrScoretype( scip ) );
      //            SeeedPtr bestSeeed = finishedSeeeds[0];
      for( size_t i = 1; i < finishedSeeeds.size(); ++ i )
      {
          SCIP_Real score = finishedSeeeds[i]->getScore( SCIPconshdlrDecompGetCurrScoretype( scip ) );
         if( score > maxscore )
         {
            maxscore = score;
         }
      }
   }

   /** delete the seeeds */
   for( size_t c = 0; c < currSeeeds.size(); ++c )
   {
      duplicate = false;
      SCIP_CALL_ABORT( currSeeeds[c]->buildDecChainString() );
      for( size_t d = 0; d < delSeeeds.size(); ++ d )
      {
         if( currSeeeds[c] == delSeeeds[d] )
         {
            duplicate = true;
            break;
         }
      }
      if( ! duplicate )
      {
         delSeeeds.push_back( currSeeeds[c] );
      }
   }

   /** postprocess the finished seeeds */
   std::vector<SeeedPtr >postprocessed(0);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Started Postprocessing of decompositions...\n");
#pragma omp parallel for schedule( static, 1 )
   for( size_t c = 0 ; c < finishedSeeeds.size(); ++c )
   {
      SeeedPtr seeedPtr = finishedSeeeds[c];
      for( int d = 0; d < getNPostprocessingDetectors(); ++d )
      {
         DEC_DETECTOR* detector = detectorToPostprocessingScipDetector[d];
         SCIP_RESULT result = SCIP_DIDNOTFIND;
         SEEED_PROPAGATION_DATA* seeedPropData;
         seeedPropData = new SEEED_PROPAGATION_DATA();
         seeedPropData->seeedpool = this;
         seeedPropData->nNewSeeeds = 0;

#pragma omp critical ( seeedcount )
         seeedPropData->seeedToPropagate = new gcg::Seeed( seeedPtr );

         if( verboseLevel > 2 )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "check if postprocessor of detector %s is enabled\n",
               DECdetectorGetName( detectorToPostprocessingScipDetector[d] ) );

         /** if the postprocessing of the detector is not enabled go on with the next detector */
         if( ! detector->enabledPostprocessing )
            continue;


         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "call finisher for detector %s \n ", DECdetectorGetName( detectorToPostprocessingScipDetector[d] ) );

         SCIP_CALL_ABORT(
            detectorToPostprocessingScipDetector[d]->postprocessSeeed( scip, detectorToPostprocessingScipDetector[d], seeedPropData,
               & result ) );

         for( int finished = 0; finished < seeedPropData->nNewSeeeds; ++ finished )
         {
            SeeedPtr seeed = seeedPropData->newSeeeds[finished];
#pragma omp critical ( seeedcount )
            seeed->setID( getNewIdForSeeed() );

            seeed->calcHashvalue();
            seeed->addDecChangesFromAncestor( seeedPtr );
            seeed->setFinishedByFinisher( true );
            seeed->setSeeedpool(this);

            if( seeedIsNoDuplicateOfSeeeds( seeed, finishedSeeeds, false ) && seeedIsNoDuplicateOfSeeeds( seeed, postprocessed, false )  )
            {
               if( verboseLevel > 2 )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, " Seeed %d is finished from next round seeeds!\n", seeed->getID() );
                  seeed->showVisualisation( this );
               }
#pragma omp critical ( seeedptrstore )
               {
                  assert( seeed->getID() >= 0 );
                  postprocessed.push_back( seeed);
               }
            }
            else
               delete seeed;

            SCIPfreeMemoryArrayNull( scip, & seeedPropData->newSeeeds );
            seeedPropData->newSeeeds = NULL;
            seeedPropData->nNewSeeeds = 0;
         }

         delete seeedPropData->seeedToPropagate;
         delete seeedPropData;
      }
//#pragma omp critical ( seeedptrstore )
//       addSeeedToAncestor(seeedPtr);

   } // end for postprocessing finished seeeds
   for( size_t c = 0; c < postprocessed.size(); ++c)
      addSeeedToFinishedUnchecked(postprocessed[c]);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "...finished postprocessing of decompositions. Added %d new decomps. \n", postprocessed.size());


   sortAllRelevantSeeeds();
   return finishedSeeeds;
 }

/* sorts seeeds in finished seeeds data structure according to their score */
 void Seeedpool::sortFinishedForScore()
{
   if( SCIPconshdlrDecompGetCurrScoretype(scip) == scoretype::MAX_WHITE )
      std::sort(finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsMaxWhite);

   if( SCIPconshdlrDecompGetCurrScoretype(scip) == scoretype::BORDER_AREA )
      std::sort(finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsBorderArea);

   if( SCIPconshdlrDecompGetCurrScoretype(scip) == scoretype::CLASSIC )
      std::sort(finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsClassic);

   if( SCIPconshdlrDecompGetCurrScoretype(scip) == scoretype::MAX_FORESSEEING_WHITE )
      std::sort(finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsFWhite);

   if( SCIPconshdlrDecompGetCurrScoretype(scip) == scoretype::MAX_FORESSEEING_AGG_WHITE )
         std::sort(finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsAggFWhite);

   if( SCIPconshdlrDecompGetCurrScoretype(scip) == scoretype::SETPART_FWHITE )
      std::sort(finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsPPCfWhite);

   if( SCIPconshdlrDecompGetCurrScoretype(scip) == scoretype::SETPART_AGG_FWHITE )
      std::sort(finishedSeeeds.begin(), finishedSeeeds.end(), cmpSeeedsPPCaggFWhite);

}



/** method to complete a set of incomplete seeeds with the help of all included detectors that implement a finishing method
 *  @return set of completed decomposition */
std::vector<SeeedPtr> Seeedpool::finishIncompleteSeeeds(
   std::vector<SeeedPtr> incompleteseeeds
   )
{
   std::vector<SeeedPtr> finisheds( 0, NULL );
   int verboseLevel = 1;

#pragma omp parallel for schedule( static, 1 )
   for( size_t i = 0; i < incompleteseeeds.size(); ++ i )
   {
      SeeedPtr seeedPtr = incompleteseeeds[i];

      for( int d = 0; d < nFinishingDetectors; ++ d )
      {
         DEC_DETECTOR* detector = detectorToFinishingScipDetector[d];
         SCIP_RESULT result = SCIP_DIDNOTFIND;
         SEEED_PROPAGATION_DATA* seeedPropData;
         seeedPropData = new SEEED_PROPAGATION_DATA();
         seeedPropData->seeedpool = this;
         seeedPropData->nNewSeeeds = 0;

#pragma omp critical ( seeedcount )
         seeedPropData->seeedToPropagate = new gcg::Seeed( seeedPtr );

         if( verboseLevel > 2 )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "check if finisher of detector %s is enabled\n", DECdetectorGetName( detectorToScipDetector[d] ) );

         /** if the finishing of the detector is not enabled go on with the next detector */
         if( ! detector->enabledFinishing )
            continue;

         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "call finisher for detector %s\n", DECdetectorGetName( detectorToFinishingScipDetector[d] ) );

         SCIP_CALL_ABORT(
            detectorToFinishingScipDetector[d]->finishSeeed( scip, detectorToFinishingScipDetector[d], seeedPropData,
               & result ) );

         for( int finished = 0; finished < seeedPropData->nNewSeeeds; ++ finished )
         {
            SeeedPtr seeed = seeedPropData->newSeeeds[finished];
#pragma omp critical ( seeedcount )
            seeed->setID( getNewIdForSeeed() );

            seeed->calcHashvalue();
            seeed->addDecChangesFromAncestor( seeedPtr );
            seeed->setFinishedByFinisher( true );

            if( seeedIsNoDuplicateOfSeeeds( seeed, finishedSeeeds, false ) )
            {
               if( verboseLevel > 2 )
               {
                  SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "seeed %d is finished from next round seeeds! \n", seeed->getID() );
                  seeed->showVisualisation( this );
               }
#pragma omp critical ( seeedptrstore )
               {
                  assert( seeed->getID() >= 0 );
                  finisheds.push_back( seeed );
               }
            }

            SCIPfreeMemoryArrayNull( scip, & seeedPropData->newSeeeds );
            seeedPropData->newSeeeds = NULL;
            seeedPropData->nNewSeeeds = 0;
         }

         delete seeedPropData->seeedToPropagate;
         delete seeedPropData;
      }
   } // end for finishing curr seeeds
   return finisheds;
}

/** calls findSeeeds method and translates the resulting seeeds into decompositions */
void Seeedpool::findDecompositions()
{
   std::vector<int> successDetectors;

   successDetectors = std::vector<int>( nDetectors, 0 );

   finishedSeeeds = findSeeeds();

   /* sort the seeeds according to maximum white measure */
   sortFinishedForScore();

}

/*SCIP_RETCODE DECdecompCheckConsistency(DEC_DECOMP* decomp)
{
   int c;
   int b;
   int v;

   for( v = 0; v < SCIPgetNVars(scip); ++v )
   {
      assert(SCIPhashmapExists(DECdecompGetVartoblock(decomp), SCIPgetVars(scip)[v]));
   }
}*/

/** returns seeed with the corresponding id or NULL if there is no finished seeed with that id */
  gcg::Seeed* Seeedpool::findFinishedSeeedByID(
     int      seeedid
     )
  {
     for( size_t fs = 0; fs < finishedSeeeds.size(); ++fs )
     {
        if ( finishedSeeeds[fs]->getID() == seeedid )
           return finishedSeeeds[fs];
     }

     return NULL;
  }



/** adds a seeed to ancestor seeeds */
void Seeedpool::addSeeedToAncestor(
   SeeedPtr seeed
   )
{
   ancestorseeeds.push_back( seeed );
}

/** adds a seeed to current seeeds */
void Seeedpool::addSeeedToCurr(
   SeeedPtr seeed
   )
{
   if( seeedIsNoDuplicateOfSeeeds(seeed, currSeeeds, false) )
      currSeeeds.push_back( seeed );
}

/** adds a seeed to finished seeeds */
void Seeedpool::addSeeedToFinished(
   SeeedPtr seeed,
   SCIP_Bool* success
   )
{
   if( seeedIsNoDuplicateOfSeeeds(seeed, finishedSeeeds, false) )
   {
      finishedSeeeds.push_back( seeed );
      *success = TRUE;
   }
   else
   {
      *success = FALSE;
   }
   return;
}

/** adds a seeed to finished seeeds withou checking for duplicates, dev has to check this on his own*/
void Seeedpool::addSeeedToFinishedUnchecked(
   SeeedPtr seeed
   )
{
   finishedSeeeds.push_back( seeed );

   return;
}



/** adds a seeed to incomplete seeeds */
void Seeedpool::addSeeedToIncomplete(
   SeeedPtr seeed,
   SCIP_Bool* success
   )
{
   if( seeedIsNoDuplicateOfSeeeds(seeed, incompleteSeeeds, false) )
   {
      incompleteSeeeds.push_back( seeed );
      *success = TRUE;
   }
   *success = FALSE;
   return;

}

/** clears ancestor seeed data structure */
void Seeedpool::clearAncestorSeeeds()
{
   ancestorseeeds.clear();
}

/** clears current seeed data structure */
void Seeedpool::clearCurrentSeeeds()
{
   currSeeeds.clear();
}

/** clears finished seeed data structure */
void Seeedpool::clearFinishedSeeeds()
{
   finishedSeeeds.clear();
}

/** clears incomplete seeed data structure */
void Seeedpool::clearIncompleteSeeeds()
{
   incompleteSeeeds.clear();
}

/** returns a seeed from ancestor seeed data structure */
SeeedPtr Seeedpool::getAncestorSeeed(
   int seeedindex
   )
{
   assert( 0 <= seeedindex && seeedindex < (int) ancestorseeeds.size() );

   return ancestorseeeds[seeedindex];
}

/** returns a seeed from current (open) seeed data structure */
SeeedPtr Seeedpool::getCurrentSeeed(
   int seeedindex
   )
{
   assert( 0 <= seeedindex && seeedindex < (int) currSeeeds.size() );

   return currSeeeds[seeedindex];
}

/** returns a seeed from finished seeed data structure */
SeeedPtr Seeedpool::getFinishedSeeed(
   int seeedindex
   )
{
   assert( 0 <= seeedindex && seeedindex < (int) finishedSeeeds.size() );

   return finishedSeeeds[seeedindex];
}

/** returns a seeed from incomplete seeed data structure */
SeeedPtr Seeedpool::getIncompleteSeeed(
   int seeedindex
   )
{
   assert( 0 <= seeedindex && seeedindex < (int) incompleteSeeeds.size() );

   return incompleteSeeeds[seeedindex];
}

/** returns size of ancestor seeed data structure */
int Seeedpool::getNAncestorSeeeds()
{
   return ancestorseeeds.size();
}

/** returns size of current (open) seeed data structure */
int Seeedpool::getNCurrentSeeeds()
{
   return currSeeeds.size();
}

/** returns size of finished seeed data structure */
int Seeedpool::getNFinishedSeeeds()
{
   return finishedSeeeds.size();
}

/** returns size of incomplete seeed data structure */
int Seeedpool::getNIncompleteSeeeds()
{
   return incompleteSeeeds.size();
}

/** returns true if the given seeed is a duplicate of a seeed that is already contained in
 *  finished seeeds or current seeeds data structure */
bool Seeedpool::hasDuplicate(
   SeeedPtr seeed
   )
{
   assert( seeed != NULL );

   return !seeedIsNoDuplicate( seeed, currSeeeds, finishedSeeeds, true );
}

/** translates seeeds and classifiers if the index structure of the problem has changed, e.g. due to presolving */
void Seeedpool::translateSeeedData(
   Seeedpool* origpool,
   std::vector<Seeed*> origseeeds,
   std::vector<Seeed*>& newseeeds,
   std::vector<ConsClassifier*> otherconsclassifiers,
   std::vector<ConsClassifier*>& newconsclassifiers,
   std::vector<VarClassifier*> othervarclassifiers,
   std::vector<VarClassifier*>& newvarclassifiers
   )
{
   assert( newseeeds.empty() );
   assert( newconsclassifiers.empty() );
   assert( newvarclassifiers.empty() );

   std::vector<int> rowothertothis;
   std::vector<int> rowthistoother;
   std::vector<int> colothertothis;
   std::vector<int> colthistoother;
   std::vector<int> missingrowinthis;

   SCIP_Bool presolvingdisabled;
   int presolvingrounds;

   presolvingdisabled = FALSE;

   SCIPgetIntParam(scip, "presolving/maxrounds", &presolvingrounds);

   if ( presolvingrounds == 0)
      presolvingdisabled = TRUE;

   SCIPverbMessage( this->scip, SCIP_VERBLEVEL_HIGH, NULL, "started translate seeed method: presolving is %s \n", (presolvingdisabled ? "disabled, try short method." : "enabled, has to do long version. " ) );

   if( presolvingdisabled )
   {
      missingrowinthis = std::vector<int>(0);
      rowothertothis = std::vector<int>(0);
      rowthistoother = std::vector<int>(0);
      colothertothis = std::vector<int>(0);
      colthistoother = std::vector<int>(0);
      for( int i = 0; i < nConss ; ++i )
      {
         rowothertothis.push_back(i);
         rowthistoother.push_back(i);
      }
      for( int j = 0; j < nVars ; ++j )
      {
         colthistoother.push_back(j);
         colothertothis.push_back(j);
      }
   } else
      calcTranslationMapping( origpool, rowothertothis, rowthistoother, colothertothis, colthistoother, missingrowinthis );

   SCIPverbMessage( this->scip, SCIP_VERBLEVEL_HIGH, NULL,
      " calculated translation; number of missing constraints: %d; number of other seeeds: %d \n", missingrowinthis.size(),
      origseeeds.size() );

   newseeeds = getTranslatedSeeeds( origseeeds, rowothertothis, rowthistoother, colothertothis, colthistoother );
   newconsclassifiers = getTranslatedConsClassifiers( otherconsclassifiers, rowothertothis, rowthistoother );
   newvarclassifiers = getTranslatedVarClassifiers( othervarclassifiers, colothertothis, colthistoother );
}

/** translates seeeds if the index structure of the problem has changed, e.g. due to presolving */
void Seeedpool::translateSeeeds(
   Seeedpool* origpool,
   std::vector<Seeed*> origseeeds,
   std::vector<Seeed*>& newseeeds
   )
{
   assert( newseeeds.empty() );

   SCIP_Bool presolvingdisabled;
   int presolvingrounds;

   std::vector<int> rowothertothis( 0 );
   std::vector<int> rowthistoother( 0 );
   std::vector<int> colothertothis( 0 );
   std::vector<int> colthistoother( 0 );
   std::vector<int> missingrowinthis( 0 );


   presolvingdisabled = FALSE;

   SCIPgetIntParam(scip, "presolving/maxrounds", &presolvingrounds);

   if ( presolvingrounds == 0)
      presolvingdisabled = TRUE;


   SCIPverbMessage( this->scip, SCIP_VERBLEVEL_HIGH, NULL, "started translate seeed method: presolving is %s \n", (presolvingdisabled ? "disabled, try short method." : "enabled, has to do long version. " ) );

   if( presolvingdisabled )
   {
      missingrowinthis = std::vector<int>(0);
      rowothertothis = std::vector<int>(0);
      rowthistoother = std::vector<int>(0);
      colothertothis = std::vector<int>(0);
      colthistoother = std::vector<int>(0);
      for( int i = 0; i < nConss ; ++i )
      {
         rowothertothis.push_back(i);
         rowthistoother.push_back(i);
      }
      for( int j = 0; j < nVars ; ++j )
      {
         colthistoother.push_back(j);
         colothertothis.push_back(j);
      }
   } else
      calcTranslationMapping( origpool, rowothertothis, rowthistoother, colothertothis, colthistoother, missingrowinthis );

   SCIPverbMessage( this->scip, SCIP_VERBLEVEL_HIGH, NULL,
      " calculated translation; number of missing constraints: %d; number of other seeeds: %d \n", missingrowinthis.size(),
      origseeeds.size() );

   newseeeds = getTranslatedSeeeds( origseeeds, rowothertothis, rowthistoother, colothertothis, colthistoother );
}

/** calculates necessary data for translating seeeds and classifiers */
void Seeedpool::calcTranslationMapping(
   Seeedpool* origpool,
   std::vector<int>& rowothertothis,
   std::vector<int>& rowthistoother,
   std::vector<int>& colothertothis,
   std::vector<int>& colthistoother,
   std::vector<int>& missingrowinthis
   )
{
   int nrowsother = origpool->nConss;
   int nrowsthis = nConss;
   int ncolsother = origpool->nVars;
   int ncolsthis = nVars;

   std::vector<SCIP_CONS*> origscipconss = origpool->consToScipCons;
   std::vector<SCIP_CONS*> thisscipconss = consToScipCons;
   std::vector<SCIP_VAR*> origscipvars = origpool->varToScipVar;
   std::vector<SCIP_VAR*> thisscipvars = varToScipVar;

   assert(nrowsother == (int) origscipconss.size() );
   assert(nrowsthis == (int) thisscipconss.size() );

   assert(ncolsother == (int) origscipvars.size() );
   assert(ncolsthis == (int) thisscipvars.size() );


//   std::vector<SCIP_CONS*>::const_iterator origiter = origscipconss.begin();
//   std::vector<SCIP_CONS*>::const_iterator origiterend = origscipconss.end();
//
//   std::vector<SCIP_CONS*>::const_iterator thisiter = thisscipconss.begin();
//   std::vector<SCIP_CONS*>::const_iterator thisiterend = thisscipconss.end();
//
//   std::vector<SCIP_VAR*>::const_iterator origitervars = origscipvars.begin();
//   std::vector<SCIP_VAR*>::const_iterator origiterendvars = origscipvars.end();
//
//   std::vector<SCIP_VAR*>::const_iterator thisitervars = thisscipvars.begin();
//   std::vector<SCIP_VAR*>::const_iterator thisiterendvars = thisscipvars.end();


   rowothertothis.assign( nrowsother, - 1 );
   rowthistoother.assign( nrowsthis, - 1 );
   colothertothis.assign( ncolsother, - 1 );
   colthistoother.assign( ncolsthis, - 1 );

   missingrowinthis.clear();

   /* identify new and deleted rows and cols; and identify bijection between maintained variables */
   for( int i = 0; i < nrowsother ; ++i )
   {
      SCIP_CONS* otherrow = origscipconss[i];
      assert( otherrow != NULL );
      SCIP_Bool foundmaintained = false;
 //     thisiter = thisscipconss.begin();
      for( int j2 = i; j2 < nrowsthis + i; ++j2 )
      {
         int j = j2 % nrowsthis;
         SCIP_CONS* thisrow = thisscipconss[j];
         assert( SCIPconsIsTransformed( thisrow ) );
         char buffer[SCIP_MAXSTRLEN];
         assert( this->scip != NULL );
         strcpy( buffer, SCIPconsGetName( thisrow ) + 2 );
         assert( this->scip != NULL );
         if( strcmp( SCIPconsGetName( otherrow ), SCIPconsGetName( thisrow ) ) == 0 )
         {
            rowothertothis[i] = j;
            rowthistoother[j] = i;
            foundmaintained = true;
            break;
         }
      }
      if( ! foundmaintained )
         missingrowinthis.push_back( i );
   }

   for( int i = 0; i < ncolsother; ++i )
   {
      SCIP_VAR* othervar = origscipvars[i];
      for( int j2 = i; j2 < ncolsthis + i; ++j2 )
      {
         int j = j2 % ncolsthis;
         if( othervar == thisscipvars[j] )
         {
            colothertothis[i] = j;
            colthistoother[j] = i;
            break;
         }
      }
   }

   if ( FALSE )
   {
      for ( int i  = 0; i < (int) rowothertothis.size(); ++i )
         std::cout << (rowothertothis[i] == i) << " " ;

      std::cout << std::endl;

      for ( int i  = 0; i < (int) colothertothis.size(); ++i )
         std::cout << ( colothertothis[i] == i ) << " " ;
      std::cout << std::endl;
   }
}

/** returns translated seeeds derived from given mapping data */
std::vector<Seeed*> Seeedpool::getTranslatedSeeeds(
   std::vector<Seeed*>& origseeeds,
   std::vector<int>& rowothertothis,
   std::vector<int>& rowthistoother,
   std::vector<int>& colothertothis,
   std::vector<int>& colthistoother
   )
{
   std::vector<Seeed*> newseeeds( 0 );

   for( size_t s = 0; s < origseeeds.size(); ++ s )
   {
      SeeedPtr otherseeed;
      SeeedPtr newseeed;

      otherseeed = origseeeds[s];

      /** ignore seeeds with one block or no block, they are supposed to be found anyway */
      if( otherseeed->getNBlocks() == 1 || otherseeed->getNBlocks() == 0 )
         continue;

      SCIPverbMessage( this->scip, SCIP_VERBLEVEL_FULL, NULL, " transform seeed %d \n", otherseeed->getID() );

      newseeed = new Seeed( scip, this->getNewIdForSeeed(), this->getNConss(), this->getNVars() );

      /** prepare new seeed */
      newseeed->setNBlocks( otherseeed->getNBlocks() );

      newseeed->setUsergiven( otherseeed->getUsergiven() );

      /** set all (which have representative in the unpresolved seeed) constraints according to their representatives in the unpresolved seeed */
      for( int b = 0; b < otherseeed->getNBlocks(); ++ b )
      {
         for( int i = 0; i < otherseeed->getNConssForBlock( b ); i ++ )
         {
            int thiscons = rowothertothis[otherseeed->getConssForBlock( b )[i]];
            if( thiscons != - 1 )
            {
               newseeed->bookAsBlockCons( thiscons, b );
            }
         }
      }

      for( int i = 0; i < otherseeed->getNMasterconss(); i ++ )
      {
         int thiscons = rowothertothis[otherseeed->getMasterconss()[i]];
         if( thiscons != - 1 )
         {
            newseeed->bookAsMasterCons( thiscons );
         }
      }

      /** set linking and master vars according to their representatives in the unpresolved seeed */

      for( int j = 0; j < otherseeed->getNLinkingvars(); j ++ )
      {
         int thisvar = colothertothis[otherseeed->getLinkingvars()[j]];
         if( thisvar != - 1 )
         {
            newseeed->bookAsLinkingVar(thisvar);
         }
      }

      for( int j = 0; j < otherseeed->getNMastervars(); j ++ )
      {
         int thisvar = colothertothis[otherseeed->getMastervars()[j]];
         if( thisvar != - 1 )
         {
            newseeed->bookAsMasterVar( thisvar );
         }
      }

      newseeed->flushBooked();

      newseeed->setDetectorchain( otherseeed->getDetectorchainVector() );
      newseeed->setAncestorList( otherseeed->getAncestorList() );

      newseeed->addAncestorID( otherseeed->getID() );

      newseeed->copyClassifierStatistics( otherseeed );

      for( int i = 0; i < otherseeed->getNDetectors(); ++i )
      {
         newseeed->addClockTime( otherseeed->getDetectorClockTime( i ) );
         newseeed->addPctConssFromFree( otherseeed->getPctConssFromFree( i ) );
         newseeed->addPctConssToBlock( otherseeed->getPctConssToBlock( i ) );
         newseeed->addPctConssToBorder( otherseeed->getPctConssToBorder( i ) );
         newseeed->addPctVarsFromFree( otherseeed->getPctVarsFromFree( i ) );
         newseeed->addPctVarsToBlock( otherseeed->getPctVarsToBlock( i ) );
         newseeed->addPctVarsToBorder( otherseeed->getPctVarsToBorder( i ) );
         newseeed->addNNewBlocks( otherseeed->getNNewBlocks( i ) );
      }

      newseeed->setDetectorChainString( otherseeed->getDetectorChainString() );
      newseeed->setStemsFromUnpresolved( true );
      newseeed->setFinishedByFinisherUnpresolved( otherseeed->getFinishedByFinisher() );

      if( otherseeed->getFinishedByFinisher() )
         newseeed->setFinishedUnpresolvedBy( otherseeed->getDetectorchain()[otherseeed->getNDetectors() - 1] );

      newseeed->setFinishedByFinisher( otherseeed->getFinishedByFinisher() );
      newseeed->sort();
      newseeed->considerImplicits( this );
      newseeed->deleteEmptyBlocks(false);
      newseeed->setSeeedpool(this);
      newseeed->getScore( SCIPconshdlrDecompGetCurrScoretype( scip ) ) ;

      if( newseeed->checkConsistency( this ) )
         newseeeds.push_back( newseeed );
      else
      {
         delete newseeed;
         newseeed = NULL;
      }
   }

   return newseeeds;
}

/** returns translated ConsClassifiers derived from given mapping data */
std::vector<ConsClassifier*> Seeedpool::getTranslatedConsClassifiers(
   std::vector<ConsClassifier*>& otherclassifiers,
   std::vector<int>& rowothertothis,
   std::vector<int>& rowthistoother
   )
{
   std::vector<ConsClassifier*> newclassifiers( 0 );

   for( size_t i = 0; i < otherclassifiers.size(); ++ i )
   {
      ConsClassifier* oldclassifier = otherclassifiers[i];
      ConsClassifier* newclassifier;
      std::stringstream newname;

      newname << oldclassifier->getName() << "-origp";
      newclassifier = new ConsClassifier( scip, newname.str().c_str(), oldclassifier->getNClasses(),
         (int) rowthistoother.size() );
      int bufferclassindex = - 1;

      /** copy class information */
      for( int j = 0; j < oldclassifier->getNClasses(); ++ j )
      {
         newclassifier->setClassName( j, oldclassifier->getClassName( j ) );
         newclassifier->setClassDescription( j, oldclassifier->getClassDescription( j ) );
         newclassifier->setClassDecompInfo( j, oldclassifier->getClassDecompInfo( j ) );
      }

      /** assign new conss to classes */
      for( int c = 0; c < (int) rowthistoother.size(); ++ c )
      {
         if( rowthistoother[c] != - 1 )
         {
            newclassifier->assignConsToClass( c, oldclassifier->getClassOfCons( rowthistoother[c] ) );
         }
         else
         {
            if( bufferclassindex == - 1 )
            {
               bufferclassindex = newclassifier->addClass( "buffer",
                  "This class contains constraints which are new in the presolved problem.", BOTH );
            }
            newclassifier->assignConsToClass( c, bufferclassindex );
         }
      }

      /** remove empty classes */
      newclassifier->removeEmptyClasses();

      newclassifiers.push_back( newclassifier );
   }

   return newclassifiers;
}

/** returns translated VarClassifiers derived from given mapping data */
std::vector<VarClassifier*> Seeedpool::getTranslatedVarClassifiers(
   std::vector<VarClassifier*>& otherclassifiers,
   std::vector<int>& colothertothis,
   std::vector<int>& colthistoother
   )
{
   std::vector<VarClassifier*> newclassifiers( 0 );

   for( size_t i = 0; i < otherclassifiers.size(); ++ i )
   {
      VarClassifier* oldclassifier = otherclassifiers[i];
      VarClassifier* newclassifier;
      std::stringstream newname;

      newname << oldclassifier->getName() << "-origp";
      newclassifier = new VarClassifier( scip, newname.str().c_str(), oldclassifier->getNClasses(),
         (int) colthistoother.size() );
      int bufferclassindex = - 1;

      /** copy class information */
      for( int j = 0; j < oldclassifier->getNClasses(); ++ j )
      {
         newclassifier->setClassName( j, oldclassifier->getClassName( j ) );
         newclassifier->setClassDescription( j, oldclassifier->getClassDescription( j ) );
         newclassifier->setClassDecompInfo( j, oldclassifier->getClassDecompInfo( j ) );
      }

      /** assign new vars to classes */
      for( int c = 0; c < (int) colthistoother.size(); ++ c )
      {
         if( colthistoother[c] != - 1 )
         {
            newclassifier->assignVarToClass( c, oldclassifier->getClassOfVar( colthistoother[c] ) );
         }
         else
         {
            if( bufferclassindex == - 1 )
            {
               bufferclassindex = newclassifier->addClass( "buffer",
                  "This class contains variables which are new in the presolved problem.", ALL );
            }
            newclassifier->assignVarToClass( c, bufferclassindex );
         }
      }

      /** remove empty classes */
      newclassifier->removeEmptyClasses();

      newclassifiers.push_back( newclassifier );
   }

   return newclassifiers;
}

/** registers translated seeeds from the original problem */
void Seeedpool::populate(
   std::vector<SeeedPtr> seeeds
   )
{
   seeedstopopulate = seeeds;
}

/** sorts the seeed and calculates a its implicit assignments, hashvalue and evaluation */
SCIP_RETCODE Seeedpool::prepareSeeed(
   SeeedPtr seeed
   )
{
   seeed->considerImplicits( this );
   seeed->calcHashvalue();
   seeed->setSeeedpool(this);
   //seeed->evaluate( this, SCIPconshdlrDecompGetCurrScoretype( scip ) );

   return SCIP_OKAY;
}

bool Seeedpool::isConsCardinalityCons(
      int  consindexd
      )
{
   SCIP_CONS* cons;

   cons = consToScipCons[consindexd];

   assert(cons != NULL);

   return GCGgetConsIsCardinalityCons(scip, cons);


}

/** is cons with specified indec partitioning, packing, or covering constraint?*/
bool Seeedpool::isConsSetppc(
      int  consindexd
      )
{
   SCIP_CONS* cons;

   cons = consToScipCons[consindexd];

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "setppc") == 0 )
         {
            switch( SCIPgetTypeSetppc(scip, cons) )
            {
            case SCIP_SETPPCTYPE_COVERING:
               return true;
               break;
            case SCIP_SETPPCTYPE_PARTITIONING:
               return true;

            case SCIP_SETPPCTYPE_PACKING:
               return true;
            }
         }
         else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "logicor") == 0 )
         {
            return true;
         }
         else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0 )
         {
            SCIP_SETPPCTYPE type;

            if( GCGgetConsIsSetppc(scip, cons, &type) )
            {
               switch( type )
               {
               case SCIP_SETPPCTYPE_COVERING:
                  return true;
               case SCIP_SETPPCTYPE_PARTITIONING:
               return true;
               case SCIP_SETPPCTYPE_PACKING:
               return true;
               }
            }

         }

   return false;
}

/** is cons with specified indec partitioning, or packing  constraint?*/
bool Seeedpool::isConsSetpp(
      int  consindexd
      )
{
   SCIP_CONS* cons;

   cons = consToScipCons[consindexd];

   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "setppc") == 0 )
         {
            switch( SCIPgetTypeSetppc(scip, cons) )
            {
            case SCIP_SETPPCTYPE_PARTITIONING:
               return true;

            case SCIP_SETPPCTYPE_PACKING:
               return true;
            case SCIP_SETPPCTYPE_COVERING:
               return false;

            }
         }
         else if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "linear") == 0 )
         {
            SCIP_SETPPCTYPE type;

            if( GCGgetConsIsSetppc(scip, cons, &type) )
            {
               switch( type )
               {
               case SCIP_SETPPCTYPE_PARTITIONING:
               return true;
               case SCIP_SETPPCTYPE_PACKING:
               return true;
               case SCIP_SETPPCTYPE_COVERING:
                  return false;

               }
            }

         }

   return false;
}




/** sorts seeeds in allrelevantseeeds data structure by ascending id */
void Seeedpool::sortAllRelevantSeeeds()
{
   int maxid = 0;
   std::vector<SeeedPtr> tmpAllRelevantSeeeds( 0 );

   for( size_t i = 0; i < ancestorseeeds.size(); ++ i )
   {
      if( ancestorseeeds[i]->getID() > maxid )
         maxid = ancestorseeeds[i]->getID();
   }

   tmpAllRelevantSeeeds = std::vector < SeeedPtr > ( maxid + 1, NULL );

   for( size_t i = 0; i < ancestorseeeds.size(); ++ i )
   {
      if( ancestorseeeds[i]->getID() < 0 )
         continue;
      tmpAllRelevantSeeeds[ancestorseeeds[i]->getID()] = ancestorseeeds[i];
   }

   ancestorseeeds = tmpAllRelevantSeeeds;
}

/** returns the variable indices of the matrix for a constraint */
const int* Seeedpool::getVarsForCons(
   int cons
   )
{
   return & varsForConss[cons][0];
}

/** returns the coefficients of the matrix for a constraint */
const SCIP_Real * Seeedpool::getValsForCons(
   int cons
   )
{
   return & valsForConss[cons][0];
}

/** returns the constraint indices of the coefficient matrix for a variable */
const int* Seeedpool::getConssForVar(
   int var
   )
{
   return & conssForVars[var][0];
}

/** returns the constraint indices of the coefficient matrix for a constraint */
const int* Seeedpool::getConssForCons(
   int cons
   )
{
   return & conssadjacencies[cons][0];
}



/** returns the number of variables for a given constraint */
int Seeedpool::getNVarsForCons(
   int cons
   )
{
   return varsForConss[cons].size();
}

/** returns the number of constraints for a given variable */
int Seeedpool::getNConssForVar(
   int var
   )
{
   return conssForVars[var].size();
}

/** returns the number of constraints for a given variable */
int Seeedpool::getNConssForCons(
   int cons
   )
{
   return conssadjacencies[cons].size();
}


/** returns the SCIP variable related to a variable index */
SCIP_VAR* Seeedpool::getVarForIndex(
   int varIndex
   )
{
   return varToScipVar[varIndex];
}

/** returns the SCIP constraint related to a constraint index */
SCIP_CONS* Seeedpool::getConsForIndex(
   int consIndex
   )
{
   return consToScipCons[consIndex];
}

/** returns the SCIP detector related to a detector index */
DEC_DETECTOR* Seeedpool::getDetectorForIndex(
   int detectorIndex
   )
{
   return detectorToScipDetector[detectorIndex];
}

/** returns the SCIP detector related to a finishing detector index */
DEC_DETECTOR* Seeedpool::getFinishingDetectorForIndex(
   int detectorIndex
   )
{
   return detectorToFinishingScipDetector[detectorIndex];
}


/** returns the SCIP detector related to a postprocessing detector index */
DEC_DETECTOR* Seeedpool::getPostprocessingDetectorForIndex(
   int detectorIndex
   )
{
   return detectorToPostprocessingScipDetector[detectorIndex];
}


/** returns a coefficient from the coefficient matrix */
SCIP_Real Seeedpool::getVal(
   int row,
   int col
   )
{
   std::tr1::unordered_map<std::pair<int, int>, SCIP_Real, pair_hash>::const_iterator iter = valsMap.find(
      std::pair<int, int>( row, col ) );

   if( iter == valsMap.end() )
      return 0.;

   return iter->second;
}

/** returns the variable index related to a SCIP variable */
int Seeedpool::getIndexForVar(
   SCIP_VAR* var
   )
{
   return scipVarToIndex[var];
}

/** returns the constraint index related to a SCIP constraint */
int Seeedpool::getIndexForCons(
   SCIP_CONS* cons
   )
{
   return scipConsToIndex[cons];
}

/** returns the detector index related to a detector */
int Seeedpool::getIndexForDetector(
   DEC_DETECTOR* detector
   )
{
   return scipDetectorToIndex[detector];
}

/** returns the finishing detector index related to a detector */
int Seeedpool::getIndexForFinishingDetector(
   DEC_DETECTOR* detector
   )
{
   return scipFinishingDetectorToIndex[detector];
}

/** returns the postprocessing detector index related to a detector */
int Seeedpool::getIndexForPostprocessingDetector(
   DEC_DETECTOR* detector
   )
{
   return scipPostprocessingDetectorToIndex[detector];
}



/** returns a new unique id for a seeed */
int Seeedpool::getNewIdForSeeed()
{
   return SCIPconshdlrDecompGetNextSeeedID( scip );
}

/** returns the number of detectors used in the seeedpool */
int Seeedpool::getNDetectors()
{
   return nDetectors;
}

/** returns the number of nonzero entries in the coefficient matrix */
int Seeedpool::getNNonzeros()
{
   return nnonzeros;
}

/** returns the number of finishing detectors used in the seeedpool */
int Seeedpool::getNFinishingDetectors()
{
   return nFinishingDetectors;
}

/** returns the number of postprocessing detectors used in the seeedpool */
int Seeedpool::getNPostprocessingDetectors()
{
   return nPostprocessingDetectors;
}


/** returns the number of variables considered in the seeedpool */
int Seeedpool::getNVars()
{
   return nVars;
}

/** returns the number of constraints considered in the seeedpool */
int Seeedpool::getNConss()
{
   return nConss;
}

/* returns associated scip */

SCIP* Seeedpool::getScip()
{
   return scip;
}

/** returns the candidates for block size sorted in descending order by how often a candidate was added */
std::vector<int> Seeedpool::getSortedCandidatesNBlocks()
{
   std::vector<int> toreturn( 0 );
   SCIP_Bool output = false;

   /** first: get the block number candidates directly given by the user */
   if( output && !usercandidatesnblocks.empty() )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "number of user block number candidates: %d  \n", usercandidatesnblocks.size() );
   }

   for( size_t i = 0; i < usercandidatesnblocks.size() ; ++i)
   {
      toreturn.push_back( usercandidatesnblocks[i]);

      if( output )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%d  \n", usercandidatesnblocks[i] );
      }

   }

   /** second: sort the current candidates */
   std::sort( candidatesNBlocks.begin(), candidatesNBlocks.end(), sort_decr() );

   /** optional: print sorted candidates */
   if( output )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "nCandidates: %d ", candidatesNBlocks.size() );


      for( size_t i = 0; i < candidatesNBlocks.size(); ++ i )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "nblockcandides: %d , %d times prop", candidatesNBlocks[i].first, candidatesNBlocks[i].second );
   }

   /** secondly: push candidates to output vector */
   for( size_t i = 0; i < candidatesNBlocks.size(); ++ i )
      toreturn.push_back( candidatesNBlocks[i].first );

   return toreturn;
}

/** adds a candidate for block size and counts how often a candidate is added */
void Seeedpool::addCandidatesNBlocks(
   int candidate
   )
{
   if( candidate > 1 )
   {
      bool alreadyIn = false;
      for( size_t i = 0; i < candidatesNBlocks.size(); ++ i )
      {
         if( candidatesNBlocks[i].first == candidate )
         {
            alreadyIn = true;
            ++ candidatesNBlocks[i].second;
            break;
         }
      }
      if( ! alreadyIn )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "added block number candidate : %d \n ", candidate );
         candidatesNBlocks.push_back( std::pair<int, int>( candidate, 1 ) );
      }
   }
}

/** adds a candidate for block size given by the user */
void Seeedpool::addUserCandidatesNBlocks(
   int candidate
   )
{
   bool alreadyIn = false;
   for( size_t i = 0; i < candidatesNBlocks.size(); ++i )
   {
      if( usercandidatesnblocks[i] == candidate )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "is already given by the user as a block number candidate, there is no advantage in adding it twice \n " );
         return;
      }
   }
   if( !alreadyIn )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "added user block number candidate : %d \n ", candidate );
      usercandidatesnblocks.push_back(candidate);
   }
}

/** returns number of user-given block size candidates */
int Seeedpool::getNUserCandidatesNBlocks()
{
   return (int) usercandidatesnblocks.size();
}

/** calculates and adds block size candidates using constraint classifications and variable classifications */
void Seeedpool::calcCandidatesNBlocks()
{
   /* strategy: for every subset of constraint classes and variable classes calculate gcd (greatest common divisors)
    * of the corresponding number of constraints/variables assigned to this class */

   /* if  distribution of classes exceeds this number it is skipped */
   int maximumnclasses = 18;

   /** firstly, iterate over all consclassifiers */
   for( size_t classifier = 0; classifier < consclassescollection.size(); ++ classifier )
   {
      /** check if there are too many classes in this distribution and skip it if so */
      if( consclassescollection[classifier]->getNClasses() > maximumnclasses )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " the current consclass distribution includes %d classes but only %d are allowed for calcCandidatesNBlocks()\n ", consclassescollection[classifier]->getNClasses(), maximumnclasses );
         continue;
      }

      /** get necessary data of current classifier */
      std::vector < std::vector<int> > subsetsOfConstypes = consclassescollection[classifier]->getAllSubsets( true, true,
         true );
      std::vector<int> nConssOfClasses = consclassescollection[classifier]->getNConssOfClasses();

      /** start with the cardinalities of the consclasses as candidates */
      for( size_t i = 0; i < nConssOfClasses.size(); ++ i )
      {
         addCandidatesNBlocks( nConssOfClasses[i] );
      }

      /** continue with gcd of all cardinalities in this subset */
      for( size_t subset = 0; subset < subsetsOfConstypes.size(); ++ subset )
      {
         int greatestCD = 1;

         if( subsetsOfConstypes[subset].size() == 0 || subsetsOfConstypes[subset].size() == 1 )
            continue;

         greatestCD = gcd( nConssOfClasses[subsetsOfConstypes[subset][0]], nConssOfClasses[subsetsOfConstypes[subset][1]] );

         for( size_t i = 2; i < subsetsOfConstypes[subset].size(); ++ i )
         {
            greatestCD = gcd( greatestCD, nConssOfClasses[subsetsOfConstypes[subset][i]] );
         }

         addCandidatesNBlocks( greatestCD );
      }
   }

   /** secondly, iterate over all varclassifiers */
   for( size_t classifier = 0; classifier < varclassescollection.size(); ++ classifier )
   {
      /** check if there are too many classes in this distribution and skip it if so */
      if( varclassescollection[classifier]->getNClasses() > maximumnclasses )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " the current varclass distribution includes %d classes but only %d are allowed for calcCandidatesNBlocks()\n ", varclassescollection[classifier]->getNClasses(), maximumnclasses );
         continue;
      }

      /** get necessary data of current classifier */
      std::vector < std::vector<int> > subsetsOfVartypes = varclassescollection[classifier]->getAllSubsets( true, true, true,
         true );
      std::vector<int> nVarsOfClasses = varclassescollection[classifier]->getNVarsOfClasses();

      /** start with the cardinalities of the varclasses as candidates */
      for( size_t i = 0; i < nVarsOfClasses.size(); ++ i )
      {
         addCandidatesNBlocks( nVarsOfClasses[i] );
      }

      /** continue with gcd of all cardinalities in this subset */
      for( size_t subset = 0; subset < subsetsOfVartypes.size(); ++ subset )
      {
         int greatestCD = 1;

         if( subsetsOfVartypes[subset].size() == 0 || subsetsOfVartypes[subset].size() == 1 )
            continue;

         greatestCD = gcd( nVarsOfClasses[subsetsOfVartypes[subset][0]], nVarsOfClasses[subsetsOfVartypes[subset][1]] );

         for( size_t i = 2; i < subsetsOfVartypes[subset].size(); ++ i )
         {
            greatestCD = gcd( greatestCD, nVarsOfClasses[subsetsOfVartypes[subset][i]] );
         }

         addCandidatesNBlocks( greatestCD );
      }
   }
}

/** adds a constraint classifier if it is no duplicate of an existing classifier */
void Seeedpool::addConsClassifier(
   ConsClassifier* givenClassifier
   )
{
   SCIP_Bool detectionstatistics;

   SCIPgetBoolParam(scip, "detection/allowclassifierduplicates/enabled", &detectionstatistics);

   if( givenClassifier != NULL )
   {
      /** check whether there already exists an equivalent consclassifier */
      ConsClassifier* equiv = NULL;

      for( size_t i = 0; !detectionstatistics && i < consclassescollection.size(); ++ i )
      {
         if( givenClassifier->classifierIsDuplicateOfClassifier( consclassescollection[i] ) )
         {
            equiv = consclassescollection[i];
            break;
         }
      }

      if( equiv == NULL )
         consclassescollection.push_back( givenClassifier );
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " consclassifier %s is not considered since it offers the same structure as  %s  consclassifier\n ", givenClassifier->getName(), equiv->getName() );
         delete givenClassifier;
      }
   }
}

/** returns a new constraint classifier
 *  where all constraints with identical SCIP constype are assigned to the same class */
ConsClassifier* Seeedpool::createConsClassifierForSCIPConstypes()
{
   std::vector<consType> foundConstypes( 0 );
   std::vector<int> constypesIndices( 0 );
   std::vector<int> classForCons = std::vector<int>( getNConss(), - 1 );
   ConsClassifier* classifier;

   /** firstly, assign all constraints to classindices */
   for( int i = 0; i < getNConss(); ++ i )
   {
      SCIP_CONS* cons;
      bool found = false;
      cons = getConsForIndex( i );
      consType cT = GCGconsGetType( cons );
      size_t constype;

      /** check whether the constraint's constype is new */
      for( constype = 0; constype < foundConstypes.size(); ++ constype )
      {
         if( foundConstypes[constype] == cT )
         {
            found = true;
            break;
         }
      }
      /** if it is new, create a new classindex */
      if( ! found )
      {
         foundConstypes.push_back( GCGconsGetType( cons ) );
         classForCons[i] = foundConstypes.size() - 1;
      }
      else
         classForCons[i] = constype;
   }

   /** secondly, use these information to create a ConsClassifier */
   classifier = new ConsClassifier( scip, "constypes", (int) foundConstypes.size(), getNConss() );

   /** set class names and descriptions of every class */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::string name;
      std::stringstream text;
      switch( foundConstypes[c] )
      {
         case linear:
            name = "linear";
            break;
         case knapsack:
            name = "knapsack";
            break;
         case varbound:
            name = "varbound";
            break;
         case setpacking:
            name = "setpacking";
            break;
         case setcovering:
            name = "setcovering";
            break;
         case setpartitioning:
            name = "setpartitioning";
            break;
         case logicor:
            name = "logicor";
            break;
         case sos1:
            name = "sos1";
            break;
         case sos2:
            name = "sos2";
            break;
         case unknown:
            name = "unknown";
            break;
         case nconsTypeItems:
            name = "nconsTypeItems";
            break;
         default:
            name = "newConstype";
            break;
      }
      classifier->setClassName( c, name.c_str() );
      text << "This class contains all constraints that are of (SCIP) constype \"" << name << "\".";
      classifier->setClassDescription( c, text.str().c_str() );
   }

   /** copy the constraint assignment information found in first step */
   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass( i, classForCons[i] );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier %s yields a classification with %d  different constraint classes \n", classifier->getName(), (int) foundConstypes.size() );
   return classifier;
}


/** returns a new constraint classifier
 *  where all constraints with identical SCIP constype are assigned to the same class */
ConsClassifier* Seeedpool::createConsClassifierForMiplibConstypes()
{
   std::vector<int> nfoundconstypesrangedsinglecount( (int) SCIP_CONSTYPE_GENERAL + 1, 0 );
   std::vector<int> nfoundconstypesrangeddoublecount( (int) SCIP_CONSTYPE_GENERAL + 1, 0 );

//   std::vector<int> constypesIndices( 0 );
   std::vector<int> classforcons = std::vector<int>( getNConss(), -1 );
   ConsClassifier* classifier;

   /** firstly, assign all constraints to classindices */
   for( int c = 0; c < getNConss(); ++ c )
   {
      SCIP_CONS* cons;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real* vals;
      SCIP_VAR** vars;
      int nvars;
      int i;

      cons = getConsForIndex( c );

      nvars =  GCGconsGetNVars(scip, cons );

      lhs = GCGconsGetLhs(scip, cons);
      rhs = GCGconsGetRhs(scip, cons);
      if( nvars != 0 )
      {
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vals, nvars));
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vars, nvars));
         SCIP_CALL_ABORT( GCGconsGetVals(scip, cons, vals, nvars ) );
         SCIP_CALL_ABORT( GCGconsGetVars(scip, cons, vars, nvars ) );
      }

      for( i = 0; i < nvars; i++ )
      {
         assert(!SCIPisZero(scip, vals[i]) );
      }


      /* is constraint of type SCIP_CONSTYPE_EMPTY? */
      if( nvars == 0 )
      {
         SCIPdebugMsg(scip, "classified as EMPTY: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_EMPTY]++;
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_EMPTY]++;
         classforcons[c] = SCIP_CONSTYPE_EMPTY;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_FREE? */
      if( SCIPisInfinity(scip, rhs) && SCIPisInfinity(scip, -lhs) )
      {
         SCIPdebugMsg(scip, "classified as FREE: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_FREE]++;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_FREE]++;
         classforcons[c] = SCIP_CONSTYPE_FREE;
         SCIPfreeBufferArray(scip, &vars);
         SCIPfreeBufferArray(scip, &vals);
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_SINGLETON? */
      if( nvars == 1 )
      {
         SCIPdebugMsg(scip, "classified as SINGLETON: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SINGLETON] += 2 ;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SINGLETON]++;
         classforcons[c] = SCIP_CONSTYPE_SINGLETON;
         SCIPfreeBufferArray(scip, &vars) ;
         SCIPfreeBufferArray(scip, &vals) ;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_AGGREGATION? */
      if( nvars == 2 && SCIPisEQ(scip, lhs, rhs) )
      {
         SCIPdebugMsg(scip, "classified as AGGREGATION: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_AGGREGATION]++;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_AGGREGATION]++;
         classforcons[c] = SCIP_CONSTYPE_AGGREGATION;
         SCIPfreeBufferArray(scip, &vars) ;
         SCIPfreeBufferArray(scip, &vals) ;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{VARBOUND}? */
      if( nvars == 2 )
      {
         SCIPdebugMsg(scip, "classified as VARBOUND: ");
         SCIPdebugPrintCons(scip, cons, NULL);
         nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_VARBOUND] += 2 ;
         nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_VARBOUND]++;
         classforcons[c] = SCIP_CONSTYPE_VARBOUND;
         SCIPfreeBufferArray(scip, &vars) ;
         SCIPfreeBufferArray(scip, &vals) ;
         continue;
      }

      /* is constraint of type SCIP_CONSTYPE_{SETPARTITION, SETPACKING, SETCOVERING, CARDINALITY, INVKNAPSACK}? */
      {
         SCIP_Real scale;
         SCIP_Real b;
         SCIP_Bool unmatched;
         int nnegbinvars;

         unmatched = FALSE;
         nnegbinvars = 0;

         scale = REALABS(vals[0]);
         for( i = 0; i < nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(scip, SCIPvarGetLbGlobal(vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(scip, SCIPvarGetUbGlobal(vars[i]), 2.0);
            unmatched = unmatched || !SCIPisEQ(scip, REALABS(vals[i]), scale);

            if( vals[i] < 0.0 )
               nnegbinvars++;
         }

         if( !unmatched )
         {
            if( SCIPisEQ(scip, lhs, rhs) )
            {
               b = rhs/scale + nnegbinvars;
               if( SCIPisEQ(scip, 1.0, b) )
               {
                  SCIPdebugMsg(scip, "classified as SETPARTITION: ");
                  SCIPdebugPrintCons(scip, cons, NULL);
                  nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SETPARTITION] += 1 ;
                  nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SETPARTITION]++;
                  classforcons[c] = SCIP_CONSTYPE_SETPARTITION;
                  SCIPfreeBufferArray(scip, &vars) ;
                  SCIPfreeBufferArray(scip, &vals) ;
                  continue;
               }
               else if( SCIPisIntegral(scip, b) && !SCIPisNegative(scip, b) )
               {
                  SCIPdebugMsg(scip, "classified as CARDINALITY: ");
                  SCIPdebugPrintCons(scip, cons, NULL);
                  nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_CARDINALITY] += 1 ;
                  nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_CARDINALITY]++;
                  classforcons[c] = SCIP_CONSTYPE_CARDINALITY;
                  SCIPfreeBufferArray(scip, &vars);
                  SCIPfreeBufferArray(scip, &vals);
                  continue;
               }
            }

            b = rhs/scale + nnegbinvars;
            if( SCIPisEQ(scip, 1.0, b) )
            {
               SCIPdebugMsg(scip, "classified as SETPACKING: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SETPACKING] += 1 ;
               nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SETPACKING]++;
               classforcons[c] = SCIP_CONSTYPE_SETPACKING;
               rhs = SCIPinfinity(scip);
            }
            else if( SCIPisIntegral(scip, b) && !SCIPisNegative(scip, b) )
            {
               SCIPdebugMsg(scip, "classified as INVKNAPSACK: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_INVKNAPSACK] += 1 ;
                nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_INVKNAPSACK]++;
                classforcons[c] = SCIP_CONSTYPE_INVKNAPSACK;
               rhs = SCIPinfinity(scip);
            }

            b = lhs/scale + nnegbinvars;
            if( SCIPisEQ(scip, 1.0, b) )
            {
               SCIPdebugMsg(scip, "classified as SETCOVERING: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_SETCOVERING] += 1 ;
               nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_SETCOVERING]++;
               classforcons[c] = SCIP_CONSTYPE_SETCOVERING;
               lhs = -SCIPinfinity(scip);
            }

            if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
            {
               SCIPfreeBufferArray(scip, &vars);
               SCIPfreeBufferArray(scip, &vals);
               continue;
            }
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{EQKNAPSACK, BINPACKING, KNAPSACK}? */
      /* @todo If coefficients or rhs are not integral, we currently do not check
       * if the constraint could be scaled (finitely), such that they are.
       */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         b = rhs;
         unmatched = FALSE;
         for( i = 0; i < nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisLE(scip, SCIPvarGetLbGlobal(vars[i]), -1.0);
            unmatched = unmatched || SCIPisGE(scip, SCIPvarGetUbGlobal(vars[i]), 2.0);
            unmatched = unmatched || !SCIPisIntegral(scip, vals[i]);

            if( SCIPisNegative(scip, vals[i]) )
               b -= vals[i];
         }
         unmatched = unmatched || !isFiniteNonnegativeIntegral(scip, b);

         if( !unmatched )
         {
            if( SCIPisEQ(scip, lhs, rhs) )
            {
               SCIPdebugMsg(scip, "classified as EQKNAPSACK: ");
               SCIPdebugPrintCons(scip, cons, NULL);
               nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_EQKNAPSACK] += 1 ;
               nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_EQKNAPSACK]++;
               classforcons[c] = SCIP_CONSTYPE_EQKNAPSACK;
               SCIPfreeBufferArray(scip, &vars);
               SCIPfreeBufferArray(scip, &vals);
               continue;
            }
            else
            {
               SCIP_Bool matched;

               matched = FALSE;
               for( i = 0; i < nvars && !matched; i++ )
               {
                  matched = matched || SCIPisEQ(scip, b, REALABS(vals[i]));
               }

               SCIPdebugMsg(scip, "classified as %s: ", matched ? "BINPACKING" : "KNAPSACK");
               SCIPdebugPrintCons(scip, cons, NULL);
               nfoundconstypesrangeddoublecount[matched ? SCIP_CONSTYPE_BINPACKING : SCIP_CONSTYPE_KNAPSACK] += 1 ;
               nfoundconstypesrangedsinglecount[matched ? SCIP_CONSTYPE_BINPACKING : SCIP_CONSTYPE_KNAPSACK]++;
               classforcons[c] = matched ? SCIP_CONSTYPE_BINPACKING : SCIP_CONSTYPE_KNAPSACK;

            }

            if( SCIPisInfinity(scip, -lhs) )
            {
               SCIPfreeBufferArray(scip, &vars);
               SCIPfreeBufferArray(scip, &vals);
               continue;
            }
            else
               rhs = SCIPinfinity(scip);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{INTKNAPSACK}? */
      {
         SCIP_Real b;
         SCIP_Bool unmatched;

         unmatched = FALSE;

         b = rhs;
         unmatched = unmatched || !isFiniteNonnegativeIntegral(scip, b);

         for( i = 0; i < nvars && !unmatched; i++ )
         {
            unmatched = unmatched || SCIPvarGetType(vars[i]) == SCIP_VARTYPE_CONTINUOUS;
            unmatched = unmatched || SCIPisNegative(scip, SCIPvarGetLbGlobal(vars[i]));
            unmatched = unmatched || !SCIPisIntegral(scip, vals[i]);
            unmatched = unmatched || SCIPisNegative(scip, vals[i]);
         }

         if( !unmatched )
         {
            SCIPdebugMsg(scip, "classified as INTKNAPSACK: ");
            SCIPdebugPrintCons(scip, cons, NULL);
            nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_INTKNAPSACK] += 1 ;
            nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_INTKNAPSACK]++;
            classforcons[c] = SCIP_CONSTYPE_INTKNAPSACK;

            if( SCIPisInfinity(scip, -lhs) )
            {
               SCIPfreeBufferArray(scip, &vars);
               SCIPfreeBufferArray(scip, &vals);
               continue;
            }
            else
               rhs = SCIPinfinity(scip);
         }
      }

      /* is constraint of type SCIP_CONSTYPE_{MIXEDBINARY}? */
      {
         SCIP_Bool unmatched;

         unmatched = FALSE;
         for( i = 0; i < nvars && !unmatched; i++ )
         {
            if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS
               && (SCIPisLE(scip, SCIPvarGetLbGlobal(vars[i]), -1.0)
                  || SCIPisGE(scip, SCIPvarGetUbGlobal(vars[i]), 2.0)) )
               unmatched = TRUE;
         }

         if( !unmatched )
         {
            SCIPdebugMsg(scip, "classified as MIXEDBINARY (%d): ", isRangedRow(scip, lhs, rhs) ? 2 : 1);
            SCIPdebugPrintCons(scip, cons, NULL);
            nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_MIXEDBINARY] += 1 ;
            nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_MIXEDBINARY]++;
            classforcons[c] = SCIP_CONSTYPE_MIXEDBINARY;
            SCIPfreeBufferArray(scip, &vars) ;
            SCIPfreeBufferArray(scip, &vals) ;
            continue;

         }
      }

      /* no special structure detected */
      SCIPdebugMsg(scip, "classified as GENERAL: ");
      SCIPdebugPrintCons(scip, cons, NULL);
      nfoundconstypesrangeddoublecount[SCIP_CONSTYPE_GENERAL] += 1 ;
      nfoundconstypesrangedsinglecount[SCIP_CONSTYPE_GENERAL]++;
      classforcons[c] = SCIP_CONSTYPE_GENERAL;
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals);
   }




   classifier = new ConsClassifier( scip, "constypes according to miplip", (int) SCIP_CONSTYPE_GENERAL + 1, getNConss() );

#ifdef WRITE_ORIG_CONSTYPES
   std::ofstream myfile;
   myfile.open ("origconstypes.csv", std::ios::app );
   myfile << SCIPgetProbName(scip) << ", ";
#endif



   /** set class names and descriptions of every class */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::string name;
      std::stringstream text;
      switch( c )
      {
         case (int) SCIP_CONSTYPE_EMPTY:
            name = "empty";
            break;
         case SCIP_CONSTYPE_FREE:
            name = "free";
            break;
         case SCIP_CONSTYPE_SINGLETON:
            name = "singleton";
            break;
         case SCIP_CONSTYPE_AGGREGATION:
            name = "aggregation";
            break;
         case SCIP_CONSTYPE_VARBOUND:
            name = "varbound";
            break;
         case SCIP_CONSTYPE_SETPARTITION:
            name = "setpartition";
            break;
         case SCIP_CONSTYPE_SETPACKING:
            name = "setpacking";
            break;
         case SCIP_CONSTYPE_SETCOVERING:
            name = "setcovering";
            break;
         case SCIP_CONSTYPE_CARDINALITY:
            name = "cardinality";
            break;
         case SCIP_CONSTYPE_INVKNAPSACK:
            name = "invknapsack";
            break;
         case SCIP_CONSTYPE_EQKNAPSACK:
            name = "eqknapsack";
            break;
         case SCIP_CONSTYPE_BINPACKING:
            name = "binpacking";
            break;
         case SCIP_CONSTYPE_KNAPSACK:
            name = "knapsack";
            break;
         case SCIP_CONSTYPE_INTKNAPSACK:
            name = "intknapsack";
            break;
         case SCIP_CONSTYPE_MIXEDBINARY:
            name = "mixed binary";
            break;
         case SCIP_CONSTYPE_GENERAL:
            name = "general";
            break;
         default:
            name = "unknown";
            break;


      }


#ifdef WRITE_ORIG_CONSTYPES
         myfile << " " <<  nfoundconstypesrangeddoublecount[c] << ",";
#endif

      classifier->setClassName( c, name.c_str() );
      text << "This class contains all constraints that are of (miplib) constype \"" << name << "\".";
      classifier->setClassDescription( c, text.str().c_str() );
   }

#ifdef WRITE_ORIG_CONSTYPES
      myfile << std::endl;
      myfile.close();
#endif



   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass( i, classforcons[i] );
   }



   classifier->removeEmptyClasses();
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier %s yields a classification with %d  different constraint classes \n", classifier->getName(), classifier->getNClasses() );

   return classifier;
}







/** returns a new constraint classifier
 *  where all constraints with identical consname (ignoring digits) are assigned to the same class */
ConsClassifier* Seeedpool::createConsClassifierForConsnamesDigitFreeIdentical()
{
   std::vector < std::string > consnamesToCompare( getNConss(), "" );
   std::vector<int> nConssConstype( 0 );
   std::vector<int> classForCons = std::vector<int>( getNConss(), - 1 );
   std::vector < std::string > nameClasses( 0 );
   ConsClassifier* classifier;

   /** firstly, remove all digits from the consnames */
   for( int i = 0; i < getNConss(); ++ i )
   {
      int nremoved;
      char consname[SCIP_MAXSTRLEN];
      strcpy( consname, SCIPconsGetName( getConsForIndex( i ) ) );

      removeDigits( consname, & nremoved );
      consnamesToCompare[i] = std::string( consname );
   }

   for( int i = 0; i < getNConss(); ++ i )
   {
      /** check if string belongs to an existing name class */
      bool belongstoexistingclass = false;

      for( size_t j = 0; j < nameClasses.size(); ++ j )
      {
         if( nameClasses[j].compare( consnamesToCompare[i] ) == 0 )
         {
            belongstoexistingclass = true;
            classForCons[i] = j;
            nConssConstype[j] ++;
            break;
         }
      }
      /** if not, create a new class */
      if( ! belongstoexistingclass )
      {
         nameClasses.push_back( consnamesToCompare[i] );
         nConssConstype.push_back( 1 );
         classForCons[i] = nameClasses.size() - 1;

      }
   }

   /** secondly, use these information to create a ConsClassifier */
   classifier = new ConsClassifier( scip, "consnames", (int) nameClasses.size(), getNConss() );

   /** set all class names and descriptions */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::stringstream text;
      classifier->setClassName( c, nameClasses[c].c_str() );
      text << "This class contains all constraints with name \"" << nameClasses[c] << "\".";
      classifier->setClassDescription( c, text.str().c_str() );
   }

   /** copy the constraint assignment information found in first step */
   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass( i, classForCons[i] );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier %s yields a classification with %d  different constraint classes \n", classifier->getName(), classifier->getNClasses() );

   return classifier;
}

/** returns a new constraint classifier
 *  where all constraints whose consnames do not a have levenshtein distance to each other
 *  higher than a given connectivity are assigned to the same class */
ConsClassifier* Seeedpool::createConsClassifierForConsnamesLevenshteinDistanceConnectivity(
   int connectivity
   )
{
   std::vector < std::string > consnamesToCompare( getNConss(), "" );
   std::vector<int> nConssConstype( 0 );
   std::vector<int> classForCons = std::vector<int>( getNConss(), - 1 );
   std::vector<bool> alreadyReached( getNConss(), false );
   std::queue<int> helpqueue = std::queue<int>();
   int nUnreachedConss = getNConss();
   int currentClass = - 1;
   int nmaxconss = 5000;

   std::stringstream classifierName;
   classifierName << "lev-dist-" << connectivity;
   ConsClassifier* classifier = new ConsClassifier( scip, classifierName.str().c_str(), 0, getNConss() );

   /** if number of conss exceeds this number, skip calculating such a classifier */
   if( getNConss() > nmaxconss )
   {

      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " skipped levenshtein distance based constraint classes calculating since number of constraints  %d  exceeds limit %d \n", getNConss(), nmaxconss );
      delete classifier;
      return NULL;
   }

   std::vector < std::vector<int> > levenshteindistances( getNConss(), std::vector<int>( getNConss(), - 1 ) );

   /** read consnames */
   for( int i = 0; i < getNConss(); ++ i )
   {
      consnamesToCompare[i] = std::string( SCIPconsGetName( getConsForIndex( i ) ) );
   }

   /** calculate levenshtein distances pairwise */
   for( int i = 0; i < getNConss(); ++ i )
   {
      for( int j = i + 1; j < getNConss(); ++ j )
      {
         levenshteindistances[i][j] = calcLevenshteinDistance( consnamesToCompare[i], consnamesToCompare[j] );
         levenshteindistances[j][i] = levenshteindistances[i][j];
      }
   }

   /** repeat doing breadth first search until every constraint is assigned to a class */
   while( nUnreachedConss > 0 )
   {
      int firstUnreached = - 1;
      currentClass ++;
      assert( helpqueue.empty() );
      for( int i = 0; i < getNConss(); ++ i )
      {
         if( classForCons[i] == - 1 )
         {
            firstUnreached = i;
            break;
         }
      }

      helpqueue.push( firstUnreached );
      alreadyReached[firstUnreached] = true;
      classForCons[firstUnreached] = currentClass;
      -- nUnreachedConss;

      /** consider all constraints which are connected to the current constraint by means of levenshtein distance */
      while( ! helpqueue.empty() )
      {
         int nodecons = helpqueue.front();
         helpqueue.pop();
         for( int j = 0; j < getNConss(); ++ j )
         {

            if( alreadyReached[j] )
               continue;

            if( j == nodecons )
               continue;

            if( levenshteindistances[j][nodecons] > connectivity )
               continue;

            alreadyReached[j] = true;
            classForCons[j] = currentClass;
            -- nUnreachedConss;
            helpqueue.push( j );
         }
      }

      /** create a new class with found constraints in ConsClassifier*/
      std::stringstream text;
      text << "This class contains all constraints with a name similar to \"" << consnamesToCompare[firstUnreached] << "\".";
      classifier->addClass( consnamesToCompare[firstUnreached].c_str(), text.str().c_str(), BOTH );
   }

   /** assign constraint indices to classes */
   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass( i, classForCons[i] );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " consclassifier levenshtein: connectivity of %d yields a classification with %d different constraint classes. \n", connectivity, currentClass + 1);

   return classifier;
}

/** returns a new constraint classifier
 *  where all constraints with identical number of nonzero coefficients are assigned to the same class */
ConsClassifier* Seeedpool::createConsClassifierForNNonzeros()
{
   std::vector<int> nconssforclass( 0 );
   std::vector<int> differentNNonzeros( 0 );
   std::vector<int> classForCons( getNConss(), - 1 );
   int counterClasses = 0;

   /** firstly, assign all constraints to classindices */
   for( int i = 0; i < getNConss(); ++ i )
   {
      int consnnonzeros = getNVarsForCons( i );
      bool nzalreadyfound = false;

      /** check if number of nonzeros belongs to an existing class index */
      for( size_t nzid = 0; nzid < differentNNonzeros.size(); ++ nzid )
      {
         if( consnnonzeros == differentNNonzeros[nzid] )
         {
            nzalreadyfound = true;
            classForCons[i] = nzid;
            ++ nconssforclass[nzid];
            break;
         }
      }

      /** if not, create a new class index */
      if( ! nzalreadyfound )
      {
         classForCons[i] = counterClasses;
         ++ counterClasses;
         differentNNonzeros.push_back( consnnonzeros );
         nconssforclass.push_back( 1 );
      }
   }

   /** secondly, use these information to create a ConsClassifier */
   ConsClassifier* classifier = new ConsClassifier( scip, "nonzeros", (int) differentNNonzeros.size(), getNConss() );

   /** set class names and descriptions of every class */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::stringstream text;
      text << differentNNonzeros[c];
      classifier->setClassName( c, text.str().c_str() );
      text.str( "" );
      text.clear();
      text << "This class contains all constraints with " << differentNNonzeros[c] << " nonzero coefficients.";
      classifier->setClassDescription( c, text.str().c_str() );
   }

   /** copy the constraint assignment information found in first step */
   for( int i = 0; i < classifier->getNConss(); ++ i )
   {
      classifier->assignConsToClass( i, classForCons[i] );
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Consclassifier %s yields a classification with %d  different constraint classes \n", classifier->getName(), classifier->getNClasses() );

   return classifier;
}

/** returns pointer to a constraint classifier */
ConsClassifier* Seeedpool::getConsClassifier(
   int givenClassifierIndex
   )
{
   assert( 0 <= givenClassifierIndex && givenClassifierIndex < (int) consclassescollection.size() );

   return consclassescollection[givenClassifierIndex];
}

/** returns the assignment of constraints to classes of a classifier as integer array */
int* Seeedpool::getConsClassifierArray(
   int givenClassifierIndex
   )
{
   int nconss = consclassescollection[givenClassifierIndex]->getNConss();
   int* output = new int[nconss];
   for( int i = 0; i < nconss; ++ i )
      output[i] = consclassescollection[givenClassifierIndex]->getClassOfCons( i );
   return & output[0];
}

/** returns number of different constraint classifiers */
int Seeedpool::getNConsClassifiers()
{
   return (int) consclassescollection.size();
}

/** adds constraint classifiers with a reduced number of classes */
void Seeedpool::reduceConsclasses()
{
   /** set the number of classes the classifiers should be reduced to */
   int maxnclasses = 0;

   if( getNConss() + getNVars() >= 50000 )
      SCIPgetIntParam(scip, "detection/maxnclassesperclassifierforlargeprobs", &maxnclasses);
   else
      SCIPgetIntParam(scip, "detection/maxnclassesperclassifier", &maxnclasses);

   for( size_t classifierid = 0; classifierid < consclassescollection.size(); ++ classifierid )
   {
      ConsClassifier* newclassifier = consclassescollection[classifierid]->reduceClasses( maxnclasses );

      if( newclassifier != NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Added reduced version of consclassifier %s with %d  different constraint classes \n", consclassescollection[classifierid]->getName(), maxnclasses );
         addConsClassifier( newclassifier );
      }
   }
}

/** adds a variable classifier if it is no duplicate of an existing variable classifier */
void Seeedpool::addVarClassifier(
   VarClassifier* givenClassifier
   )
{
   SCIP_Bool detectionstatistics;

   SCIPgetBoolParam(scip, "detection/allowclassifierduplicates/enabled", &detectionstatistics);

   if( givenClassifier != NULL )
   {
      /** check whether there already exists an equivalent varclassifier */
      VarClassifier* equiv = NULL;

      for( size_t i = 0; !detectionstatistics && i < varclassescollection.size(); ++ i )
      {
         if( givenClassifier->classifierIsDuplicateOfClassifier( varclassescollection[i] ) )
         {
            equiv = varclassescollection[i];
            break;
         }
      }

      if( equiv == NULL )
         varclassescollection.push_back( givenClassifier );
      else
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier %s not considered since it offers the same structure as  %s.\n", givenClassifier->getName(), equiv->getName() );
         delete givenClassifier;
      }

   }
}

/** returns a new variable classifier
 *  where all variables with identical objective function value are assigned to the same class */
VarClassifier* Seeedpool::createVarClassifierForObjValues()
{
   std::vector<SCIP_Real> foundobjvals( 0 ); /** all found objective function values */
   std::vector<int> classforvars( getNVars(), -1 ); /** vector assigning a class index to each variable */
   int curclassindex; /** stores a var's classindex if the objective value of a var has already been found for another var */
   SCIP_Real curobjval;
   VarClassifier* classifier; /** new VarClassifier */

   for( int v = 0; v < getNVars(); ++v )
   {
      assert( getVarForIndex( v ) != NULL );
      curobjval = SCIPvarGetObj( getVarForIndex( v ) );
      curclassindex = -1;

      /** check whether current objective funtion value already exists */
      for( size_t c = 0; c < foundobjvals.size(); ++c )
      {
         if( SCIPisEQ( scip, curobjval, foundobjvals[c] ) )
         {
            curclassindex = c;
            break;
         }
      }

      /** assign var to class and save objective function value, if it is new */
      if( curclassindex == -1 )
      {
         foundobjvals.push_back( curobjval );
         classforvars[v] = foundobjvals.size() - 1;
      }
      else
      {
         classforvars[v] = curclassindex;
      }
   }

   classifier = new VarClassifier( scip, "varobjvals", (int) foundobjvals.size(), getNVars() );

   /** set up class information */
   for ( int c = 0; c < classifier->getNClasses(); ++c )
   {
      std::stringstream name;
      std::stringstream text;

      name << std::setprecision( 5 ) << foundobjvals[c];
      text << "This class contains all variables with objective function value " << name.str() << ".";

      classifier->setClassName( c, name.str().c_str() );
      classifier->setClassDescription( c, text.str().c_str() );
   }

   /** assign vars according to classforvars vactor */
   for ( int v = 0; v < classifier->getNVars(); ++v )
   {
      classifier->assignVarToClass( v, classforvars[v] );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier %s yields a classification with %d different variable classes.\n", classifier->getName(), classifier->getNClasses() ) ;

   return classifier;
}

/** returns a new variable classifier
 *  where all variables are assigned to class zero, positive or negative according to their objective function value sign
 *  all class zero variables are assumed to be only master variables (set via DECOMPINFO)
 *  @todo correct? */
VarClassifier* Seeedpool::createVarClassifierForObjValueSigns()
{
   VarClassifier* classifier= new VarClassifier( scip, "varobjvalsigns", 3, getNVars() ); /** new VarClassifier */
   SCIP_Real curobjval;

   /** set up class information */
   classifier->setClassName( 0, "zero" );
   classifier->setClassDescription( 0, "This class contains all variables with objective function value zero." );
   classifier->setClassDecompInfo( 0, MASTER );
   classifier->setClassName( 1, "positive" );
   classifier->setClassDescription( 1, "This class contains all variables with positive objective function value." );
   classifier->setClassDecompInfo( 1, ALL );
   classifier->setClassName( 2, "negative" );
   classifier->setClassDescription( 2, "This class contains all variables with negative objective function value." );
   classifier->setClassDecompInfo( 2, ALL );

   /** assign vars */
   for( int v = 0; v < getNVars(); ++v )
   {
      assert( getVarForIndex( v ) != NULL );
      curobjval = SCIPvarGetObj( getVarForIndex( v ) );

      if( SCIPisZero( scip, curobjval ) )
      {
         classifier->assignVarToClass( v, 0 );
      }
      else if ( SCIPisPositive( scip, curobjval ) )
      {
         classifier->assignVarToClass( v, 1 );
      }
      else
      {
         classifier->assignVarToClass( v, 2 );
      }
   }

   /* remove a class if there is no variable with the respective sign */
   classifier->removeEmptyClasses();

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier %s yields a classification with %d different variable classes.\n", classifier->getName(), classifier->getNClasses() ) ;

   return classifier;
}

/** returns a new variable classifier
 *  where all variables with identical SCIP vartype are assigned to the same class */
VarClassifier* Seeedpool::createVarClassifierForSCIPVartypes()
{
   std::vector < SCIP_VARTYPE > foundVartypes( 0 );
   std::vector<int> classForVars = std::vector<int>( getNVars(), - 1 );
   VarClassifier* classifier;

   /** firstly, assign all variables to classindices */
   for( int i = 0; i < getNVars(); ++ i )
   {
      SCIP_VAR* var;
      bool found = false;
      var = getVarForIndex( i );
      SCIP_VARTYPE vT = SCIPvarGetType( var );
      size_t vartype;

      /** check whether the variable's vartype is new */
      for( vartype = 0; vartype < foundVartypes.size(); ++ vartype )
      {
         if( foundVartypes[vartype] == vT )
         {
            found = true;
            break;
         }
      }
      /** if it is new, create a new class index */
      if( ! found )
      {
         foundVartypes.push_back( vT );
         classForVars[i] = foundVartypes.size() - 1;
      }
      else
         classForVars[i] = vartype;
   }

   /** secondly, use these information to create a VarClassifier */
   classifier = new VarClassifier( scip, "vartypes", (int) foundVartypes.size(), getNVars() );

   /** set class names and descriptions of every class */
   for( int c = 0; c < classifier->getNClasses(); ++ c )
   {
      std::string name;
      std::stringstream text;
      switch( foundVartypes[c] )
      {
         case SCIP_VARTYPE_BINARY:
            name = "bin";
            break;
         case SCIP_VARTYPE_INTEGER:
            name = "int";
            break;
         case SCIP_VARTYPE_IMPLINT:
            name = "impl";
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            name = "cont";
            break;
         default:
            name = "newVartype";
            break;
      }
      classifier->setClassName( c, name.c_str() );
      text << "This class contains all variables that are of (SCIP) vartype \"" << name << "\".";
      classifier->setClassDescription( c, text.str().c_str() );
   }

   /** copy the variable assignment information found in first step */
   for( int i = 0; i < classifier->getNVars(); ++ i )
   {
      classifier->assignVarToClass( i, classForVars[i] );
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Varclassifier %s yields a classification with %d different variable classes.\n", classifier->getName(), classifier->getNClasses() ) ;

   return classifier;
}

/** returns pointer to a variable classifier */
VarClassifier* Seeedpool::getVarClassifier(
   int givenClassifierIndex
   )
{
   assert( 0 <= givenClassifierIndex && givenClassifierIndex < (int) varclassescollection.size() );

   return varclassescollection[givenClassifierIndex];
}

/** returns the assignment of variables to classes of a classifier as integer array */
int* Seeedpool::getVarClassifierArray(
   int givenClassifierIndex
   )
{
   int nvars = varclassescollection[givenClassifierIndex]->getNVars();
   int* output = new int[nvars];
   for( int i = 0; i < nvars; ++ i )
      output[i] = varclassescollection[givenClassifierIndex]->getClassOfVar( i );
   return & output[0];
}

/** returns number of different variable classifiers */
int Seeedpool::getNVarClassifiers()
{
   return (int) varclassescollection.size();
}

/** adds variable classifiers with a reduced number of classes */
void Seeedpool::reduceVarclasses()
{
   /** set the number of classes the classifiers should be reduced to */
   int maxnclasses = 0;

   if( getNConss() + getNVars() >= 50000 )
      SCIPgetIntParam(scip, "detection/maxnclassesperclassifierforlargeprobs", &maxnclasses);
   else
      SCIPgetIntParam(scip, "detection/maxnclassesperclassifier", &maxnclasses);

   for( size_t classifierid = 0; classifierid < varclassescollection.size(); ++ classifierid )
   {
      VarClassifier* newclassifier = varclassescollection[classifierid]->reduceClasses( maxnclasses );

      if( newclassifier != NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Added reduced version of varclassifier %s with %d different variable classes.\n", varclassescollection[classifierid]->getName(), maxnclasses ) ;
         addVarClassifier( newclassifier );
      }
   }
}

/** returns a vector of seeeds where all seeeds of the given seeeds having only one block are removed
 *  except for the two seeeds with the lowest numbers of masterconss */
std::vector<SeeedPtr> Seeedpool::removeSomeOneblockDecomps(
   std::vector<SeeedPtr> seeeds
   )
{
   std::vector<SeeedPtr> remainingSeeeds( 0 );
   std::vector<SeeedPtr> oneBlockSeeeds( 0 );

   int nmasterconssfirst = 1000;
   int nmasterconsssecond = 1001;

   for( size_t i = 0; i < seeeds.size(); ++ i )
   {
      /** calculate lowest and second lowest number of masterconss of all one block seeeds */
      if( seeeds[i]->getNBlocks() == 1 )
      {
         if( seeeds[i]->getNMasterconss() < nmasterconssfirst )
         {
            nmasterconsssecond = nmasterconssfirst;
            nmasterconssfirst = seeeds[i]->getNMasterconss();
         }
         else if( seeeds[i]->getNMasterconss() < nmasterconsssecond )
            nmasterconsssecond = seeeds[i]->getNMasterconss();

      }
      else
         remainingSeeeds.push_back( seeeds[i] );
   }

   /** the two one block seeeds with lowest number of masterconss remain */
   for( int i = 0; i < (int) seeeds.size(); ++ i )
   {
      if( seeeds[i]->getNBlocks() == 1
         && ( seeeds[i]->getNMasterconss() == nmasterconssfirst || seeeds[i]->getNMasterconss() == nmasterconsssecond ) )
         remainingSeeeds.push_back( seeeds[i] );
      else if( seeeds[i]->getNBlocks() == 1 )
         oneBlockSeeeds.push_back( seeeds[i] );
   }

   /** all other one block seeeds are removed */
   for( int i = 0; i < (int) oneBlockSeeeds.size(); ++ i )
   {
      delete oneBlockSeeeds[i];
      oneBlockSeeeds[i] = NULL;
   }

   return remainingSeeeds;
}

/** creates a decomposition for a given seeed */
SCIP_RETCODE Seeedpool::createDecompFromSeeed(
   SeeedPtr seeed,
   DEC_DECOMP** newdecomp
   )
{
   char detectorchaininfo[SCIP_MAXSTRLEN];
   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;
   SCIP_HASHMAP* varindex;
   SCIP_HASHMAP* consindex;
   SCIP_VAR*** stairlinkingvars;
   SCIP_VAR*** subscipvars;
   SCIP_VAR** linkingvars;
   SCIP_CONS** linkingconss;
   SCIP_CONS*** subscipconss;
   int* nsubscipconss;
   int* nsubscipvars;
   int* nstairlinkingvars;
   int nlinkingvars;
   int varcounter = 1; /* in varindex counting starts with 1 */
   int conscounter = 1; /* in consindex counting starts with 1 */
   int counterstairlinkingvars = 0;
   int size;
   int modifier;
   int nlinkingconss;
   assert( seeed->checkConsistency( this ) );




   /* create decomp data structure */
   SCIP_CALL_ABORT( DECdecompCreate( scip, newdecomp ) );

   /** set nblocks */
   DECdecompSetNBlocks( * newdecomp, seeed->getNBlocks() );

   //detectorchaininfo ;
   /** set constraints */
   if( seeed->getNMasterconss() != 0 )
      SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & linkingconss, seeed->getNMasterconss() ) );
   else
      linkingconss = NULL;

   SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & nsubscipconss, seeed->getNBlocks() ) );
   SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & subscipconss, seeed->getNBlocks() ) );

   SCIP_CALL_ABORT( SCIPhashmapCreate( & constoblock, SCIPblkmem( scip ), seeed->getNConss() ) );
   SCIP_CALL_ABORT( SCIPhashmapCreate( & consindex, SCIPblkmem( scip ), seeed->getNConss() ) );

   /* set linking constraints */
   modifier = 0;
   nlinkingconss = seeed->getNMasterconss();
   for( int c = 0; c < seeed->getNMasterconss(); ++ c )
   {
      int consid = seeed->getMasterconss()[c];
      SCIP_CONS* scipcons = consToScipCons[consid];
      if( SCIPconsIsDeleted( scipcons) || scipcons == NULL || SCIPconsIsObsolete(scipcons))
      {
         --nlinkingconss;
         ++modifier;
      }
      else
      {
         linkingconss[c-modifier] = scipcons;
         SCIP_CALL_ABORT( SCIPhashmapInsert( constoblock, scipcons, (void*) ( size_t )( seeed->getNBlocks() + 1 ) ) );
         SCIP_CALL_ABORT( SCIPhashmapInsert( consindex, scipcons, (void*) (size_t) conscounter ) );
         conscounter ++;
      }
   }

   if( nlinkingconss != 0 )
      DECdecompSetLinkingconss( scip, * newdecomp, linkingconss, nlinkingconss );
   else
      linkingconss = NULL;
   /* set block constraints */
   for( int b = 0; b < seeed->getNBlocks(); ++ b )
   {
      modifier = 0;
      SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & subscipconss[b], seeed->getNConssForBlock( b ) ) );
      nsubscipconss[b] = seeed->getNConssForBlock( b );
      for( int c = 0; c < seeed->getNConssForBlock( b ); ++ c )
      {
         int consid = seeed->getConssForBlock( b )[c];
         SCIP_CONS* scipcons = consToScipCons[consid];

         if( SCIPconsIsDeleted( scipcons) || scipcons == NULL )
            {
               --nsubscipconss[b];
               ++modifier;
            }
         else
         {
            assert( scipcons != NULL );
            subscipconss[b][c-modifier] = scipcons;
            SCIP_CALL_ABORT( SCIPhashmapInsert( constoblock, scipcons, (void*) ( size_t )( b + 1 ) ) );
            SCIP_CALL_ABORT( SCIPhashmapInsert( consindex, scipcons, (void*) (size_t) conscounter ) );
            conscounter ++;
         }
      }
   }

   DECdecompSetSubscipconss( scip, * newdecomp, subscipconss, nsubscipconss );

   DECdecompSetConstoblock( * newdecomp, constoblock );
   DECdecompSetConsindex( * newdecomp, consindex );

   /* finished setting constraint data structures */
   /** now: set variables */
   SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & nsubscipvars, seeed->getNBlocks() ) );
   SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & subscipvars, seeed->getNBlocks() ) );
   SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & nstairlinkingvars, seeed->getNBlocks() ) );
   SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & stairlinkingvars, seeed->getNBlocks() ) );

   SCIP_CALL_ABORT( SCIPhashmapCreate( & vartoblock, SCIPblkmem( scip ), seeed->getNVars() ) );
   SCIP_CALL_ABORT( SCIPhashmapCreate( & varindex, SCIPblkmem( scip ), seeed->getNVars() ) );

   /** set linkingvars */
   nlinkingvars = seeed->getNLinkingvars() + seeed->getNMastervars() + seeed->getNTotalStairlinkingvars();

   if( nlinkingvars != 0 )
      SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & linkingvars, nlinkingvars ) );
   else
      linkingvars = NULL;

   for( int v = 0; v < seeed->getNLinkingvars(); ++ v )
   {
      int var = seeed->getLinkingvars()[v];
      SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
      assert( scipvar != NULL );

      linkingvars[v] = scipvar;
      SCIP_CALL_ABORT( SCIPhashmapInsert( vartoblock, scipvar, (void*) ( size_t )( seeed->getNBlocks() + 2 ) ) );
      SCIP_CALL_ABORT( SCIPhashmapInsert( varindex, scipvar, (void*) (size_t) varcounter ) );
      varcounter ++;
   }

   for( int v = 0; v < seeed->getNMastervars(); ++ v )
   {
      int var = seeed->getMastervars()[v];
      SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
      linkingvars[v + seeed->getNLinkingvars()] = scipvar;
      SCIP_CALL_ABORT( SCIPhashmapInsert( vartoblock, scipvar, (void*) ( size_t )( seeed->getNBlocks() + 1 ) ) );
      SCIP_CALL_ABORT( SCIPhashmapInsert( varindex, scipvar, (void*) (size_t) varcounter ) );
      varcounter ++;
   }

   /* set block variables */
   for( int b = 0; b < seeed->getNBlocks(); ++ b )
   {
      if( seeed->getNVarsForBlock( b ) > 0 )
         SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & subscipvars[b], seeed->getNVarsForBlock( b ) ) );
      else
         subscipvars[b] = NULL;

      if( seeed->getNStairlinkingvars( b ) > 0 )
         SCIP_CALL_ABORT( SCIPallocBufferArray( scip, & stairlinkingvars[b], seeed->getNStairlinkingvars( b ) ) );
      else
         stairlinkingvars[b] = NULL;

      nsubscipvars[b] = seeed->getNVarsForBlock( b );
      nstairlinkingvars[b] = seeed->getNStairlinkingvars( b );

      for( int v = 0; v < seeed->getNVarsForBlock( b ); ++ v )
      {
         int var = seeed->getVarsForBlock( b )[v];
         SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
         assert( scipvar != NULL );

         subscipvars[b][v] = scipvar;
         SCIP_CALL_ABORT( SCIPhashmapInsert( vartoblock, scipvar, (void*) ( size_t )( b + 1 ) ) );
         SCIP_CALL_ABORT( SCIPhashmapInsert( varindex, scipvar, (void*) (size_t) varcounter ) );
         varcounter ++;
      }

      for( int v = 0; v < seeed->getNStairlinkingvars( b ); ++ v )
      {
         int var = seeed->getStairlinkingvars( b )[v];
         SCIP_VAR* scipvar = SCIPvarGetProbvar( varToScipVar[var] );
         assert( scipvar != NULL );

         stairlinkingvars[b][v] = scipvar;
         linkingvars[seeed->getNLinkingvars() + seeed->getNMastervars() + counterstairlinkingvars] = scipvar;
         SCIP_CALL_ABORT( SCIPhashmapInsert( vartoblock, scipvar, (void*) ( size_t )( seeed->getNBlocks() + 2 ) ) );
         SCIP_CALL_ABORT( SCIPhashmapInsert( varindex, scipvar, (void*) (size_t) varcounter ) );
         varcounter ++;
         counterstairlinkingvars ++;
      }
   }

   DECdecompSetSubscipvars( scip, * newdecomp, subscipvars, nsubscipvars );
   DECdecompSetStairlinkingvars( scip, * newdecomp, stairlinkingvars, nstairlinkingvars );
   DECdecompSetLinkingvars( scip, * newdecomp, linkingvars, nlinkingvars, seeed->getNMastervars() );
   DECdecompSetVarindex( * newdecomp, varindex );
   DECdecompSetVartoblock( * newdecomp, vartoblock );

   /** free stuff */


   /** free vars stuff */
   SCIPfreeBufferArrayNull( scip, & ( linkingvars ) );
   for( int b = seeed->getNBlocks() - 1; b >= 0; -- b )
   {
      if( nstairlinkingvars[b] != 0 )
      {
         SCIPfreeBufferArrayNull( scip, & ( stairlinkingvars[b] ) );
      }
   }


   SCIPfreeBufferArrayNull( scip, & ( stairlinkingvars ) );
   SCIPfreeBufferArrayNull( scip, & ( nstairlinkingvars ) );

   for( int b = seeed->getNBlocks() - 1; b >= 0; -- b )
   {
      if( nsubscipvars[b] != 0 )
      {
         SCIPfreeBufferArrayNull( scip, & ( subscipvars[b] ) );
      }
   }

   SCIPfreeBufferArrayNull( scip, & ( subscipvars ) );
   SCIPfreeBufferArrayNull( scip, & ( nsubscipvars ) );


   /** free constraints */
   for( int b = seeed->getNBlocks() - 1; b >= 0; -- b )
   {
      SCIPfreeBufferArrayNull( scip, & ( subscipconss[b] ) );
   }
   SCIPfreeBufferArrayNull( scip, & ( subscipconss ) );


   SCIPfreeBufferArrayNull( scip, & nsubscipconss );
   SCIPfreeBufferArrayNull( scip, & linkingconss );


   /** set detectorchain */
   int ndetectors = seeed->getNDetectors();
   ( * newdecomp )->sizedetectorchain = ndetectors;
   size = SCIPcalcMemGrowSize( scip, ( * newdecomp )->sizedetectorchain );
   SCIP_CALL_ABORT( SCIPallocBlockMemoryArray( scip, & ( * newdecomp )->detectorchain, size ) );
   for( int k = 0; k < ndetectors; ++ k )
   {
      if( k != ndetectors - 1 || ! seeed->getFinishedByFinisher() )
      {
         //          std::cout << " added detector of " << i << "-th seeed to its detetcor chain" << std::endl;
         ( * newdecomp )->detectorchain[k] = seeed->getDetectorchain()[k];
      }
      else
         ( * newdecomp )->detectorchain[k] = seeed->getDetectorchain()[k];
   }

   /** set statistical detector chain data */
   DECdecompSetSeeedID( * newdecomp, seeed->getID() );
   if( seeed->getNDetectors() > 0 )
   {
      DECdecompSetDetectorClockTimes( scip, * newdecomp, & ( seeed->getDetectorClockTimes()[0] ) );
      DECdecompSetDetectorPctVarsToBorder( scip, * newdecomp, & ( seeed->getPctVarsToBorderVector()[0] ) );
      DECdecompSetDetectorPctVarsToBlock( scip, * newdecomp, & ( seeed->getPctVarsToBlockVector()[0] ) );
      DECdecompSetDetectorPctVarsFromOpen( scip, * newdecomp, & ( seeed->getPctVarsFromFreeVector()[0] ) );
      DECdecompSetDetectorPctConssToBorder( scip, * newdecomp, & ( seeed->getPctConssToBorderVector()[0] ) );
      DECdecompSetDetectorPctConssToBlock( scip, * newdecomp, & ( seeed->getPctConssToBlockVector()[0] ) );
      DECdecompSetDetectorPctConssFromOpen( scip, * newdecomp, & ( seeed->getPctConssFromFreeVector()[0] ) );
      DECdecompSetNNewBlocks( scip, * newdecomp, & ( seeed->getNNewBlocksVector()[0] ) );
   }


   /** set detector chain info string */
   SCIPsnprintf( detectorchaininfo, SCIP_MAXSTRLEN, "") ;
   if( seeed->getUsergiven() == USERGIVEN::PARTIAL || seeed->getUsergiven() == USERGIVEN::COMPLETE || seeed->getUsergiven() == USERGIVEN::COMPLETED_CONSTOMASTER)
   {
      seeed->buildDecChainString();
   }
   SCIP_CALL( DECdecompSetDetectorChainString( scip, * newdecomp, seeed->getDetectorChainString() ) );

   /** set dectype */
   if( ( * newdecomp )->nlinkingvars == seeed->getNTotalStairlinkingvars() && ( * newdecomp )->nlinkingconss == 0
      && DECdecompGetNLinkingvars( * newdecomp ) > 0 )
   {
      ( * newdecomp )->type = DEC_DECTYPE_STAIRCASE;
   }
   else if( ( * newdecomp )->nlinkingvars > 0 || seeed->getNTotalStairlinkingvars() > 0 )
   {
      ( * newdecomp )->type = DEC_DECTYPE_ARROWHEAD;
   }
   else if( ( * newdecomp )->nlinkingconss > 0 )
   {
      ( * newdecomp )->type = DEC_DECTYPE_BORDERED;
   }
   else if( ( * newdecomp )->nlinkingconss == 0 && seeed->getNTotalStairlinkingvars() == 0 )
   {
      ( * newdecomp )->type = DEC_DECTYPE_DIAGONAL;
   }
   else
   {
      ( * newdecomp )->type = DEC_DECTYPE_UNKNOWN;
   }


   //SCIP_CALL( DECevaluateDecomposition( scip, * newdecomp, & scores ) );

   std::cout <<" seeed maxwhitescore: " << seeed->getMaxWhiteScore() << std::endl;

   DECsetMaxWhiteScore(scip, *newdecomp, seeed->getMaxWhiteScore() );


   assert( DECdecompCheckConsistency( scip, ( * newdecomp ) ) );
   assert( ! SCIPhashmapIsEmpty( ( * newdecomp )->constoblock ) );
   assert( ! SCIPhashmapIsEmpty( ( * newdecomp )->vartoblock ) );

   return SCIP_OKAY;
}

/** creates a seeed for a given decomposition
 *  the resulting seeed will not have a detectorchaininfo or any ancestor or finishing detector data
 *  only use this method if the seeedpool is for the transformed problem
 *  the resulting seeed may only be added to the seeedpool for the presolved problem */
SCIP_RETCODE Seeedpool::createSeeedFromDecomp(
   DEC_DECOMP* decomp,
   SeeedPtr* newseeed
   )
{
   assert( decomp != NULL );
   assert( DECdecompCheckConsistency( scip, decomp ) );
   assert( nConss == DECdecompGetNConss( decomp ) );
   assert( DECdecompGetPresolved( decomp ) );
   assert( transformed );

//   std::cout << "Linkingvars decomp: " << DECdecompGetNLinkingvars( decomp ) << "\tStairlinkingvars decomp: " << DECdecompGetNTotalStairlinkingvars( decomp ) << "\n";

   /* create new seeed and initialize its data */
   SeeedPtr seeed = new Seeed( scip, getNewIdForSeeed(), nConss, nVars );
   seeed->setNBlocks( DECdecompGetNBlocks( decomp ) );

   assert( seeed->getNOpenconss() == nConss );
   assert( seeed->getNOpenvars() == nVars );

   SCIP_CONS** linkingconss = DECdecompGetLinkingconss( decomp );
   int nlinkingconss = DECdecompGetNLinkingconss( decomp );
   SCIP_HASHMAP* constoblock = DECdecompGetConstoblock( decomp );
   int nblock;

   /* set linking conss */
   for( int c = 0; c < nlinkingconss; ++c )
   {
      seeed->bookAsMasterCons( getIndexForCons( linkingconss[c] ) );
   }

   /* set block conss */
   for( int c = 0; c < nConss; ++c )
   {
      nblock = (int) (size_t) SCIPhashmapGetImage( constoblock, (void*) (size_t) getConsForIndex( c ) );
      if( nblock >= 1 && nblock <= seeed->getNBlocks() )
      {
         seeed->bookAsBlockCons( c, nblock - 1 );
      }
   }

   SCIP_VAR*** stairlinkingvars = DECdecompGetStairlinkingvars( decomp );
   SCIP_HASHMAP* vartoblock = DECdecompGetVartoblock(decomp);
   assert( vartoblock != NULL );

   /* currently, stairlinkingvars of the decomposition are ignored
    * alternatively, a stairlinkingvar detection is done with the newly created seeed */
   if( false && stairlinkingvars != NULL )
   {
      int* nstairlinkingvars = DECdecompGetNStairlinkingvars(decomp);
      int varindex;

      /* set stairlinkingvars */
      for( int b = 0; b < seeed->getNBlocks(); ++b )
      {
         for( int v = 0; v < nstairlinkingvars[b]; ++v )
         {
            if( stairlinkingvars[b][v] != NULL )
            {
               varindex = getIndexForVar(stairlinkingvars[b][v]);
               seeed->bookAsStairlinkingVar(varindex, b);
            }
         }
      }
   }

   /* flush booked conss and vars in order to be able to check whether a var is already assigned to stairlinking */
   seeed->flushBooked();

   /* set other vars */
   for( int v = 0; v < getNVars(); ++v )
   {
      nblock = (int) (size_t) SCIPhashmapGetImage( vartoblock, (void*) (size_t) SCIPvarGetProbvar( getVarForIndex( v ) ) );
      if( nblock == seeed->getNBlocks() + 2 && !seeed->isVarStairlinkingvar( v ) )
      {
         seeed->bookAsLinkingVar( v );
      }
      else if( nblock == seeed->getNBlocks() + 1 )
      {
         seeed->bookAsMasterVar( v );
      }
      else if( nblock >= 1 && nblock <= seeed->getNBlocks() )
      {
         seeed->bookAsBlockVar( v, nblock - 1 );
      }
   }

   seeed->flushBooked();

   /* now all conss and vars should be assigned */
   assert( seeed->isComplete() );
   /*set all detector-related information*/
   for( int i = 0; i < DECdecompGetDetectorChainSize( decomp ); ++i )
   {
      seeed->setDetectorPropagated( DECdecompGetDetectorChain( decomp )[i] );
      seeed->addClockTime( DECdecompGetDetectorClockTimes( decomp )[i] );
      seeed->addPctConssFromFree( 1 - *(DECdecompGetDetectorPctConssFromOpen( decomp )) );
      seeed->addPctConssToBlock( *(DECdecompGetDetectorPctConssToBlock( decomp )) );
      seeed->addPctConssToBorder( *(DECdecompGetDetectorPctConssToBorder( decomp )) );
      seeed->addPctVarsFromFree( 1 - *(DECdecompGetDetectorPctVarsFromOpen( decomp )) );
      seeed->addPctVarsToBlock( *(DECdecompGetDetectorPctVarsToBlock( decomp )) );
      seeed->addPctVarsToBorder( *(DECdecompGetDetectorPctVarsToBorder( decomp )) );
      seeed->addNNewBlocks( *(DECdecompGetNNewBlocks( decomp )) );
   }

   if ( DECdecompGetDetectorChainString( scip, decomp ) != NULL )
      seeed->setDetectorChainString( DECdecompGetDetectorChainString( scip, decomp ) );

   /* detectorchaininfo cannot be set in the seeed as the detectors do not store the corresponding strings */

   /* calc maxwhitescore and hashvalue */
   prepareSeeed( seeed );

   seeed->setIsFromUnpresolved( false );

//   SCIPdebugMessagePrint(scip, "Check. DEC: %f, seeed: %f\n", DECgetMaxWhiteScore( scip, decomp ), seeed->getMaxWhiteScore() );

//   assert( DECgetMaxWhiteScore( scip, decomp ) == seeed->getMaxWhiteScore() );

   assert( seeed->checkConsistency( this ) );

   seeed->calcStairlinkingVars( this );

//   SCIPdebugMessagePrint( scip, "Reassigned %d of %d linking vars to stairlinking.\n",
//      seeed->getNTotalStairlinkingvars(), seeed->getNTotalStairlinkingvars() + seeed->getNLinkingvars() );

   *newseeed = seeed;

   return SCIP_OKAY;
}

/** returns true if the matrix structure corresponds to the transformed problem */
SCIP_Bool Seeedpool::getTransformedInfo()
{
   return transformed;
}

SCIP_RETCODE Seeedpool::printBlockcandidateInformation(
 SCIP*                 givenscip,               /**< SCIP data structure */
 FILE*                 file                /**< output file or NULL for standard output */
)
{

   std::sort( candidatesNBlocks.begin(), candidatesNBlocks.end(), sort_decr() );
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d  \n", (int) candidatesNBlocks.size() );
   for( size_t i  = 0; i  < candidatesNBlocks.size(); ++i )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d : %d  \n", candidatesNBlocks[i].first, candidatesNBlocks[i].second );
   }

   return SCIP_OKAY;
}

SCIP_RETCODE Seeedpool::printClassifierInformation(
 SCIP*                 givenscip,               /**< SCIP data structure */
 FILE*                 file                /**< output file or NULL for standard output */
)
{

   /** NCLASSIFIER */
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d  \n", (int) consclassescollection.size()  );

   for( size_t c = 0; c < consclassescollection.size() ; ++c )
   {
      gcg::ConsClassifier* classifier = consclassescollection[c];

      std::vector<std::vector<int> > conssofclasses = std::vector<std::vector<int> >(classifier->getNClasses()) ;
      for( int cons = 0; cons < getNConss(); ++cons )
         conssofclasses[classifier->getClassOfCons(cons)].push_back(cons);

      /** CLASSIFIERNAME */
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%s  \n",  classifier->getName() );


      /** NCLASSES */
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d  \n",  classifier->getNClasses() );

      for( int cl = 0; cl < classifier->getNClasses(); ++cl )
      {
         /** CLASSNAME */
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%s: %s\n", classifier->getClassName(cl), classifier->getClassDescription(cl) );
         /** NMEMBERS */
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d\n",  conssofclasses[cl].size() );
      }
   }

   /** NCLASSIFIER */
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d  \n", (int) varclassescollection.size()  );

      for( size_t c = 0; c < varclassescollection.size() ; ++c )
      {
         gcg::VarClassifier* classifier = varclassescollection[c];

         std::vector<std::vector<int> > varsofclasses = std::vector<std::vector<int> >(classifier->getNClasses()) ;
         for( int var = 0; var < getNVars(); ++var )
            varsofclasses[classifier->getClassOfVar(var)].push_back(var);

         /** CLASSIFIERNAME */
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%s  \n",  classifier->getName() );


         /** NCLASSES */
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d  \n",  classifier->getNClasses() );

         for( int cl = 0; cl < classifier->getNClasses(); ++cl )
         {
            /** CLASSNAME */
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%s: %s\n", classifier->getClassName(cl), classifier->getClassDescription(cl) );
            /** NMEMBERS */
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%d\n",  varsofclasses[cl].size() );
         }
      }


   return SCIP_OKAY;
}




} /* namespace gcg */
