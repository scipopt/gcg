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

/**@file   class_seeed.cpp
 * @brief  class with functions for seeed (aka incomplete decomposition )
 * @author Michael Bastubbe
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*#define SCIP_DEBUG*/

#include "class_seeed.h"
#include "gcg.h"
#include "class_seeedpool.h"
#include "scip/cons_setppc.h"
#include "scip/scip.h"
#include "scip_misc.h"
#include "decomp.h"
#include "struct_detector.h"
#include "struct_decomp.h"
#include "cons_decomp.h"
#include "params_visu.h"
#include "class_miscvisualization.h"
#include "reader_gp.h"

#include <sstream>
#include <iostream>
#include <exception>
#include <algorithm>
#include <queue>
#include <fstream>
#include <stdlib.h>

#ifdef WITH_BLISS
#include "pub_bliss.h"
#include "bliss_automorph.h"
#endif


#define SCIP_CALL_EXC( x ) do                                                                                 \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( ( _restat_ = ( x ) ) !=  SCIP_OKAY )                                            \
                          {                                                                                   \
                             SCIPerrorMessage( "Error <%d> in function call\n", _restat_ );                   \
                             throw std::exception();                                                          \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

namespace gcg {

const int Seeed::primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
   101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227,
   229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349};

const int Seeed::nPrimes = 70;

/** constructor
 *  initially, all conss and vars are open */
Seeed::Seeed(
   SCIP* _scip,
   int givenid,
   Seeedpool* givenseeedpool
   ) :
   scip( _scip ), id( givenid ), nBlocks( 0 ), nVars( givenseeedpool->getNVars() ), nConss( givenseeedpool->getNConss() ), masterConss( 0 ),
   masterVars( 0 ), conssForBlocks( 0 ), varsForBlocks( 0 ), linkingVars( 0 ), stairlinkingVars( 0 ), isvaropen( givenseeedpool->getNVars(), true ),
   isconsopen( givenseeedpool->getNConss(), true ), isvarmaster( givenseeedpool->getNVars(), false ),   isconsmaster( givenseeedpool->getNConss(), false ),
   ncoeffsforblock(std::vector<int>(0)), calculatedncoeffsforblock(FALSE), ncoeffsforblockformastercons(0),
   varsforblocksorted(true), stairlinkingvarsforblocksorted(true),
   conssforblocksorted(true), linkingvarssorted(true), mastervarssorted(true),
   masterconsssorted(true), hashvalue( 0 ), changedHashvalue( false ), isselected( false ), isagginfoalreadytoexpensive(false), isFinishedByFinisher( false ),
   agginfocalculated(FALSE), nrepblocks(0), reptoblocks(std::vector<std::vector<int>>(0)), blockstorep(std::vector<int>(0) ), pidtopidvarmaptofirst(std::vector<std::vector<std::vector<int> > >(0)),
   detectorChain( 0 ), detectorChainFinishingUsed( 0 ), detectorClockTimes( 0 ), pctVarsToBorder( 0 ),
   pctVarsToBlock( 0 ), pctVarsFromFree( 0 ), pctConssToBorder( 0 ), pctConssToBlock( 0 ), pctConssFromFree( 0 ),
   nNewBlocks( 0 ), usedClassifier( 0 ), classesToMaster( 0 ), classesToLinking( 0 ), listofancestorids( 0 ),
   usergiven( USERGIVEN::NOT ), isfromlegacymode( false ), score( -1. ), maxwhitescore( -1. ), bendersscore(-1.), benderareascore(-1.),  strongdecompositionscore(-1.), borderareascore( -1. ),
   maxwhitescoreagg(-1.), blockareascore(-1.), blockareascoreagg(-1.), maxforeseeingwhitescore(-1.),
   maxforeseeingwhitescoreagg(-1.), setpartfwhitescore(-1.), setpartfwhitescoreagg(-1.),
   detectorchainstring( NULL ), stemsFromUnpresolved( false ), isfromunpresolved( FALSE ),
   isFinishedByFinisherUnpresolved( false ), finishedUnpresolvedBy( NULL ), seeedpool(givenseeedpool)
{

   for( int i = 0; i < nConss; ++i )
      openConss.push_back(i);

   for( int i = 0; i < nVars; ++i )
      openVars.push_back(i);

}

/** copy constructor */
Seeed::Seeed(
   const Seeed *seeedtocopy
   )
{
   scip = ( seeedtocopy->scip );
   id = seeedtocopy->id;
   nBlocks = seeedtocopy->nBlocks;
   nVars = seeedtocopy->nVars;
   nConss = seeedtocopy->nConss;
   masterConss = seeedtocopy->masterConss;
   masterVars = seeedtocopy->masterVars;
   conssForBlocks = seeedtocopy->conssForBlocks;
   varsForBlocks = seeedtocopy->varsForBlocks;
   linkingVars = seeedtocopy->linkingVars;
   stairlinkingVars = seeedtocopy->stairlinkingVars;
   openVars = seeedtocopy->openVars;
   openConss = seeedtocopy->openConss;

   isvaropen = seeedtocopy->isvaropen;
   masterconsssorted = seeedtocopy->masterconsssorted;

   isconsopen = seeedtocopy->isconsopen;

   isvarmaster = seeedtocopy->isvarmaster;
   isconsmaster = seeedtocopy->isconsmaster;

   detectorChain = seeedtocopy->detectorChain;
   detectorChainFinishingUsed = seeedtocopy->detectorChainFinishingUsed;
   detectorchaininfo = seeedtocopy->detectorchaininfo;
   hashvalue = seeedtocopy->hashvalue;
   usergiven = seeedtocopy->usergiven;
   isfromlegacymode = seeedtocopy->isfromlegacymode;
   score = seeedtocopy->score;
   borderareascore = seeedtocopy->borderareascore;
   maxwhitescore = seeedtocopy->maxwhitescore;
   bendersscore = -1.;
   benderareascore = -1.;
   changedHashvalue = seeedtocopy->changedHashvalue;
   detectorClockTimes = seeedtocopy->detectorClockTimes;
   pctVarsToBorder = seeedtocopy->pctVarsToBorder;
   pctVarsToBlock = seeedtocopy->pctVarsToBlock;
   pctVarsFromFree = seeedtocopy->pctVarsFromFree;
   pctConssToBorder = seeedtocopy->pctConssToBorder;
   pctConssToBlock = seeedtocopy->pctConssToBlock;
   pctConssFromFree = seeedtocopy->pctConssFromFree;
   usedClassifier = seeedtocopy->usedClassifier;
   classesToMaster = seeedtocopy->classesToMaster;
   classesToLinking = seeedtocopy->classesToLinking;
   isFinishedByFinisher = seeedtocopy->isFinishedByFinisher;
   ncoeffsforblockformastercons = seeedtocopy->ncoeffsforblockformastercons;
   changedHashvalue = seeedtocopy->changedHashvalue;
   nNewBlocks = seeedtocopy->nNewBlocks;
   stemsFromUnpresolved = seeedtocopy->stemsFromUnpresolved;
   finishedUnpresolvedBy = seeedtocopy->finishedUnpresolvedBy;
   isFinishedByFinisherUnpresolved = seeedtocopy->isFinishedByFinisherUnpresolved;
   isselected = false;
   detectorchainstring = NULL;
   isfromunpresolved = seeedtocopy->isfromunpresolved;
   listofancestorids = seeedtocopy->listofancestorids;

   varsforblocksorted = seeedtocopy->varsforblocksorted;
   stairlinkingvarsforblocksorted = seeedtocopy->stairlinkingvarsforblocksorted;
   conssforblocksorted = seeedtocopy->conssforblocksorted;
   linkingvarssorted = seeedtocopy->linkingvarssorted;
   mastervarssorted = seeedtocopy->mastervarssorted;

   agginfocalculated = FALSE;
   nrepblocks  = seeedtocopy->nrepblocks;
   reptoblocks = seeedtocopy->reptoblocks;
   blockstorep = seeedtocopy->blockstorep;
   pidtopidvarmaptofirst = seeedtocopy->pidtopidvarmaptofirst;
   ncoeffsforblock = seeedtocopy->ncoeffsforblock;
   calculatedncoeffsforblock = FALSE;

   blockareascore = -1.;
   maxwhitescoreagg = -1.;
   blockareascoreagg = -1.;
   maxforeseeingwhitescore = -1.;
   maxforeseeingwhitescoreagg = -1.;

   setpartfwhitescore = -1.;
   setpartfwhitescoreagg = -1.;

   isagginfoalreadytoexpensive = seeedtocopy->isagginfoalreadytoexpensive;

   seeedpool = seeedtocopy->seeedpool;

}

/** destructor */
Seeed::~Seeed()
{
   if ( detectorchainstring != NULL )
      SCIPfreeBlockMemoryArrayNull( scip, & detectorchainstring, SCIP_MAXSTRLEN );
}


/** checks whether two arrays of SCIP_Real's are identical */
static
SCIP_Bool realArraysAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            array1,             /**< first array */
   int                   array1length,       /**< length of first array */
   SCIP_Real*            array2,             /**< second array */
   int                   array2length        /**< length of second array */
   )
{
   int i;

   if( array1length != array2length )
      return FALSE;

   if( array1length == 0 )
      return TRUE;

   assert(array1 != NULL);
   assert(array2 != NULL);

   for( i = 0; i < array1length; i++ )
   {
      if( !SCIPisEQ(scip, array1[i], array2[i]) )
         return FALSE;
   }

   return TRUE;
}


/** returns true iff the second value of a is lower than the second value of b */
bool compare_blocks(
   std::pair<int, int> const & a,
   std::pair<int, int> const & b
   )
{
   return ( a.second < b.second );
}

/** adds a block, returns the number of the new block */
int Seeed::addBlock()
{
   std::vector<int> vector = std::vector<int>( 0 );

   changedHashvalue = true;

   assert( (int) conssForBlocks.size() == nBlocks );
   assert( (int) varsForBlocks.size() == nBlocks );
   assert( (int) stairlinkingVars.size() == nBlocks );

   conssForBlocks.push_back( vector );
   varsForBlocks.push_back( vector );
   stairlinkingVars.push_back( vector );
   nBlocks ++;
   return nBlocks - 1;
}

/** incorporates the the needed time of a certain detector in the detector chain */
void Seeed::addClockTime(
   SCIP_Real clocktime
   )
{
   detectorClockTimes.push_back( clocktime );
}

/** incorporates the changes from ancestor seeed */
void Seeed::addDecChangesFromAncestor(
   Seeed* ancestor
   )
{
   /** add number of new blocks */
   assert( ancestor != NULL );

   nNewBlocks.push_back( getNBlocks() - ancestor->getNBlocks() );
   pctConssFromFree.push_back( getNConss() != 0 ? ( ancestor->getNOpenconss() - getNOpenconss() ) / (SCIP_Real) getNConss() : 0. );
   pctVarsFromFree.push_back( getNVars() != 0 ? ( ancestor->getNOpenvars() - getNOpenvars() ) / (SCIP_Real) getNVars() : 0. );
   pctConssToBlock.push_back( getNConss() != 0 ?
      ( - getNOpenconss() - getNMasterconss() + ancestor->getNOpenconss() + ancestor->getNMasterconss() ) / getNConss() : 0. );
   pctVarsToBlock.push_back( getNVars() != 0 ? ( - getNOpenvars() - getNMastervars() - getNLinkingvars() - getNTotalStairlinkingvars() + ancestor->getNOpenvars()
         + ancestor->getNMastervars() + ancestor->getNLinkingvars() + ancestor->getNTotalStairlinkingvars() ) / getNVars() : 0. );
   pctConssToBorder.push_back( getNConss() != 0 ? ( getNMasterconss() - ancestor->getNMasterconss() ) / (SCIP_Real) getNConss() : 0. );
   pctVarsToBorder.push_back( getNVars() != 0 ? ( getNMastervars() + getNLinkingvars() + getNTotalStairlinkingvars() - ancestor->getNMastervars()
         - ancestor->getNLinkingvars() - ancestor->getNTotalStairlinkingvars() ) / (SCIP_Real) getNVars() : 0. );
   listofancestorids.push_back( ancestor->getID() );
}

/** adds a detector chain info */
void Seeed::addDetectorChainInfo(
   const char* decinfo
   )
{
   std::stringstream help;
   help << decinfo;
   detectorchaininfo.push_back( help.str() );
}

/** adds empty entries for all classifier statistics for a detector added to the detector chain */
void Seeed::addEmptyClassifierStatistics()
{
   std::vector<int> emptyVector( 0 );
   usedClassifier.push_back( NULL );
   classesToMaster.push_back( emptyVector );
   classesToLinking.push_back( emptyVector );
}

 /** adds number of new blocks created by a detector added to detector chain */
 void Seeed::addNNewBlocks(
    int nnewblocks
    )
 {
    nNewBlocks.push_back( nnewblocks );

    assert( nNewBlocks.size() <= detectorChain.size() );
 }

 /** adds fraction of constraints that are not longer open for a detector added to detector chain */
 void Seeed::addPctConssFromFree(
    SCIP_Real pct
    )
 {
    pctConssFromFree.push_back( pct );

    assert( pctConssFromFree.size() <= detectorChain.size() );
 }

 /** adds fraction of constraints assigned to a block for a detector added to detector chain */
 void Seeed::addPctConssToBlock(
    SCIP_Real pct
    )
 {
    pctConssToBlock.push_back( pct );

    assert( pctConssToBlock.size() <= detectorChain.size() );
 }

 /** adds fraction of constraints assigned to the border for a detector added to detector chain */
 void Seeed::addPctConssToBorder(
    SCIP_Real pct
    )
 {
    pctConssToBorder.push_back( pct );

    assert( pctConssToBorder.size() <= detectorChain.size() );
 }

 /** adds fraction of variables that are not longer open for a detector added to detector chain */
 void Seeed::addPctVarsFromFree(
    SCIP_Real pct
    )
 {
    pctVarsFromFree.push_back( pct );

    assert( pctVarsFromFree.size() <= detectorChain.size() );
 }

 /** adds fraction of variables assigned to a block for a detector added to detector chain */
 void Seeed::addPctVarsToBlock(
    SCIP_Real pct
    )
 {
    pctVarsToBlock.push_back( pct );

    assert( pctVarsToBlock.size() <= detectorChain.size() );
 }

 /** adds fraction of variables assigned to the border for a detector added to detector chain */
 void Seeed::addPctVarsToBorder(
    SCIP_Real pct
    )
 {
    pctVarsToBorder.push_back( pct );

    assert( pctVarsToBorder.size() <= detectorChain.size() );
 }

/** returns true if at least one constraint is assigned to a block */
bool Seeed::alreadyAssignedConssToBlocks()
{
   for( int b = 0; b < this->nBlocks; ++ b )
      if( conssForBlocks[b].size() != 0 )
         return true;
   return false;
}

/** assigns open conss to master according to the cons assignment information given in constoblock hashmap */
SCIP_RETCODE Seeed::assignBorderFromConstoblock(
   SCIP_HASHMAP* constoblock,
   int givenNBlocks
   )
{
   int cons;

   changedHashvalue = true;

   for( int i = 0; i < getNOpenconss(); ++ i )
   {
      cons = openConss[i];
      if( ! SCIPhashmapExists( constoblock, (void*) (size_t) cons ) )
         continue;
      if( (int) (size_t) SCIPhashmapGetImage( constoblock, (void*) (size_t) cons ) - 1 == givenNBlocks )
         bookAsMasterCons( cons );
   }

   flushBooked();

   sort();
   assert( checkConsistency( ) );
   return SCIP_OKAY;
}

/** assigns open vars to stairlinking if they can be found in two consecutive blocks, returns true if stairlinkingvars
 * are assigned */
bool Seeed::assignCurrentStairlinking(
   )
{
   std::vector<int> blocksOfOpenvar;
   bool assigned = false;
   int var;
   int cons;

   changedHashvalue = true;

   /** assign all vars included in two consecutive blocks to stairlinking */
   for( int i = 0; i < getNOpenvars(); ++ i )
   {
      blocksOfOpenvar.clear();
      var = openVars[i];
      for( int b = 0; b < nBlocks; ++ b )
      {
         for( int c = 0; c < getNConssForBlock( b ); ++ c )
         {
            cons = conssForBlocks[b][c];
            if( seeedpool->getVal( cons, var ) != 0 )
            {
               blocksOfOpenvar.push_back( b );
               break;
            }
         }
      }
      if( blocksOfOpenvar.size() == 2 && blocksOfOpenvar[0] + 1 == blocksOfOpenvar[1] )
      {
         bookAsStairlinkingVar( var, blocksOfOpenvar[0] );
         assigned = true;
      }
   }

   flushBooked();

   if( assigned )
      sort();
   return assigned;
}

/** assigns every open cons
 *  - to master if it hits blockvars of different blocks
 *  - to the respective block if it hits a blockvar of exactly one block and no stairlinking var
 *  - to master if it hits a stairlinking var but there is no block the cons may be assigned to
 *  - to the block with the lowest number of conss if it hits a stairlinking var and there are blocks the cons may be assigned to
 *  - leave it open if cannot be found in any block
 *  returns true if there is a cons that has been assigned */
bool Seeed::assignHittingOpenconss(
   )
{
   int cons;
   int var;
   int block;
   bool stairlinking; /** true if the cons includes stairlinkingvars */
   bool assigned = false; /** true if open conss get assigned in this function */
   std::vector<int>::iterator it;
   std::vector<int> blocksOfStairlinkingvars; /** first block of the stairlinkingvars which can be found in the conss */
   std::vector<int> blocksOfVars; /** blocks of the vars which can be found in the conss */
   std::vector<int> blocks; /** cons can be assigned to the blocks stored in this vector */
   std::vector<int> eraseBlock;

   changedHashvalue = true;

   for( size_t c = 0; c < openConss.size(); ++ c )
   {
      cons = openConss[c];
      stairlinking = false;

      blocksOfVars.clear();
      blocks.clear();
      blocksOfStairlinkingvars.clear();
      eraseBlock.clear();

      /** fill out blocksOfStairlinkingvars and blocksOfBlockvars */
      for( int b = 0; b < nBlocks; ++ b )
      {
         for( int v = 0; v < seeedpool->getNVarsForCons( cons ); ++ v )
         {
            var = seeedpool->getVarsForCons( cons )[v];
            if( isVarBlockvarOfBlock( var, b ) )
            {
               blocksOfVars.push_back( b );
               break;
            }
         }
      }

      for( int b = 0; b < nBlocks; ++ b )
      {
         for( int v = 0; v < seeedpool->getNVarsForCons( cons ); ++ v )
         {
            int var2 = seeedpool->getVarsForCons(cons)[v];
            std::vector<int>::iterator lb = lower_bound( stairlinkingVars[b].begin(), stairlinkingVars[b].end(), var2 );
            if( lb != stairlinkingVars[b].end() &&  *lb == var2 )
            {
               stairlinking = true;
               blocksOfStairlinkingvars.push_back( b );
               break;
            }
         }
      }

      /** fill out blocks */
      if( stairlinking && blocksOfVars.size() < 2 )
      {
         if( blocksOfVars.size() == 0 )
         {
            blocks.push_back( blocksOfStairlinkingvars[0] );
            blocks.push_back( blocksOfStairlinkingvars[0] + 1 );
            for( size_t i = 1; i < blocksOfStairlinkingvars.size(); ++ i )
            {
               for( it = blocks.begin(); it != blocks.end(); ++ it )
               {
                  if( * it != blocksOfStairlinkingvars[i] && * it != blocksOfStairlinkingvars[i] + 1 )
                     eraseBlock.push_back( * it );
               }
               for( size_t j = 0; j < eraseBlock.size(); ++ j )
               {
                  it = find( blocks.begin(), blocks.end(), eraseBlock[j] );
                  assert( it != blocks.end() );
                  blocks.erase( it );
               }
            }
         }
         else
         {
            blocks.push_back( blocksOfVars[0] );
            for( size_t i = 0; i < blocksOfStairlinkingvars.size(); ++ i )
            {
               if( blocks[0] != blocksOfStairlinkingvars[i] && blocks[0] != blocksOfStairlinkingvars[i] + 1 )
               {
                  blocks.clear();
                  break;
               }
            }
         }
      }

      if( blocksOfVars.size() > 1 )
      {
         bookAsMasterCons( cons );
         assigned = true;
      }
      else if( ! stairlinking && blocksOfVars.size() == 1 )
      {
         bookAsBlockCons( cons, blocksOfVars[0] );
         assigned = true;
      }
      else if( stairlinking && blocks.size() == 0 )
      {
         bookAsMasterCons( cons );
         assigned = true;
      }
      else if( stairlinking && blocks.size() == 1 )
      {
         bookAsBlockCons( cons, blocks[0] );
         assigned = true;
      }
      else if( stairlinking && blocks.size() > 1 )
      {
         block = blocks[0];
         for( size_t i = 1; i < blocks.size(); ++ i )
         {
            if( getNConssForBlock( i ) < getNConssForBlock( block ) )
               block = i;
         }
         bookAsBlockCons( cons, block );
         assigned = true;
      }
   }

   flushBooked();

   if( assigned )
      sort();

   return assigned;
}

/** assigns every open var
 *  - to the respective block if it hits blockconss of exactly one block
 *  - to linking if it hits blockconss of more than one different blocks
 *  returns true if there is a var that has been assigned */
bool Seeed::assignHittingOpenvars(
   )
{
   int cons;
   int var;
   std::vector<int> blocksOfOpenvar;
   bool found;
   bool assigned = false;

   changedHashvalue = true;

   /** set vars to linking, if they can be found in more than one block;
    * set vars to block if they can be found in only one block */
   for( size_t i = 0; i < openVars.size(); ++ i )
   {
      blocksOfOpenvar.clear();
      var = openVars[i];
      assert( var >= 0 && var < nVars );
      for( int b = 0; b < nBlocks; ++ b )
      {
         found = false;
         for( int c = 0; c < getNConssForBlock( b ) && ! found; ++ c )
         {
            cons = conssForBlocks[b][c];
            for( int v = 0; v < seeedpool->getNVarsForCons( cons ) && ! found; ++ v )
            {
               if( seeedpool->getVarsForCons( cons )[v] == var )
               {
                  blocksOfOpenvar.push_back( b );
                  found = true;
               }
            }
         }
      }
      if( blocksOfOpenvar.size() == 1 )
      {
         bookAsBlockVar( var, blocksOfOpenvar[0] );
         assigned = true;
      }
      else if( blocksOfOpenvar.size() > 1 )
      {
         bookAsLinkingVar( var );
         assigned = true;
      }
   }

   flushBooked();

   if( assigned )
      sort();

   return assigned;
}

/** assigns every open cons to master that hits
 *  - exactly one block var and at least one open var or
 *  - a master var */
SCIP_RETCODE Seeed::assignOpenPartialHittingConsToMaster(
   )
{
   int cons;
   int var;
   std::vector<int> blocksOfBlockvars; /** blocks with blockvars which can be found in the cons */
   std::vector<int> blocksOfOpenvar; /** blocks in which the open var can be found */
   bool master;
   bool hitsOpenVar;
   std::vector<bool> isblockhit;
   changedHashvalue = true;

   /** set openConss with more than two blockvars to master */
   for( size_t c = 0; c < openConss.size(); ++ c )
   {
      isblockhit= std::vector<bool>(getNBlocks(), false );
      blocksOfBlockvars.clear();
      master = false;
      hitsOpenVar = false;
      cons = openConss[c];


      for( int v = 0; v < seeedpool->getNVarsForCons( cons ) && ! master; ++ v )
      {
         var = seeedpool->getVarsForCons( cons )[v];

         if( isVarOpenvar( var ) )
         {
            hitsOpenVar = true;
            continue;
         }

         if( isVarMastervar( var ) )
         {
            master = true;
            bookAsMasterCons( cons );
            continue;
         }

         for( int b = 0; b < nBlocks; ++ b )
         {
            if( isblockhit[b] )
               continue;

            if( isVarBlockvarOfBlock( var, b ) )
            {
               blocksOfBlockvars.push_back( b );
               isblockhit[b] = true;
               break;
            }
         }
      }
      if( blocksOfBlockvars.size() == 1 && hitsOpenVar )
      {
         bookAsMasterCons( cons );
      }
   }

   flushBooked();

   return SCIP_OKAY;
}

/** assigns open conss/vars that hit exactly one block and at least one open var/cons to border */
SCIP_RETCODE Seeed::assignOpenPartialHittingToMaster(
   )
{
   changedHashvalue = true;
   assignOpenPartialHittingConsToMaster( );
   assignOpenPartialHittingVarsToMaster( );
   return SCIP_OKAY;
}

/** assigns every open var to linking that hits
 *  - exactly one block cons and at least one open cons */
SCIP_RETCODE Seeed::assignOpenPartialHittingVarsToMaster(
   )
{
   int cons;
   int var;
   std::vector<int> blocksOfBlockvars; /** blocks with blockvars which can be found in the cons */
   std::vector<int> blocksOfOpenvar; /** blocks in which the open var can be found */
   bool hitsOpenCons;
   std::vector<bool> isblockhit;
   bool benders;

   changedHashvalue = true;
   benders = seeedpool->isForBenders();

   /** set open var to linking if it can be found in one block and open constraint */
   for( size_t i = 0; i < openVars.size(); ++ i )
   {
      isblockhit= std::vector<bool>(getNBlocks(), false );
      blocksOfOpenvar.clear();
      var = openVars[i];
      hitsOpenCons = false;

      for( int c = 0; c < seeedpool->getNConssForVar( var ); ++ c )
      {
         cons = seeedpool->getConssForVar( var )[c];

         if( benders && isConsMastercons( cons ) )
         {
            continue;
         }

         if( isConsOpencons( cons ) )
         {
            hitsOpenCons = true;
            continue;
         }
         for( int b = 0; b < nBlocks; ++ b )
         {
            if ( isblockhit[b] )
               continue;

            if( isConsBlockconsOfBlock( cons, b ) )
            {
               blocksOfOpenvar.push_back( b );
               isblockhit[b] = true;
               break;
            }
         }

      }

//      if( benders && hitsmastercons )
//      {
//         bookAsLinkingVar( var );
//      }

      if(  blocksOfOpenvar.size() == 1 && hitsOpenCons )
      {
         bookAsLinkingVar( var );
      }
   }

   flushBooked();


   return SCIP_OKAY;
}

/** adds blocks and assigns open conss to such a new block or to master
 *  according to the cons assignment information given in constoblock hashmap */
SCIP_RETCODE Seeed::assignSeeedFromConstoblock(
   SCIP_HASHMAP* constoblock,
   int additionalNBlocks
)
{
   int oldNBlocks = nBlocks;
   int consblock;
   int cons;

   assert( additionalNBlocks >= 0 );

   changedHashvalue = true;

   for( int b = 0; b < additionalNBlocks; ++ b )
      addBlock();

   for( int i = 0; i < getNOpenconss(); ++ i )
   {
      cons = openConss[i];

      if( ! SCIPhashmapExists( constoblock, (void*) (size_t) cons ) )
         continue;
      consblock = oldNBlocks + ( (int) (size_t) SCIPhashmapGetImage( constoblock, (void*) (size_t) cons ) - 1 );
      assert( consblock >= oldNBlocks && consblock <= nBlocks );
      if( consblock == nBlocks )
         bookAsMasterCons( cons );
      else
         bookAsBlockCons( cons, consblock );
   }

   flushBooked();

   deleteEmptyBlocks(false);
   sort();
   assert( checkConsistency( ) );
   return SCIP_OKAY;
}

/** adds blocks and assigns open conss to such a new block or to master
 *  according to the cons assignment information given in constoblock vector */
SCIP_RETCODE Seeed::assignSeeedFromConstoblockVector(
   std::vector<int> constoblock,
   int additionalNBlocks
      )
{
   int oldNBlocks = nBlocks;
   int consblock;
   int cons;

   assert( additionalNBlocks >= 0 );

   changedHashvalue = true;

   for( int b = 0; b < additionalNBlocks; ++ b )
      addBlock();

   for( int i = 0; i < getNOpenconss(); ++ i )
   {
      cons = openConss[i];

      if( constoblock[cons] == - 1 )
         continue;

      consblock = oldNBlocks + ( constoblock[cons] - 1 );
      assert( consblock >= oldNBlocks && consblock <= nBlocks );
      if( consblock == nBlocks )
         bookAsMasterCons( cons );
      else
         bookAsBlockCons( cons, consblock );
   }

   flushBooked();

   deleteEmptyBlocks(false);
   sort();
   assert( checkConsistency( ) );
   return SCIP_OKAY;
}

/** books a constraint to be added to the block constraints of the given block (after calling flushBooked) */
SCIP_RETCODE Seeed::bookAsBlockCons(
   int consToBlock,
   int block
   )
{
   assert( consToBlock >= 0 && consToBlock < nConss );
   assert( isconsopen[consToBlock] );
   if( block >= nBlocks )
      setNBlocks(block+1);
   assert( block >= 0 && block < nBlocks );
   std::pair<int, int> pair( consToBlock, block );
   bookedAsBlockConss.push_back( pair );
   return SCIP_OKAY;
}

/** books a variable to be added to the block variables of the given block (after calling flushBooked) */
SCIP_RETCODE Seeed::bookAsBlockVar(
   int varToBlock,
   int block
   )
{
   assert( varToBlock >= 0 && varToBlock < nVars );
   assert( isvaropen[varToBlock] );
   assert( block >= 0 && block < nBlocks );
   std::pair<int, int> pair( varToBlock, block );
   bookedAsBlockVars.push_back( pair );
   return SCIP_OKAY;
}

/** books a constraint to be added to the master constraints (after calling flushBooked)*/
SCIP_RETCODE Seeed::bookAsMasterCons(
   int consToMaster
   )
{
   assert( consToMaster >= 0 && consToMaster < nConss );
   assert(isconsopen[consToMaster]);
   bookedAsMasterConss.push_back( consToMaster );
   return SCIP_OKAY;
}

/** books a variable to be added to the master variables (after calling flushBooked) */
SCIP_RETCODE Seeed::bookAsMasterVar(
   int varToMaster
   )
{
   assert( varToMaster >= 0 && varToMaster < nVars );
   assert( isvaropen[varToMaster]);
   bookedAsMasterVars.push_back( varToMaster );
   return SCIP_OKAY;
}

/** books a variable to be added to the linking variables (after calling flushBooked) */
SCIP_RETCODE Seeed::bookAsLinkingVar(
   int varToLinking
   )
{
   assert( varToLinking >= 0 && varToLinking < nVars );
   assert( isvaropen[varToLinking]);
   bookedAsLinkingVars.push_back( varToLinking );
   return SCIP_OKAY;
}

/** books a variable to be added to the stairlinking variables of the given block and the following block (after calling flushBooked) */
SCIP_RETCODE Seeed::bookAsStairlinkingVar(
   int varToStairlinking,
   int firstBlock
   )
{
   assert( isvaropen[varToStairlinking]);
   assert( varToStairlinking >= 0 && varToStairlinking < nVars );
   assert( firstBlock >= 0 && firstBlock < ( nBlocks - 1 ) );
   std::pair<int, int> pair( varToStairlinking, firstBlock );
   bookedAsStairlinkingVars.push_back( pair );
   return SCIP_OKAY;
}

/**
 * prepare the seeed such that all predecessors have the folloing property:
 * all variables in the master problem are binary variables
 * thus all other variables are assigned to a block
 *
 * reqirement: all constraints and variables are open when this method is called
 */
void Seeed::initOnlyBinMaster(){

   std::vector<int> constoblocks = std::vector<int>(nConss, -1);
   std::vector<int> vartoblocks = std::vector<int>(nVars, -1);
   std::vector<int> mergeblocks = std::vector< int>(0);
   int potentialnewblocks = 0;


   for( int openvar = 0; openvar < getNVars(); ++openvar )
   {
      /** check if var is NOT a binary variable */
      if( SCIPvarGetType(seeedpool->getScipVar(openvar) ) != SCIP_VARTYPE_BINARY )
      {
         int blockid = -1;
         std::vector<bool> hitblock = std::vector<bool>(potentialnewblocks, false);
         std::vector<int> hittenblocks = std::vector<int>(0);
         for( int c = 0; c < seeedpool->getNConssForVar(openvar); ++c )
         {
            int cons = seeedpool->getConssForVar(openvar)[c];
            if( constoblocks[cons] != -1 )
            {
               if( !hitblock[constoblocks[cons]] )
               {
                  hitblock[constoblocks[cons]] = true;
                  hittenblocks.push_back(constoblocks[cons]);
               }
            }
         }

         if( hittenblocks.size() == 0 )
         {
            ++potentialnewblocks;
            blockid = potentialnewblocks - 1;
            mergeblocks.push_back(blockid);
            vartoblocks[openvar] = blockid;
         }

         if( hittenblocks.size() == 1 )
         {
            blockid = hittenblocks[0];
            vartoblocks[openvar] = hittenblocks[0];
         }

         if( hittenblocks.size() > 1 )
         {
            int minblock = hittenblocks[0];
            for( size_t i = 1; i < hittenblocks.size(); ++i )
            {
               if( hittenblocks[i] < minblock )
                  minblock = hittenblocks[i];
            }
            int mergeto = mergeblocks[minblock];
            for( size_t i = 0; i < hittenblocks.size(); ++i )
            {
               mergeblocks[hittenblocks[i]] = mergeto;
            }
            blockid = minblock;
            vartoblocks[openvar] = blockid;
         }

         assert(blockid != -1);
         for( int c = 0; c < seeedpool->getNConssForVar(openvar); ++c )
         {
            int cons = seeedpool->getConssForVar(openvar)[c];
            constoblocks[cons] = blockid;
         }
      }
   }


   int nnewblocks = 0;
   std::vector<int> currblocks = std::vector<int>(mergeblocks.size(), -1);
   for( int b = 0; b < potentialnewblocks; ++b )
      if( mergeblocks[b] == b )
      {
         currblocks[b] = nnewblocks;
         ++nnewblocks;
      }
      else
         currblocks[b] = currblocks[mergeblocks[b]];

   setNBlocks(nnewblocks);

   for( int openvar = 0; openvar < getNVars(); ++openvar )
   {
      if( vartoblocks[openvar] != -1)
         bookAsBlockVar(openvar, currblocks[vartoblocks[openvar]]  );
   }

   for( int cons = 0; cons < getNConss(); ++cons )
   {
      if( constoblocks[cons] != -1)
         bookAsBlockCons(cons, currblocks[constoblocks[cons]]);
   }

   flushBooked();

   return;

}

SCIP_Bool Seeed::isAgginfoToExpensive()
{

   int limitfornconss;
   int limitfornvars;

   if( isagginfoalreadytoexpensive )
      return TRUE;

   SCIPgetIntParam(seeedpool->getScip(), "detection/aggregation/limitnconssperblock", &limitfornconss);
   SCIPgetIntParam(seeedpool->getScip(), "detection/aggregation/limitnvarsperblock", &limitfornvars);



   /** check if calculating aggregation information is too expensive */
   for( int b1 = 0; b1 < getNBlocks() ; ++b1 )
   {
      for( int b2 = b1+1; b2 < getNBlocks(); ++b2 )
      {
         if( getNVarsForBlock(b1) != getNVarsForBlock(b2) )
            continue;

         if( getNConssForBlock(b1) != getNConssForBlock(b2) )
            continue;

         SCIPdebugMessage("Checking  if agg info is too expensive for blocks %d and %d, nconss: %d, nvars: %d . \n", b1, b2, getNConssForBlock(b2), getNVarsForBlock(b2) );
         if( getNConssForBlock(b2) >= limitfornconss || getNVarsForBlock(b2) >= limitfornvars )
         {
            SCIPdebugMessage("Calculating agg info is too expensive, nconss: %d, nvars: %d . \n", getNConssForBlock(b2), getNVarsForBlock(b2) );
            isagginfoalreadytoexpensive = true;
            return TRUE;
         }
      }

   }

   /** check if there are too many master coeffs */

   SCIPdebugMessage("Calculated: agg info is NOT too expensive.\n");
   return FALSE;
}



/** checks if aggregation of sub problems is possible and stores the corresponding aggregation information; */
  void Seeed::calcAggregationInformation(
     )
  {
#ifdef WITH_BLISS
     SCIP_Bool tooexpensive;
#endif
     SCIP_Bool aggisnotactive;
     SCIP_Bool discretization;
     SCIP_Bool aggregation;

     int nreps = 1;

     if( agginfocalculated )
        return;

     if( !isComplete() )
        return;

#ifdef WITH_BLISS
     if( isAgginfoToExpensive() )
        tooexpensive = TRUE;
     else
        tooexpensive = FALSE;
#endif

     SCIPgetBoolParam(seeedpool->getScip(), "relaxing/gcg/aggregation", &aggregation);
     SCIPgetBoolParam(seeedpool->getScip(), "relaxing/gcg/discretization", &discretization);

     if( discretization && aggregation )
        aggisnotactive = FALSE;
     else
        aggisnotactive = TRUE;

     std::vector<std::vector<int>> identblocksforblock( getNBlocks(), std::vector<int>(0) );

     blockstorep = std::vector<int>(getNBlocks(), -1);

     for( int b1 = 0; b1 < getNBlocks() ; ++b1 )
     {
        std::vector<int> currrep = std::vector<int>(0);
        std::vector< std::vector<int> > currrepvarmapforthisrep =std::vector<std::vector<int>>(0);
        std::vector<int> identityvec = std::vector<int>(0);


        if( !identblocksforblock[b1].empty() )
           continue;

        for( int i = 0; i  < getNVarsForBlock(b1); ++i )
           identityvec.push_back(i);

        currrep.push_back(b1);
        currrepvarmapforthisrep.push_back(identityvec);


        for( int b2 = b1+1; b2 < getNBlocks(); ++b2 )
        {
           SCIP_Bool identical;
           SCIP_Bool notidentical;
           std::vector<int> varmap;
           SCIP_HASHMAP* varmap2;

           notidentical = FALSE;
           identical = FALSE;

           if( !identblocksforblock[b2].empty() )
              continue;

           if( aggisnotactive )
              continue;


           SCIP_CALL_ABORT( SCIPhashmapCreate(  &varmap2,
                          SCIPblkmem(seeedpool->getScip()),
                          5 * getNVarsForBlock(b1)+1) ); /* +1 to deal with empty subproblems */

           SCIPdebugMessage("Check identity for block %d and block %d!\n", b1, b2);

           checkIdenticalBlocksTrivial( b1, b2, &notidentical);

           if( !notidentical )
           {
              checkIdenticalBlocksBrute( b1, b2, varmap, varmap2, &identical);

#ifdef WITH_BLISS
              if( !tooexpensive && !identical )
                 checkIdenticalBlocksBliss(b1, b2, varmap, varmap2, &identical);
#endif
           }
           else
              identical = FALSE;

           if( identical )
           {
              SCIPdebugMessage("Block %d is identical to block %d!\n", b1, b2);
              identblocksforblock[b1].push_back(b2);
              identblocksforblock[b2].push_back(b1);
              currrep.push_back(b2);
              /** handle varmap */
              currrepvarmapforthisrep.push_back(varmap);

           }
           else
           {
              SCIPdebugMessage("Block %d is not identical to block %d!\n", b1, b2);
           }
           SCIPhashmapFree(&varmap2);
        }

        reptoblocks.push_back( currrep );
        pidtopidvarmaptofirst.push_back(currrepvarmapforthisrep);
        for( size_t i = 0; i < currrep.size(); ++i )
           blockstorep[currrep[i]] = nreps-1;
        ++nreps;

     }
     nrepblocks = nreps-1;

     agginfocalculated = TRUE;

     return;
  }



/** calculates the hashvalue of the seeed for comparing */
void Seeed::calcHashvalue()
{
   std::vector<std::pair<int, int>> blockorder = std::vector < std::pair<int, int> > ( 0 );
   long hashval = 0;
   long borderval = 0;

   /** find sorting for blocks (non decreasing according smallest row index) */
   for( int i = 0; i < this->nBlocks; ++ i )
   {
      if( this->conssForBlocks[i].size() > 0 )
         blockorder.push_back( std::pair<int, int>( i, this->conssForBlocks[i][0] ) );
      else
      {
         assert( this->varsForBlocks[i].size() > 0 );
         blockorder.push_back( std::pair<int, int>( i, this->getNConss() + this->varsForBlocks[i][0] ) );
      }
   }

   std::sort( blockorder.begin(), blockorder.end(), compare_blocks );

   for( int i = 0; i < nBlocks; ++ i )
   {
      long blockval = 0;
      int blockid = blockorder[i].first;

      for( size_t tau = 0; tau < conssForBlocks[blockid].size(); ++ tau )
      {
         blockval += ( 2 * conssForBlocks[blockid][tau] + 1 ) * pow( 2, tau % 16 );
      }

      hashval += primes[i % ( nPrimes - 1 )] * blockval;
   }

   for( size_t tau = 0; tau < masterConss.size(); ++ tau )
   {
      borderval += ( 2 * masterConss[tau] + 1 ) * pow( 2, tau % 16 );
   }

   hashval += primes[nBlocks % nPrimes] * borderval;

   hashval += primes[( nBlocks + 1 ) % nPrimes] * openVars.size();

   this->hashvalue = hashval;
}

/** calculates the number of nonzero coefficients for the blocks */
SCIP_RETCODE Seeed::calcNCoeffsForBlocks(
){

   if( calculatedncoeffsforblock )
      return SCIP_OKAY;

   ncoeffsforblock = std::vector<int>(getNBlocks(), 0);
   int counter;



   for( int b = 0; b < getNBlocks(); ++b )
   {
      counter = 0;
      for( int blco = 0; blco < getNConssForBlock(b); ++blco )
      {
            int consid = getConssForBlock(b)[blco];

            for( int cva = 0; cva < seeedpool->getNVarsForCons(consid) ;++cva )
               if( isVarBlockvarOfBlock(seeedpool->getVarsForCons(consid)[cva], b ) )
                  ++counter;
      }
      ncoeffsforblock[b] = counter;
   }

   counter = 0;

   for( int mco = 0; mco < getNMasterconss(); ++mco )
   {
         int consid = getMasterconss()[mco];

         counter += seeedpool->getNVarsForCons(consid);
   }
   ncoeffsformaster = counter;

   calculatedncoeffsforblock = TRUE;

   return SCIP_OKAY;
}



/** reassigns linking vars stairlinkingvars if possible
 *  potentially reorders blocks for making a maximum number of linking vars stairlinking
 *  if all vars that connect exactly two blocks have a staircase structure, all of them become stairlinkingvars
 *  otherwise, the stairlinking assignment is done greedily
 *  precondition: seeed does not have any stairlinking vars */
void Seeed::calcStairlinkingVars(
   )
{
   assert( getNTotalStairlinkingvars() == 0 );

   /* data structure containing pairs of varindices and blocknumbers */
   std::vector< std::pair< int, std::vector< int > > > blocksOfVars = findLinkingVarsPotentiallyStairlinking( );

   /* if there are no vars that are potentially stairlinking, return without further calculations */
   if( blocksOfVars.size() == 0 )
      return;

   GraphGCG* g = new GraphGCG( getNBlocks(), true );

   /* create block graph: every block is represented by a node and two nodes are adjacent if there exists a
    * var that potentially links these blocks, the edge weight is the number of such variables */
   for( int i = 0; i < (int) blocksOfVars.size(); ++i )
   {
      assert( blocksOfVars[i].second.size() == 2 );
      int v = blocksOfVars[i].second[0];
      int w = blocksOfVars[i].second[1];

      if ( g->isEdge( v, w ) )
      {
         g->setEdge( v, w, g->getEdgeWeight( v, w ) + 1 );
      }
      else
      {
         g->setEdge( v, w, 1 );
      }
   }


   bool isstaircase = true; /* maintains information whether staircase structure is still possible */
   std::vector< int > sources( 0 ); /* all nodes with degree one */
   std::vector< bool > marks( getNBlocks() ); /* a node is marked if its degree is zero or it is reachable from a source  */

   /* firstly, check whether every node has an degree of at most 2 */
   for( int b = 0; b < getNBlocks(); ++b )
   {
      if( g->getNNeighbors( b ) > 2 )
      {
         isstaircase = false;
         break;
      }
      else if( g->getNNeighbors( b ) == 1 )
      {
         sources.push_back( b );
      }
      else if ( g->getNNeighbors( b ) == 0 )
      {
         marks[b] = true;
      }
   }

   /* secondly, check whether there exists a circle in the graph by moving along all paths starting from a source */
   for( int s = 0; s < (int) sources.size() && isstaircase; ++s )
   {
      int curBlock = sources[s];
      if( marks[curBlock] )
         continue;

      marks[curBlock] = true;

      /* check whether there is an unmarked neighbor
       * if there is none, a circle is detected */
      do
      {
         std::vector< int > neighbors = g->getNeighbors( curBlock );
         if( !marks[neighbors[0]] )
         {
            marks[neighbors[0]] = true;
            curBlock = neighbors[0];
         }
         else if ( !marks[neighbors[1]] )
         {
            marks[neighbors[1]] = true;
            curBlock = neighbors[1];
         }
         else
         {
            isstaircase = false;
            break;
         }
      }
      while( g->getNNeighbors( curBlock ) != 1 );
   }

   /* thirdly, check whether all nodes with neighbors are reachable from a source,
    * since there is a circle if this is not the case */
   for( int b = 0; b < getNBlocks() && isstaircase; ++b )
   {
      if( !marks[b] )
      {
         isstaircase = false;
         break;
      }
   }

   if( isstaircase )
   {
      changeBlockOrderStaircase( g );
   }
   else
   {
      /* check if stairlinkingheuristic is activated */
      SCIP_Bool stairlinkingheur;
      SCIPgetBoolParam(scip, "detection/legacymode/stairlinkingheur", &stairlinkingheur);
      if( !stairlinkingheur )
         return;

      changeBlockOrderGreedily( g );
   }

   findVarsLinkingToStairlinking( );

   assert( checkConsistency( ) );
}

/** changes the block order in a way such that all linking vars that are potentially stairlinking
 *  may be reassigned to stairlinking
 *  precondition: all potentially stairlinking vars have a staircase structure */
void Seeed::changeBlockOrderStaircase(
   GraphGCG* g
   )
{
   int blockcounter = 0; /* counts current new block to assign an old one to */
   std::vector< int > blockmapping( getNBlocks() ); /* stores new block order */
   for( int b = 0; b < getNBlocks(); ++b )
      blockmapping[b] = -1;

   for( int b = 0; b < getNBlocks(); ++b )
   {
      if( g->getNNeighbors( b ) == 0 )
      {
         /* if block does not have a neighbor, just assign it to current blockindex */
         assert( blockmapping[b] == -1 );
         blockmapping[b] = blockcounter;
         ++blockcounter;
      }
      else if( blockmapping[b] == -1 && g->getNNeighbors( b ) == 1 )
      {
         /* if the block is the source of an yet unconsidered path, assign whole path to ascending new block ids */
         int curBlock = b;
         blockmapping[b] = blockcounter;
         ++blockcounter;

         do
         {
            std::vector< int > neighbors = g->getNeighbors( curBlock );

            if( blockmapping[neighbors[0]] == -1 )
            {
               blockmapping[neighbors[0]] = blockcounter;
               curBlock = neighbors[0];
            }
            else if ( blockmapping[neighbors[1]] == -1 )
            {
               blockmapping[neighbors[1]] = blockcounter;
               curBlock = neighbors[1];
            }
            else
            {
               assert( false );
            }
            ++blockcounter;
         }
         while( g->getNNeighbors( curBlock ) != 1 );
      }
   }

   changeBlockOrder( blockmapping );
}

/** changes the block order in a way such that some linking vars that are potentially stairlinking
 *  may be reassigned to stairlinking using a greedy method */
void Seeed::changeBlockOrderGreedily(
   GraphGCG* g
   )
{
   int blockcounter = 0; /* counts current new block to assign an old one to */
   std::vector< int > blockmapping( getNBlocks() ); /* stores new block order */
   for( int b = 0; b < getNBlocks(); ++b )
      blockmapping[b] = -1;

   for( int b = 0; b < getNBlocks(); ++b )
   {
      if( g->getNNeighbors( b ) == 0 )
      {
         /* if block does not have a neighbor, just assign it to current blockindex */
         assert( blockmapping[b] == -1 );
         blockmapping[b] = blockcounter;
         ++blockcounter;
      }
      else if( blockmapping[b] == -1 )
      {
         /* if the block is part of an yet unconsidered path, walk along this path greedily
          * and assign whole path to ascending new block ids */
         int curBlock = b;
         blockmapping[b] = blockcounter;
         int maxNeighbor;
         int maxNeighborVal;

         do
         {
            ++blockcounter;
            std::vector< int > neighbors = g->getNeighbors( curBlock );

            maxNeighbor = -1;
            maxNeighborVal = -1;

            /* find yet unassigned neighbor block with maximum number of stairlinking vars connecting it to current block */
            for( int i = 0; i < (int) neighbors.size(); ++i )
            {
               if( blockmapping[neighbors[i]] == -1 && g->getEdgeWeight( curBlock, neighbors[i] ) > maxNeighborVal )
               {
                  maxNeighbor = neighbors[i];
                  maxNeighborVal = g->getEdgeWeight( curBlock, neighbors[i] );
               }
            }

            if( maxNeighbor != -1 )
            {
               assert( blockmapping[maxNeighbor] == -1 );
               blockmapping[maxNeighbor] = blockcounter;
               curBlock = maxNeighbor;
            }
         }
         while( maxNeighbor != -1 );
      }
   }

   changeBlockOrder( blockmapping );
}

/** changes the order of the blocks according to the given mapping
 *  precondition: given mapping needs to be an adequately sized permutation */
void Seeed::changeBlockOrder(
   std::vector<int> oldToNewBlockIndex
   )
{
   assert((int ) oldToNewBlockIndex.size() == getNBlocks() );
   assert( getNTotalStairlinkingvars() == 0 );

   std::vector< std::vector< int > > newconssforblocks( getNBlocks() );
   std::vector< std::vector< int > > newvarsforblocks( getNBlocks() );

   for( int b = 0; b < getNBlocks(); ++b )
   {
      assert( 0 <= oldToNewBlockIndex[b] && oldToNewBlockIndex[b] < getNBlocks() );

      newconssforblocks[oldToNewBlockIndex[b]] = conssForBlocks[b];
      newvarsforblocks[oldToNewBlockIndex[b]] = varsForBlocks[b];
   }

   conssForBlocks = newconssforblocks;
   varsForBlocks = newvarsforblocks;
}

/** returns whether all cons are assigned and deletes the vector open cons if all are assigned */
bool Seeed::checkAllConssAssigned()
{
   for( size_t i = 0; i < openConss.size(); ++ i )
   {
      bool consfound = false;
      for( size_t k = 0; k < masterConss.size(); ++ k )
      {
         if( openConss[i] == masterConss[k] )
         {
            consfound = true;
            break;
         }
      }
      for( int b = 0; b < nBlocks && ! consfound; ++ b )
      {
         for( size_t k = 0; k < conssForBlocks[b].size(); ++ k )
         {
            if( openConss[i] == conssForBlocks[b][k] )
            {
               consfound = true;
               break;
            }
         }
      }
      if( ! consfound )
      {
         return false;
      }
   }
   openConss.clear();
   isconsopen = std::vector<bool>(nConss, false);
   return true;
}




/** returns true if the assignments in the seeed are consistent */
bool Seeed::checkConsistency(
   )
{
   std::vector<bool> openVarsBool( nVars, true );
   std::vector<int> stairlinkingvarsvec( 0 );
   std::vector<int>::const_iterator varIter = linkingVars.begin();
   std::vector<int>::const_iterator varIterEnd = linkingVars.end();

   int value;

   /** check if nblocks is set appropriately */
   if( nBlocks != (int) conssForBlocks.size() )
   {
      SCIPwarningMessage(scip, "In (seeed %d) nBlocks %d and size of conssForBlocks %d are not identical! \n" , id, nBlocks, conssForBlocks.size() );
      assert( false );
      return false;
   }

   if( nBlocks != (int) varsForBlocks.size() )
   {
      SCIPwarningMessage(scip, "In (seeed %d) nBlocks %d and size of varsForBlocks %d are not identical! \n" , id, nBlocks, varsForBlocks.size() );
      assert( false );
      return false;
   }

   /** check for empty (row- and col-wise) blocks */

   for( int b = 0; b < nBlocks; ++ b )
   {
      if( conssForBlocks[b].size() == 0 && varsForBlocks[b].size() == 0 )
      {
         SCIPwarningMessage(scip, "In (seeed %d) block %d is empty! \n" , id, b );
         this->displaySeeed();
         assert( false );
         return false;
      }
   }

   /**check variables (every variable is assigned at most once) */
   for( ; varIter != varIterEnd; ++ varIter )
   {
      if( ! openVarsBool[ * varIter] )
      {
         SCIPwarningMessage(scip, "In (seeed %d) linking variable with index %d is already assigned! \n" , id, * varIter );

         assert( false );
         return false;
      }
      openVarsBool[ * varIter] = false;
   }

   for( int b = 0; b < nBlocks; ++ b )
   {
      varIterEnd = varsForBlocks[b].end();
      for( varIter = varsForBlocks[b].begin(); varIter != varIterEnd; ++ varIter )
      {
         if( ! openVarsBool[ * varIter] )
         {
            SCIPwarningMessage(scip, "In (seeed %d) variable with index %d is already assigned but also assigned to block %d! \n" , id, * varIter, b );
            assert( false );
            return false;
         }
         openVarsBool[ * varIter] = false;
      }
   }

   varIterEnd = masterVars.end();
   for( varIter = masterVars.begin(); varIter != varIterEnd; ++ varIter )
   {
      if( ! openVarsBool[ * varIter] )
      {
         SCIPwarningMessage(scip, "In (seeed %d) variable with index %d is already assigned but also assigned to master! \n" , id, * varIter);
         assert( false );
         return false;
      }
      openVarsBool[ * varIter] = false;
   }

   for( int b = 0; b < nBlocks; ++ b )
   {
      varIter = stairlinkingVars[b].begin();
      varIterEnd = stairlinkingVars[b].end();
      for( ; varIter != varIterEnd; ++ varIter )
      {
         if( ! openVarsBool[ * varIter] )
         {
            SCIPwarningMessage(scip, "In (seeed %d) variable with index %d is already assigned but also assigned to stairlinking block %d! \n" , id, * varIter, b );
            assert( false );
            return false;
         }
         openVarsBool[ * varIter] = false;
      }
      if( ( b == nBlocks - 1 ) && ( (int) stairlinkingVars[b].size() != 0 ) )
      {
         SCIPwarningMessage(scip, "In (seeed %d) variable with index %d is is as assigned as stairlinking var of last block! \n" , id, * varIter );
         assert( false );
         return false;
      }
   }

   /** check if all not assigned variables are open vars */
   for( int v = 0; v < nVars; ++ v )
   {
      if( openVarsBool[v] == true && isVarOpenvar( v ) == false )
      {
         SCIPwarningMessage(scip, "In (seeed %d) variable with index %d is not assigned and not an open var! \n" , id, v );
         assert( false );
         return false;
      }
   }

   /** check if all open vars are not assigned */
   for( size_t i = 0; i < openVars.size(); ++ i )
   {
      if( openVarsBool[openVars[i]] == false )
      {
         SCIPwarningMessage(scip, "In (seeed %d) variable with index %d is an open var but assigned! \n" , id, openVars[i]  );
         assert( false );
         return false;
      }
   }

   for( size_t i = 0; i < openVarsBool.size(); ++ i )
   {
      if( openVarsBool[i] != isvaropen[i] )
      {
         SCIPwarningMessage(scip, "In (seeed %d) variable with index %d is causes asynchronity with isvaropen array ! \n" , id, openVars[i]  );
         assert( false );
         return false;

      }
   }

   /** check constraints (every constraint is assigned at most once) */
   std::vector<bool> openConssBool( nConss, true );
   std::vector<int> openConssVec( 0 );
   std::vector<int>::const_iterator consIter = masterConss.begin();
   std::vector<int>::const_iterator consIterEnd = masterConss.end();

   for( ; consIter != consIterEnd; ++ consIter )
   {
      if( ! openConssBool[ * consIter] )
      {
         SCIPwarningMessage(scip, "In (seeed %d) constraint with index %d is at least two times assigned as a master constraint! \n" , id, * consIter  );
         assert( false );
         return false;
      }
      openConssBool[ * consIter] = false;
   }

   for( int b = 0; b < nBlocks; ++ b )
   {
      consIterEnd = conssForBlocks[b].end();
      for( consIter = conssForBlocks[b].begin(); consIter != consIterEnd; ++ consIter )
      {
         if( ! openConssBool[ * consIter] )
         {
            SCIPwarningMessage(scip, "In (seeed %d) constraint with index %d is already assigned but also assigned to block %d! \n" , id, * consIter, b  );
            assert( false );
            return false;
         }
         openConssBool[ * consIter] = false;
      }
   }

   /** check if all not assigned constraints are open cons */
   for( int v = 0; v < nConss; ++ v )
   {
      if( openConssBool[v] == true && isConsOpencons( v ) == false )
      {
         SCIPwarningMessage(scip, "In (seeed %d) constraint with index %d is not assigned and not an open cons! \n" , id, v  );
         assert( false );
         return false;
      }
   }

   /** check if all open conss are not assigned */
   for( size_t i = 0; i < openConss.size(); ++ i )
   {
      if( openConssBool[openConss[i]] == false )
      {
         SCIPwarningMessage(scip, "In (seeed %d) constraint with index %d is an open cons but assigned! \n" , id,  openConss[i] );
         assert( false );
         return false;
      }
   }

   /** check if the seeed is sorted */
   for( int b = 0; b < nBlocks; ++ b )
   {
      value = - 1;
      for( int v = 0; v < getNVarsForBlock( b ); ++ v )
      {
         if( ! ( value < getVarsForBlock( b )[v] ) )
         {
            SCIPwarningMessage(scip, "In (seeed %d) variables of block %d are not sorted! \n" , id,  b );
            assert( false );
            return false;
         }
         value = getVarsForBlock( b )[v];
      }
   }
   for( int b = 0; b < nBlocks; ++ b )
   {
      value = - 1;
      for( int v = 0; v < getNStairlinkingvars( b ); ++ v )
      {
         if( ! ( value < getStairlinkingvars( b )[v] ) )
         {
            SCIPwarningMessage(scip, "In (seeed %d) stairlinking variables of block %d are not sorted! \n" , id,  b );
            assert( false );
            return false;
         }
         value = getStairlinkingvars( b )[v];
      }
   }
   value = - 1;
   for( int v = 0; v < getNLinkingvars(); ++ v )
   {
      if( ! ( value < getLinkingvars()[v] ) )
      {
         SCIPwarningMessage(scip, "In (seeed %d) linking variables are not sorted! \n" , id );
         assert( false );
         return false;
      }
      value = getLinkingvars()[v];
   }
   value = - 1;
   for( int v = 0; v < getNMastervars(); ++ v )
   {
      if( ! ( value < getMastervars()[v] ) )
      {
         SCIPwarningMessage(scip, "In (seeed %d) master variables are not sorted! \n" , id );
         assert( false );
         return false;
      }
      value = getMastervars()[v];
   }
   for( int b = 0; b < nBlocks; ++ b )
   {
      value = - 1;
      for( int v = 0; v < getNConssForBlock( b ); ++ v )
      {
         if( ! ( value < getConssForBlock( b )[v] ) )
         {
            SCIPwarningMessage(scip, "In (seeed %d) constraints of block %d are not sorted! \n" , id,  b );
            assert( false );
            return false;
         }
         value = getConssForBlock( b )[v];
      }
   }
   value = - 1;
   for( int v = 0; v < getNMasterconss(); ++ v )
   {
      if( ! ( value < getMasterconss()[v] ) )
      {
         SCIPwarningMessage(scip, "In (seeed %d) master constraints are not sorted! \n" , id);
         assert( false );
         return false;
      }
      value = getMasterconss()[v];
   }

   /** check if variables hitting a cons are either in the cons's block or border or still open */
   for( int b = 0; b < nBlocks; ++ b )
   {
      for( int c = 0; c < getNConssForBlock( b ); ++ c )
      {
         for( int v = 0; v < seeedpool->getNVarsForCons( getConssForBlock( b )[c] ); ++ v )
         {
            int varid = seeedpool->getVarsForCons( getConssForBlock( b )[c] )[v];

            if( ! ( isVarBlockvarOfBlock( varid, b ) || isVarLinkingvar( varid ) || isVarStairlinkingvarOfBlock( varid, b )
               || isVarOpenvar( varid ) ) )
            {
               SCIP_Bool partofblock;

               partofblock = FALSE;

               SCIPwarningMessage( scip,
                  "This should only happen during translation of (partial) decompositions from orginal to transformed problem, and means that translation has failed for this particaluar partial decomposition. Variable %d is not part of block %d or linking or open as constraint %d suggests! \n ", varid, b,
                  getConssForBlock( b )[c] );

               for( int b2 = 0; b2 < getNBlocks(); ++b2 )
               {
                  if ( isVarBlockvarOfBlock(varid, b2 ) )
                  {
                     partofblock = TRUE;
                     SCIPwarningMessage( scip,
                        "instead Variable %d is part of block %d  \n ", varid, b2 );
                     break;
                  }
               }

               if( !partofblock )
               {
                  if( isvarmaster[varid] )
                     SCIPwarningMessage( scip,
                                             "instead Variable %d is part of master  \n ", varid );
                  else
                     SCIPwarningMessage( scip,
                                                               "in fact Variable %d is completely unassigned  \n ", varid );
               }
               return false;
            }
         }
      }
   }

   return true;
}


#ifdef WITH_BLISS
/** checks blocks for identity by graph automorphism check done by bliss, identity is only found if variables are in correct order */
void Seeed::checkIdenticalBlocksBliss(
   int                  b1,
   int                  b2,
   std::vector<int>&    varmap,         /**< maps variable indices (corresponding to  seeedpool indices) of prob2 to prob1 */
   SCIP_HASHMAP*        varmap2,
   SCIP_Bool*           identical
   )
{
   *identical = FALSE;
   SCIP_HASHMAP* consmap;
   SCIP_Result result;

   varmap = std::vector<int>(getNVarsForBlock(b1), -1);

   SCIP_CALL_ABORT( SCIPhashmapCreate(&consmap,
      SCIPblkmem(seeedpool->getScip() ),
      getNConssForBlock(b1)+1) ); /* +1 to deal with empty subproblems */


   SCIPdebugMessage("obvious test fails, start building graph \n");

   cmpGraphPairNewdetection(seeedpool->getScip(), (SEEED_WRAPPER*) this, b1, b2, &result, varmap2, consmap );
   if ( result == SCIP_SUCCESS )
   {
      *identical = TRUE;
      /** TODO translate varmaps */
      for( int var2idinblock = 0; var2idinblock < getNVarsForBlock(b2) ; ++var2idinblock )
      {
         SCIP_VAR* var2;
         SCIP_VAR* var1;
         int var1idinblock;
         int var1id;

         var2 = seeedpool->getVarForIndex(getVarsForBlock(b2)[var2idinblock]);
         var1 = (SCIP_VAR*) SCIPhashmapGetImage(varmap2, (void*) var2);
         var1id = seeedpool->getIndexForVar(var1);
         var1idinblock = getVarProbindexForBlock(var1id, b1);
         varmap[var2idinblock] = var1idinblock;
      }

   }
   else
      *identical = FALSE;

   SCIPhashmapFree(&consmap);

   return;

}
#endif


/** checks blocks for identity by brute force, identity is only found if variables are in correct order */
void Seeed::checkIdenticalBlocksBrute(
   int                  b1,
   int                  b2,
   std::vector<int>&    varmap,         /**< maps variable indices (corresponding to  seeedpool indices) of prob2 to prob1 */
   SCIP_HASHMAP*        varmap2,
   SCIP_Bool*           identical
   )
{


   *identical = FALSE;
   SCIPdebugMessage("check block %d and block %d for identity...\n", b1, b2);
   varmap = std::vector<int>(getNVars(), -1);


   /** check variables */
   for( int i = 0; i < getNVarsForBlock(b1); ++i )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      var1 = seeedpool->getVarForIndex( getVarsForBlock(b1)[i] );
      var2 = seeedpool->getVarForIndex( getVarsForBlock(b2)[i] );


      if( !SCIPisEQ(scip, SCIPvarGetObj(var1), SCIPvarGetObj(var2) ) )
      {
         SCIPdebugMessage("--> obj differs for var %s and var %s!\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
             return;
      }
      if( !SCIPisEQ(scip, SCIPvarGetLbGlobal(var1), SCIPvarGetLbGlobal(var2) ) )
      {
         SCIPdebugMessage("--> lb differs for var %s and var %s!\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
             return;
      }
      if( !SCIPisEQ(scip, SCIPvarGetUbGlobal(var1), SCIPvarGetUbGlobal(var2) ) )
      {
         SCIPdebugMessage("--> ub differs for var %s and var %s!\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
             return;
      }
      if( SCIPvarGetType(var1) != SCIPvarGetType(var2) )
      {
         SCIPdebugMessage("--> type differs for var %s and var %s!\n", SCIPvarGetName(var1), SCIPvarGetName(var2));
             return;
      }

      for( int mc = 0; mc < getNMasterconss(); ++mc )
      {

         if( !SCIPisEQ(scip, seeedpool->getVal(getMasterconss()[mc], getVarsForBlock(b1)[i]), seeedpool->getVal(getMasterconss()[mc], getVarsForBlock(b2)[i])  ))
         {
            SCIPdebugMessage("--> master coefficients differ for var %s (%f) and var %s  (%f) !\n", SCIPvarGetName(  seeedpool->getVarForIndex(getVarsForBlock(b1)[i]) ), seeedpool->getVal(getMasterconss()[mc], getVarsForBlock(b1)[i]), SCIPvarGetName( seeedpool->getVarForIndex(getVarsForBlock(b2)[i])), seeedpool->getVal(getMasterconss()[mc], getVarsForBlock(b2)[i])  );
            return;
         }
      }

      /** variables seem to be identical so far */
      varmap[getVarsForBlock(b2)[i]] = getVarsForBlock(b1)[i];
   }

   for( int i = 0; i < getNConssForBlock(b1); ++i )
   {
      int cons1id;
      int cons2id;
      SCIP_CONS* cons1;
      SCIP_CONS* cons2;
      SCIP_Real* vals1;
      SCIP_Real* vals2;
      int nvals1;
      int nvals2;

      cons1id = getConssForBlock(b1)[i];
      cons2id = getConssForBlock(b2)[i];

      cons1 = seeedpool->getConsForIndex(cons1id);
      cons2 = seeedpool->getConsForIndex(cons2id);

      if( seeedpool->getNVarsForCons(cons1id) != seeedpool->getNVarsForCons(cons2id) )
      {
         SCIPdebugMessage("--> nvars differs for cons %s and cons %s!\n", SCIPconsGetName(cons1), SCIPconsGetName(cons2));
         return;
      }

      if( !SCIPisEQ(scip, GCGconsGetLhs(scip, cons1), GCGconsGetLhs(scip, cons2) ) )
      {
         SCIPdebugMessage("--> lhs differs for cons %s and cons %s!\n", SCIPconsGetName(cons1), SCIPconsGetName(cons2));
         return;
      }

      if( !SCIPisEQ(scip, GCGconsGetRhs(scip, cons1), GCGconsGetRhs(scip, cons2) ) )
      {
         SCIPdebugMessage("--> rhs differs for cons %s and cons %s!\n", SCIPconsGetName(cons1), SCIPconsGetName(cons2));
         return;
      }

      nvals1 = GCGconsGetNVars(scip, cons1);
      nvals2 = GCGconsGetNVars(scip, cons2);
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vals1, nvals1) );
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &vals2, nvals2) );
      GCGconsGetVals(scip, cons1, vals1, nvals1);
      GCGconsGetVals(scip, cons2, vals2, nvals2);


      if( !realArraysAreEqual(scip, vals1, nvals1, vals2, nvals2) )
       {
          SCIPdebugMessage("--> coefs differ for cons %s and cons %s!\n", SCIPconsGetName(cons1), SCIPconsGetName(cons2));
          SCIPfreeBufferArray(scip, &vals1);
           SCIPfreeBufferArray(scip, &vals2);
          return;
       }

      for( int v = 0; v < seeedpool->getNVarsForCons(cons1id) ; ++v )
      {
         if( varmap[seeedpool->getVarsForCons(cons2id)[v]] != seeedpool->getVarsForCons(cons1id)[v])
         {
            SCIPfreeBufferArray(scip, &vals1);
             SCIPfreeBufferArray(scip, &vals2);
            SCIPdebugMessage("--> vars differ for cons %s and cons %s!\n", SCIPconsGetName(cons1), SCIPconsGetName(cons2));
            return;
         }
      }



      SCIPfreeBufferArray(scip, &vals1);
      SCIPfreeBufferArray(scip, &vals2);

   }


   varmap = std::vector<int>(getNVarsForBlock(b1), -1);
   for( int i = 0; i < getNVarsForBlock(b1); ++i )
      varmap[i] = i;

   *identical = TRUE;
   return;
}

void Seeed::calcNCoeffsForBlockForMastercons(
   )
{
   ncoeffsforblockformastercons = std::vector<std::vector<int>>(getNBlocks());

   for( int b = 0; b < getNBlocks(); ++b )
      ncoeffsforblockformastercons[b] = std::vector<int>(getNMasterconss(), 0);

   for( int mc = 0; mc < getNMasterconss(); ++mc )
   {
      int cons = getMasterconss()[mc];
      for ( int vmc = 0; vmc < seeedpool->getNVarsForCons(cons); ++vmc )
      {
         int var = seeedpool->getVarsForCons(cons)[vmc];
         for( int b = 0; b < getNBlocks(); ++b )
         {
            if( isVarBlockvarOfBlock(var, b) )
               ++ncoeffsforblockformastercons[b][mc];
         }
      }
   }
   return;
}


SCIP_RETCODE Seeed::checkIdenticalBlocksTrivial(
   int                  b1,
   int                  b2,
   SCIP_Bool*           notidentical)
{

   if( getNConssForBlock(b1) != getNConssForBlock(b2) )
     {
        SCIPdebugMessage("--> number of constraints differs!\n");
        *notidentical = TRUE;
        return SCIP_OKAY;
     }


     if( getNVarsForBlock(b1) != getNVarsForBlock(b2) )
     {
        SCIPdebugMessage("--> number of variables differs!\n");
        *notidentical = TRUE;
        return SCIP_OKAY;
     }

     if( getNCoeffsForBlock(b1) != getNCoeffsForBlock( b2) )
     {
        SCIPdebugMessage("--> number of nonzero coeffs differs!\n");
        *notidentical = TRUE;
        return SCIP_OKAY;
     }

     if( ncoeffsforblockformastercons.size() == 0 )
         calcNCoeffsForBlockForMastercons();

     for( int mc = 0; mc < getNMasterconss(); ++mc )
     {
        if ( ncoeffsforblockformastercons[b1][mc] != ncoeffsforblockformastercons[b2][mc] )
        {
           SCIPdebugMessage("--> number of nonzero coeffs in %d-th master cons differs!\n", mc);
           *notidentical = TRUE;
           return SCIP_OKAY;
        }
     }


   return SCIP_OKAY;
}



/** assigns all open constraints and open variables
 *  strategy: assigns all conss and vars to the same block if they are indirectly connected
 *  a cons and a var are directly connected if the var appears in the cons */
SCIP_RETCODE Seeed::completeByConnected(
     )
{

   int cons;
   int var;

   changedHashvalue = true;

   /** tools to check if the openVars can still be found in a constraint yet */
   std::vector<int> varInBlocks; /** stores, in which block the variable can be found */

   /** tools to update openVars */
   std::vector<int> openvarsToDelete( 0 );
   std::vector<int> oldOpenconss;

   std::vector<bool> isConsOpen( nConss, false );
   std::vector<bool> isConsVisited( nConss, false );

   std::vector<bool> isVarOpen( nVars, false );
   std::vector<bool> isVarVisited( nVars, false );

   std::queue<int> helpqueue = std::queue<int>();
   std::vector<int> neighborConss( 0 );
   std::vector<int> neighborVars( 0 );

   assert( (int) conssForBlocks.size() == nBlocks );
   assert( (int) varsForBlocks.size() == nBlocks );
   assert( (int) stairlinkingVars.size() == nBlocks );

   SCIP_CALL( refineToMaster( ) );


   if( nBlocks < 0 )
      nBlocks = 0;

   /** initialize data structures */
   for( size_t c = 0; c < openConss.size(); ++ c )
   {
      cons = openConss[c];
      isConsOpen[cons] = true;
   }

   for( size_t v = 0; v < openVars.size(); ++ v )
   {
      var = openVars[v];
      isVarOpen[var] = true;
   }

   /** do breadth first search to find connected conss and vars */
   while( ! openConss.empty() )
   {
      int newBlockNr;

      assert( helpqueue.empty() );
      helpqueue.push( openConss[0] );
      neighborConss.clear();
      neighborConss.push_back( openConss[0] );
      isConsVisited[openConss[0]] = true;
      neighborVars.clear();

      while( ! helpqueue.empty() )
      {
         int nodeCons = helpqueue.front();
         assert( isConsOpencons( nodeCons ) );
         helpqueue.pop();
         for( int v = 0; v < seeedpool->getNVarsForCons( nodeCons ); ++ v )
         {
            var = seeedpool->getVarsForCons( nodeCons )[v];
            assert( isVarOpenvar( var ) || isVarLinkingvar( var ) );

            if( isVarVisited[var] || isVarLinkingvar( var ) )
               continue;

            for( int c = 0; c < seeedpool->getNConssForVar( var ); ++ c )
            {
               int otherNodeCons = seeedpool->getConssForVar( var )[c];
               if( ! isConsOpen[otherNodeCons] || isConsVisited[otherNodeCons] )
               {
                  continue;
               }
               assert( isConsOpencons( otherNodeCons ) );
               isConsVisited[otherNodeCons] = true;
               neighborConss.push_back( otherNodeCons );
               helpqueue.push( otherNodeCons );
            }
            isVarVisited[var] = true;
            neighborVars.push_back( var );
         }
      }

      /** assign found conss and vars to a new block */
      newBlockNr = getNBlocks() + 1;
      setNBlocks( newBlockNr );
      for( size_t i = 0; i < neighborConss.size(); ++ i )
      {
         cons = neighborConss[i];

         assert( isConsOpencons( cons ) );
         setConsToBlock( cons, newBlockNr - 1 );

         deleteOpencons( cons );
      }
      for( size_t i = 0; i < neighborVars.size(); ++ i )
      {
         var = neighborVars[i];
         setVarToBlock( var, newBlockNr - 1 );
         assert( isVarOpenvar( var ) );
         deleteOpenvar( var );
      }
   }

   /** assign left open vars to block 0, if it exists, and to master, otherwise */
   for( size_t i = 0; i < openVars.size(); ++ i )
   {
      var = openVars[i];
      if( getNBlocks() != 0 )
         setVarToBlock( var, 0 );
      else
         setVarToMaster( var );
      openvarsToDelete.push_back( var );
   }

   for( size_t i = 0; i < openvarsToDelete.size(); ++ i )
   {
      var = openvarsToDelete[i];
      deleteOpenvar( var );
   }

   assert( openConss.empty() );
   assert( openVars.empty() );

   sort();

   assert( checkConsistency( ) );

   return SCIP_OKAY;
}


/** try to reassign each  mastercons to one block without inducing conflicts  */
 SCIP_RETCODE Seeed::postprocessMasterToBlocks(
    SCIP_Bool* success
    )
 {
    *success = FALSE;
    return SCIP_OKAY;
 }


 /** try to reassign each  mastercons to one block without inducing conflicts  */
 SCIP_RETCODE Seeed::postprocessMasterToBlocksConssAdjacency(
    SCIP_Bool* success
    )
 {
    *success = FALSE;
    std::vector<int> constoreassign(0);
    std::vector<int> blockforconstoreassign(0);

    sort();

    std::vector<int> blockforvar(getNVars(), -1 );


    /**  */
    for( int b = 0; b < getNBlocks(); ++b )
    {
       for( size_t j  = 0; j < (size_t) getNVarsForBlock(b); ++j )
       {
          blockforvar[getVarsForBlock(b)[j] ] = b;
       }
    }


    for( int mc = 0; mc < getNMasterconss(); ++mc )
    {
       int masterconsid = getMasterconss()[mc];
       int hittenblock  = -1;

       SCIP_Bool hitsmastervar = FALSE;
       SCIP_Bool varhitsotherblock = FALSE;

       for( int var = 0; var < seeedpool->getNVarsForCons(masterconsid); ++var )
       {
          int varid = seeedpool->getVarsForCons(masterconsid)[var];
          if( isvarmaster[varid] )
          {
             hitsmastervar = TRUE;
             break;
          }

          if ( blockforvar[varid] != -1 )
          {
             if( hittenblock == -1 )
                hittenblock = blockforvar[varid];
             else if( hittenblock != blockforvar[varid] )
             {
                varhitsotherblock = TRUE;
                break;
             }
          }
       }

       if( hitsmastervar || varhitsotherblock )
          continue;

       if ( hittenblock != -1 )
       {
          constoreassign.push_back(masterconsid);
          blockforconstoreassign.push_back(hittenblock);
       }
    }


    for( size_t i = 0; i < constoreassign.size() ; ++i )
    {
          std::vector<int>::iterator todelete = lower_bound( masterConss.begin(), masterConss.end(), constoreassign[i] );
          masterConss.erase(todelete);

          conssForBlocks[blockforconstoreassign[i]].push_back( constoreassign[i] );
          conssforblocksorted = false;

    }

    if( constoreassign.size() > 0 )
       *success = SCIP_SUCCESS;

    sort();

    getScore( SCIPconshdlrDecompGetCurrScoretype( scip ) ) ;
    calcHashvalue();

    return SCIP_OKAY;
 }


/** assigns all open constraints and open variables
  *  strategy: assigns all conss same block if they are connected
  *  two constraints are adjacent if there is a common variable
  *  this relies on the consadjacency structure of the seeedpool
  *  hence it cannot be applied in presence of linking variables */
 SCIP_RETCODE Seeed::completeByConnectedConssAdjacency(
    ){

    int cons;
    int var;

    changedHashvalue = true;


    /** tools to check if the openVars can still be found in a constraint yet */
    std::vector<int> varInBlocks; /** stores, in which block the variable can be found */

    /** tools to update openVars */
    std::vector<int> oldOpenconss;
    std::vector<int> openvarsToDelete;

    if( getNLinkingvars() != 0 )
       return completeByConnected();

    std::vector<bool> isConsOpen( nConss, false );
    std::vector<bool> isConsVisited( nConss, false );

    varInBlocks = std::vector<int>(nVars, -1);

    std::queue<int> helpqueue = std::queue<int>();
    std::vector<int> neighborConss( 0 );

    assert( (int) conssForBlocks.size() == nBlocks );
    assert( (int) varsForBlocks.size() == nBlocks );
    assert( (int) stairlinkingVars.size() == nBlocks );

    SCIP_CALL( refineToMaster( ) );

    assert(checkConsistency() );

    if( nBlocks < 0 )
       nBlocks = 0;

    /** initialize data structures */
    for( size_t c = 0; c < openConss.size(); ++ c )
    {
       cons = openConss[c];
       isConsOpen[cons] = true;
    }

    /** do breadth first search to find connected conss */
    while( ! openConss.empty() )
    {
       int newBlockNr;

       assert( helpqueue.empty() );
       helpqueue.push( openConss[0] );
       neighborConss.clear();
       neighborConss.push_back( openConss[0] );
       isConsVisited[openConss[0]] = true;

       while( ! helpqueue.empty() )
       {
          int nodeCons = helpqueue.front();
          assert( isConsOpencons( nodeCons ) );
          helpqueue.pop();
          for( int c = 0; c < seeedpool->getNConssForCons( nodeCons ); ++ c )
          {
             int othercons = seeedpool->getConssForCons( nodeCons )[c];

             if( isConsVisited[othercons] || isConsMastercons( othercons ) || ! isConsOpen[othercons] )
                continue;

             assert( isConsOpencons( othercons ) );
             isConsVisited[othercons] = true;
             neighborConss.push_back( othercons );
             helpqueue.push( othercons );
          }
       }

       /** assign found conss and vars to a new block */
       newBlockNr = getNBlocks() + 1;
       setNBlocks( newBlockNr );
       for( size_t i = 0; i < neighborConss.size(); ++ i )
       {
          cons = neighborConss[i];

          assert( isConsOpencons( cons ) );
          setConsToBlock( cons, newBlockNr - 1 );
          deleteOpencons( cons );

          for( int j = 0; j < seeedpool->getNVarsForCons(cons); ++ j )
          {
             int newvar = seeedpool->getVarsForCons(cons)[j];

             if( isVarLinkingvar(newvar) || varInBlocks[newvar] != -1 )
                continue;

             assert(! isVarMastervar( newvar) );
             setVarToBlock( newvar, newBlockNr - 1 );
             varInBlocks[newvar] = newBlockNr - 1;
             if( isVarOpenvar(newvar) )
                deleteOpenvar( newvar );
          }

       }

    }

    /** assign left open vars to block 0, if it exists, and to master, otherwise */
    for( size_t i = 0; i < openVars.size(); ++ i )
    {
       var = openVars[i];
       if( getNBlocks() != 0 )
          setVarToBlock( var, 0 );
       else
          setVarToMaster( var );
       openvarsToDelete.push_back( var );
    }

    for( size_t i = 0; i < openvarsToDelete.size(); ++ i )
    {
       var = openvarsToDelete[i];
       deleteOpenvar( var );
    }

    assert( openConss.empty() );
    assert( openVars.empty() );

    sort();

    assert( checkConsistency( ) );

    return SCIP_OKAY;
 }


 /** assigns all open constraints and open variables
   *  strategy: assigns all conss same block if they are connected
   *  two constraints are adjacent if there is a common variable
   *  this relies on the consadjacency structure of the seeedpool
   *  hence it cannot be applied in presence of linking variables */
  SCIP_RETCODE Seeed::assignSmallestComponentsButOneConssAdjacency(
     ){

     int cons;
     SCIP_Bool conssadjcalculated;

     changedHashvalue = true;


     /** tools to check if the openVars can still be found in a constraint yet */
     std::vector<int> varInBlocks; /** stores, in which block the variable can be found */

     /** tools to update openVars */
     std::vector<int> oldOpenconss;
     std::vector<int> openvarsToDelete;

     if( getNLinkingvars() != 0 )
        return completeByConnected();

     SCIPgetBoolParam(scip, "detection/conssadjcalculated", &conssadjcalculated);

     if( !conssadjcalculated )
     {
        seeedpool->createConssAdjacency();
        SCIPsetBoolParam(scip, "detection/conssadjcalculated", TRUE);
     }

     std::vector<bool> isConsOpen( nConss, false );
     std::vector<bool> isConsVisited( nConss, false );

     std::vector<std::vector<int>> conssfornewblocks(0);
     std::vector<std::vector<int>> varsfornewblocks(0);

     std::vector<int> constoconsider;
     int nnewblocks;
     int largestcomponent;
     int sizelargestcomponent;

     constoconsider = openConss;

     varInBlocks = std::vector<int>(nVars, -1);
     nnewblocks = 0;
     largestcomponent = -1;
     sizelargestcomponent = 0;

     std::queue<int> helpqueue = std::queue<int>();
     std::vector<int> neighborConss( 0 );

     assert( (int) conssForBlocks.size() == nBlocks );
     assert( (int) varsForBlocks.size() == nBlocks );
     assert( (int) stairlinkingVars.size() == nBlocks );

     assert(checkConsistency() );

     if( nBlocks < 0 )
        nBlocks = 0;

     /** initialize data structures */
     for( size_t c = 0; c < openConss.size(); ++ c )
     {
        cons = openConss[c];
        isConsOpen[cons] = true;
     }

     /** do breadth first search to find connected conss */
     while( ! constoconsider.empty() )
     {
        std::vector<int> newconss(0);
        std::vector<int> newvars(0);

        assert( helpqueue.empty() );
        helpqueue.push( constoconsider[0] );
        neighborConss.clear();
        neighborConss.push_back( constoconsider[0] );
        isConsVisited[constoconsider[0]] = true;

        while( ! helpqueue.empty() )
        {
           int nodeCons = helpqueue.front();
           assert( isConsOpencons( nodeCons ) );
           helpqueue.pop();
           for( int c = 0; c < seeedpool->getNConssForCons( nodeCons ); ++ c )
           {
              int othercons;
              othercons = seeedpool->getConssForCons( nodeCons )[c];

              if( isConsVisited[othercons] || isConsMastercons( othercons ) || ! isConsOpen[othercons] )
                 continue;

              assert( isConsOpencons( othercons ) );
              isConsVisited[othercons] = true;
              neighborConss.push_back( othercons );
              helpqueue.push( othercons );
           }
        }

        /** assign found conss and vars to a new block */
        ++nnewblocks;
        for( size_t i = 0; i < neighborConss.size(); ++ i )
        {
           std::vector<int>::iterator consiter;
           cons = neighborConss[i];
           consiter = std::lower_bound(constoconsider.begin(), constoconsider.end(), cons);
           assert(consiter != constoconsider.end() );
           constoconsider.erase(consiter);
           assert( isConsOpencons( cons ) );
           newconss.push_back(cons);

           for( int j = 0; j < seeedpool->getNVarsForCons(cons); ++ j )
           {
              int newvar = seeedpool->getVarsForCons(cons)[j];

              if( isVarLinkingvar(newvar) || varInBlocks[newvar] != -1 )
                 continue;

              assert(! isVarMastervar( newvar) );
              newvars.push_back(newvar);
              varInBlocks[newvar] = nnewblocks;
           }
        }
        conssfornewblocks.push_back(newconss);
        varsfornewblocks.push_back(newvars);
     }

     for( int i = 0; i < nnewblocks; ++i )
     {
        if( (int)conssfornewblocks[i].size() > sizelargestcomponent )
        {
           sizelargestcomponent = (int)conssfornewblocks[i].size();
           largestcomponent = i;
        }
     }

     if( nnewblocks > 1 )
     {
        int oldnblocks;
        bool largestdone = false;
        oldnblocks = getNBlocks();
        setNBlocks(nnewblocks - 1 + getNBlocks());

        for( int i = 0; i < nnewblocks; ++i)
        {
           if( i == largestcomponent )
           {
              largestdone = true;
              continue;
           }
           for( int c = 0; c < (int) conssfornewblocks[i].size() ; ++c)
           {
              bookAsBlockCons(conssfornewblocks[i][c], oldnblocks + i - (largestdone ? 1 : 0) );
           }

           for( int v = 0; v < (int) varsfornewblocks[i].size() ; ++v )
           {
              bookAsBlockVar(varsfornewblocks[i][v], oldnblocks + i - (largestdone ? 1 : 0) );
           }
        }

        flushBooked();

        sort();
     }

     assert( checkConsistency( ) );

     return SCIP_OKAY;
  }




/** assigns all open constraints and open variables
 *  strategy: assigns a cons (and related vars) to any block if possible by means of prior var assignments
 *  and to master, if there does not exist such a block */
SCIP_RETCODE Seeed::completeGreedily(
   )
{
   bool checkVar;
   bool varInBlock;
   bool notassigned;

   changedHashvalue = true;

   /** tools to check if the openVars can still be found in a constraint yet*/
   std::vector<int> varInBlocks; /** stores, in which block the variable can be found */

   assert( (int) conssForBlocks.size() == nBlocks );
   assert( (int) varsForBlocks.size() == nBlocks );
   assert( (int) stairlinkingVars.size() == nBlocks );

   if( nBlocks == 0 && openConss.size() > 0 )
   {
      addBlock();
      if( openConss.size() != 0 )
      {
         setConsToBlock( openConss[0], 0 );
         openConss.erase( openConss.begin() );
      }
      else if( masterConss.size() != 0 )
      {
         setConsToBlock( masterConss[0], 0 );
         masterConss.erase( masterConss.begin() );
         isconsmaster[masterConss[0]] = false;
      }
      else
         assert( ! ( openConss.size() == 0 && masterConss.size() == 0 ) );
   }

   /** check if the openVars can already be found in a constraint */
   for( size_t i = 0; i < openVars.size(); ++ i )
   {/** assigns all open constraints and open variables
    *  strategy: assign all conss and vars that are indirectly connected to the same block
    *  a cons and a var are directly connected if the var appears in the cons */
      varInBlocks.clear();

      /** test if the variable can be found in blocks */
      for( int b = 0; b < nBlocks; ++ b )
      {
         varInBlock = false;
         for( size_t k = 0; k < conssForBlocks[b].size() && ! varInBlock; ++ k )
         {
            for( int l = 0; l < seeedpool->getNVarsForCons( conssForBlocks[b][k] ); ++ l )
            {
               if( openVars[i] == seeedpool->getVarsForCons( conssForBlocks[b][k] )[l] )
               {
                  varInBlocks.push_back( b );
                  varInBlock = true;
                  break;
               }
            }
         }
      }
      if( varInBlocks.size() == 1 ) /** if the variable can be found in one block set the variable to a variable of the block*/
      {
         bookAsBlockVar( openVars[i], varInBlocks[0] );
         continue; /** the variable does'nt need to be checked any more */
      }
      else if( varInBlocks.size() == 2 ) /** if the variable can be found in two blocks check if it is a linking var or a stairlinking var*/
      {
         if( varInBlocks[0] + 1 == varInBlocks[1] )
         {
            bookAsStairlinkingVar( openVars[i], varInBlocks[0] );
            continue; /** the variable does'nt need to be checked any more */
         }
         else
         {
            bookAsLinkingVar( openVars[i] );
            continue; /** the variable does'nt need to be checked any more */
         }
      }
      else if( varInBlocks.size() > 2 ) /** if the variable can be found in more than two blocks it is a linking var */
      {
         bookAsLinkingVar( openVars[i] );
         continue; /** the variable does'nt need to be checked any more */
      }

      checkVar = true;

      /** if the variable can be found in an open constraint it is still an open var */
      for( size_t j = 0; j < openConss.size(); ++ j )
      {
         checkVar = true;
         for( int k = 0; k < seeedpool->getNVarsForCons( j ); ++ k )
         {
            if( openVars[i] == seeedpool->getVarsForCons( j )[k] )
            {
               checkVar = false;
               break;
            }
         }
         if( ! checkVar )
         {
            break;
         }
      }

      /** test if the variable can be found in a master constraint yet */
        for( int k = 0; k < seeedpool->getNConssForVar( openVars[i] ) && checkVar; ++ k )
        {
           if( isconsmaster[seeedpool->getConssForVar(openVars[i])[k]] )
           {
              bookAsMasterVar( openVars[i] );
              checkVar = false; /** the variable does'nt need to be checked any more */
              break;
           }
        }

//      for( size_t j = 0; j < masterConss.size() && checkVar; ++ j )
//      {
//         for( int k = 0; k < seeedpool->getNVarsForCons( masterConss[j] ); ++ k )
//         {
//            if( openVars[i] == seeedpool->getVarsForCons( masterConss[j] )[k] )
//            {
//               bookAsMasterVar( openVars[i] );
//               checkVar = false; /** the variable does'nt need to be checked any more */
//               break;
//            }
//         }
//      }
   }

   flushBooked();

   /** assign open conss greedily */
   for( size_t i = 0; i < openConss.size(); ++ i )
   {
      std::vector<int> vecOpenvarsOfBlock; /** stores the open vars of the blocks */
      bool consGotBlockcons = false; /** if the constraint can be assigned to a block */

      /** check if the constraint can be assigned to a block */
      for( int j = 0; j < nBlocks; ++ j )
      {
         /** check if all vars of the constraint are a block var of the current block, an open var, a linkingvar or a mastervar*/
         consGotBlockcons = true;
         for( int k = 0; k < seeedpool->getNVarsForCons( openConss[i] ); ++ k )
         {
            if( isVarBlockvarOfBlock( seeedpool->getVarsForCons( openConss[i] )[k], j )
               || isVarOpenvar( seeedpool->getVarsForCons( openConss[i] )[k] )
               || isVarLinkingvar( seeedpool->getVarsForCons( openConss[i] )[k] )
               || isVarStairlinkingvarOfBlock( seeedpool->getVarsForCons( openConss[i] )[k], j )
               || ( j != 0 && isVarStairlinkingvarOfBlock( seeedpool->getVarsForCons( openConss[i] )[k], j - 1 ) ) )
            {
               if( isVarOpenvar( seeedpool->getVarsForCons( openConss[i] )[k] ) )
               {
                  vecOpenvarsOfBlock.push_back( seeedpool->getVarsForCons( openConss[i] )[k] );
               }
            }
            else
            {
               vecOpenvarsOfBlock.clear(); /** the open vars don't get vars of the block */
               consGotBlockcons = false; /** the constraint can't be constraint of the block, check the next block */
               break;
            }
         }
         if( consGotBlockcons ) /** the constraint can be assigned to the current block */
         {
            bookAsBlockCons( openConss[i], j );
            for( size_t k = 0; k < vecOpenvarsOfBlock.size(); ++ k ) /** the openvars in the constraint get block vars */
            {
               setVarToBlock( vecOpenvarsOfBlock[k], j );
               deleteOpenvar( vecOpenvarsOfBlock[k] );
            }
            vecOpenvarsOfBlock.clear();

            break;
         }
      }

      if( ! consGotBlockcons ) /** the constraint can not be assigned to a block, set it to master */
         bookAsMasterCons( openConss[i] );
   }

   flushBooked();

   /** assign open vars greedily */
   for( size_t i = 0; i < openVars.size(); ++ i )
   {
      notassigned = true;
      for( size_t j = 0; j < masterConss.size() && notassigned; ++ j )
      {
         for( int k = 0; k < seeedpool->getNVarsForCons( masterConss[j] ); ++ k )
         {
            if( openVars[i] == seeedpool->getVarsForCons( masterConss[j] )[k] )
            {
               bookAsMasterVar( openVars[i] );
               notassigned = false;
               break;
            }
         }
      }
   }

   flushBooked();

   /** check if the open conss are all assigned */
   if( ! checkAllConssAssigned() )
   {
      SCIPwarningMessage(scip,"ERROR: Something went wrong, there are still open cons, although all should have been assigned\n" );
      assert( false ); /** assigns all open constraints and open variables
       *  strategy: assign all conss and vars to the same block if they are indirectly connected
       *  a cons and a var are directly connected if the var appears in the cons */
   }

   /** check if the open vars are all assigned */
   if( ! openVars.empty() )
   {
      SCIPwarningMessage(scip,"ERROR: Something went wrong, there are still open vars, although all should have been assigned\n" );
      assert( false );
   }

   assert( checkConsistency( ) );

   return SCIP_OKAY;
}

/** returns true if the given detector used a consclassifier */
bool Seeed::consClassifierUsed(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) usedClassifier.size() );

   return ( usedClassifier[detectorchainindex] != NULL )
      && ( dynamic_cast<ConsClassifier*>( usedClassifier[detectorchainindex] ) != NULL );
}

/** assigns every open cons/var
 *  - to the respective block if it hits exactly one blockvar/blockcons and no open vars/conss
 *  - to master/linking if it hits blockvars/blockconss assigned to different blocks
 *  - and every cons to master that hits a master var
 *  - and every var to master if it does not hit any blockcons and has no open cons */
SCIP_RETCODE Seeed::considerImplicits(
   )
{
   int cons;
   int var;
   std::vector<int> blocksOfBlockvars; /** blocks with blockvars which can be found in the cons */
   std::vector<int> blocksOfOpenvar; /** blocks in which the open var can be found */
   bool master;
   bool hitsOpenVar;
   bool hitsOpenCons;
   bool benders;

   benders = seeedpool->isForBenders();

   changedHashvalue = true;

   sort();

   /** set openConss with more than two blockvars to master */
   for( size_t c = 0; c < openConss.size(); ++ c )
   {
      std::vector<bool> hitsblock = std::vector<bool>(nBlocks, false);
      blocksOfBlockvars.clear();
      master = false;
      hitsOpenVar = false;
      cons = openConss[c];

      for( int v = 0; v < seeedpool->getNVarsForCons( cons ) && ! master; ++ v )
      {
         var = seeedpool->getVarsForCons( cons )[v];

         if( isVarMastervar( var ) )
         {
            master = true;
            bookAsMasterCons( cons );
            continue;
         }

         if( isVarOpenvar( var ) )
         {
            hitsOpenVar = true;
            if( !benders )
               continue;
         }


         for( int b = 0; b < nBlocks && ! master; ++ b )
         {
            if( isVarBlockvarOfBlock( var, b ) && !hitsblock[b] )
            {
               hitsblock[b] = true;
               blocksOfBlockvars.push_back( b );
               break;
            }
         }
      }

      if ( benders && blocksOfBlockvars.size() == 1 && !master )
         bookAsBlockCons( cons, blocksOfBlockvars[0] );

//      if( benders && master )
//         bookAsMasterCons( cons );

      if( !benders && blocksOfBlockvars.size() > 1 )
         bookAsMasterCons( cons );

      /* also assign open constraints that only have vars assigned to one single block and no open vars*/
      if( blocksOfBlockvars.size() == 1 && ! hitsOpenVar && ! master && ! benders )
         bookAsBlockCons( cons, blocksOfBlockvars[0] );
   }

   flushBooked();

   /** set open var to linking, if it can be found in more than one block or set it to a block if it has only constraints in that block and no open constraints or set it to master if it only hits master constraints */
   for( size_t i = 0; i < openVars.size(); ++ i )
   {
      std::vector<bool> hitsblock = std::vector<bool>(nBlocks, false);
      bool hitsmasterconss = false;
      bool hitonlymasterconss = true;
      bool hitonlyblockconss = true;
      blocksOfOpenvar.clear();
      var = openVars[i];
      hitsOpenCons = false;


      for( int c = 0; c < seeedpool->getNConssForVar( var ); ++ c )
      {
         cons = seeedpool->getConssForVar( var )[c];
         if ( isConsMastercons(cons) )
         {
            hitsmasterconss = true;
            hitonlyblockconss = false;
            continue;
         }

         if( isConsOpencons( cons ) )
         {
            hitsOpenCons = true;
            hitonlyblockconss = false;
            hitonlymasterconss = false;
            continue;
         }
      }
      for( int b = 0; b < nBlocks; ++ b )
      {
         for( int c = 0; c < seeedpool->getNConssForVar( var ); ++ c )
         {
            cons = seeedpool->getConssForVar( var )[c];
            if( isConsBlockconsOfBlock( cons, b ) && !hitsblock[b] )
            {
               hitsblock[b] = true;
               hitonlymasterconss = false;
               blocksOfOpenvar.push_back( b );
               break;
            }
         }
      }

      if( blocksOfOpenvar.size() > 1 )
      {
         bookAsLinkingVar( var );
         continue;
      }

      if( benders && blocksOfOpenvar.size() == 1 &&  hitsmasterconss )
      {
         bookAsLinkingVar( var );
      }

      if( benders && hitonlyblockconss && blocksOfOpenvar.size() > 0 )
      {
         bookAsBlockVar( var, blocksOfOpenvar[0] );
      }

      if( benders && hitonlymasterconss)
      {
         bookAsMasterVar( var);
      }

      if( !benders && blocksOfOpenvar.size() == 1 && ! hitsOpenCons )
      {
         bookAsBlockVar( var, blocksOfOpenvar[0] );
      }

      if( !benders && blocksOfOpenvar.size() == 0 && ! hitsOpenCons )
      {
         bookAsMasterVar( var );
      }



   }

   flushBooked();

   return SCIP_OKAY;
}

/** copies the given seeed's classifier statistics */
SCIP_RETCODE Seeed::copyClassifierStatistics(
   const Seeed* otherseeed
   )
{
   usedClassifier = otherseeed->usedClassifier;
   classesToMaster = otherseeed->classesToMaster;
   classesToLinking = otherseeed->classesToLinking;

   return SCIP_OKAY;
}

/** deletes empty blocks */
SCIP_RETCODE Seeed::deleteEmptyBlocks(
   bool variables                            /* if TRUE a block is only considered to be empty if it contains neither constraints or variables */
   )
{
   bool emptyBlocks = true;
   bool benders = false;
   int block = - 1;
   int b;

   changedHashvalue = true;

   assert( (int) conssForBlocks.size() == nBlocks );
   assert( (int) varsForBlocks.size() == nBlocks );
   assert( (int) stairlinkingVars.size() == nBlocks );

   benders =  seeedpool->isForBenders();

   while( emptyBlocks )
   {
      emptyBlocks = false;
      for( b = 0; b < nBlocks; ++ b )
      {
         if( conssForBlocks[b].size() == 0 &&  ( variables ? varsForBlocks[b].size() == 0 : true) )
         {
            emptyBlocks = true;
            block = b;
         }
         if( benders && ( conssForBlocks[b].size() == 0 ||  varsForBlocks[b].size() == 0) )
         {
            emptyBlocks = true;
            block = b;
         }

      }
      if( emptyBlocks )
      {
         nBlocks --;

         std::vector<std::vector<int>>::iterator it;

         it = stairlinkingVars.begin();
         for( b = 0; b < block; ++ b )
            it ++;
         stairlinkingVars.erase( it );

         it = conssForBlocks.begin();
         for( b = 0; b < block; ++ b )
            it ++;
         for( size_t j = 0; j < conssForBlocks[block].size(); ++j )
         {
            masterConss.push_back(conssForBlocks[block][j]);
            isconsmaster[conssForBlocks[block][j]] = true;
         }
         std::sort( masterConss.begin(), masterConss.end() );
         conssForBlocks.erase( it );

         it = varsForBlocks.begin();
         for( b = 0; b < block; ++ b )
            it ++;
         for( size_t j = 0; j < varsForBlocks[block].size(); ++ j )
         {
            masterVars.push_back( varsForBlocks[block][j] );
            isvarmaster[varsForBlocks[block][j]] = true;
         }
         varsForBlocks.erase( it );
         std::sort( masterVars.begin(), masterVars.end() );

         //set stairlinkingvars of the previous block to block vars
         if( block != 0 && (int) stairlinkingVars[block - 1].size() != 0 )
         {
            std::vector<int>::iterator iter = stairlinkingVars[block - 1].begin();
            std::vector<int>::iterator iterEnd = stairlinkingVars[block - 1].end();
            std::vector<int> stairlinkingVarsOfPreviousBlock;
            for( ; iter != iterEnd; ++ iter )
            {
               bookAsBlockVar( * iter, block - 1 );
               stairlinkingVarsOfPreviousBlock.push_back( * iter );
            }
            for( size_t i = 0; i < stairlinkingVarsOfPreviousBlock.size(); ++ i )
            {
               iter = find( stairlinkingVars[block - 1].begin(), stairlinkingVars[block - 1].end(),
                  stairlinkingVarsOfPreviousBlock[i] );
               assert( iter != stairlinkingVars[block - 1].end() );
               stairlinkingVars[block - 1].erase( iter );
            }
            flushBooked();
         }
      }
   }
   return SCIP_OKAY;
}

/** deletes a cons from list of open conss */
SCIP_RETCODE Seeed::deleteOpencons(
   int opencons
   )
{
   assert( opencons >= 0 && opencons < nConss );
   std::vector<int>::iterator it;
   it = lower_bound( openConss.begin(), openConss.end(), opencons );
   assert( it != openConss.end() && *it == opencons );
   openConss.erase( it );
   isconsopen[opencons] = false;
   changedHashvalue = true;

   return SCIP_OKAY;
}

/** deletes a var from the list of open vars */
SCIP_RETCODE Seeed::deleteOpenvar(
   int openvar
   )
{
   assert( openvar >= 0 && openvar < nVars );
   std::vector<int>::iterator it;
   it = lower_bound( openVars.begin(), openVars.end(), openvar );
   assert( it != openVars.end() && *it == openvar );
   openVars.erase( it );
   isvaropen[openvar] = false;
   changedHashvalue = true;
   return SCIP_OKAY;
}

SCIP_RETCODE Seeed::displayAggregationInformation()
{
   if( !agginfocalculated )
   {
      SCIPinfoMessage(scip, NULL, " Aggregation information is not calculated yet \n ");
      return SCIP_OKAY;
   }

   SCIPinfoMessage(scip, NULL, " number of representative blocks: %d \n", nrepblocks);
   for( int i = 0; i < nrepblocks; ++i )
   {
      SCIPinfoMessage(scip, NULL, "representative block %d : ", i);

      for( size_t b = 0; b < reptoblocks[i].size(); ++b )
         SCIPinfoMessage(scip, NULL, "%d ", reptoblocks[i][b] );

      SCIPinfoMessage(scip, NULL, "\n", i);
   }

   return SCIP_OKAY;
}

/** displays the assignments of the conss */
SCIP_RETCODE Seeed::displayConss()
{
   for( int b = 0; b < nBlocks; ++ b )
   {
      if( getNConssForBlock( b ) != 0 )
      {
         std::cout << "constraint(s) in block " << b << ": ";
         std::cout << getConssForBlock( b )[0] << "|" << SCIPconsGetName(seeedpool->getConsForIndex(getConssForBlock( b )[0]) ) ;
         for( int c = 1; c < getNConssForBlock( b ); ++ c )
            std::cout << ", " << getConssForBlock( b )[c] << "|" << SCIPconsGetName(seeedpool->getConsForIndex(getConssForBlock( b )[c]) ) ;
         std::cout << "\n";
      }
      else
         std::cout << "0 constraints in block " << b << std::endl;
   }

   if( getNMasterconss() != 0 )
   {
      std::cout << "masterconstraint(s): ";
      std::cout << masterConss[0];
      for( int c = 1; c < getNMasterconss(); ++ c )
         std::cout << ", " << masterConss[c];
      std::cout << "\n";
   }
   else
      std::cout << "0 masterconstraints" << std::endl;

   if( getNOpenconss() != 0 )
   {
      std::cout << "open constraint(s): ";
      std::cout << openConss[0];
      for( int c = 1; c < getNOpenconss(); ++ c )
         std::cout << ", " << openConss[c];
      std::cout << "\n";
   }
   else
      std::cout << "0 open constraints" << std::endl;

   return SCIP_OKAY;
}

/** displays the relevant information of the seeed */
SCIP_RETCODE Seeed::displayInfo(
   int detailLevel
   )
{
   assert( seeedpool != NULL );
   assert( 0 <= detailLevel );


   std::cout << std::endl;

   /* general information */
   std::cout << "-- General information --" << std::endl;
   std::cout << " ID: " << id << std::endl;
   std::cout << " Hashvalue: " << hashvalue << std::endl;
   std::cout << " Score: " << score << std::endl;
   if( getNOpenconss() + getNOpenconss() > 0 )
      std::cout << " Maxwhitescore >= " << maxwhitescore << std::endl;
   else
      std::cout << " Max white score: " << maxwhitescore << std::endl;
   if( getNOpenconss() + getNOpenconss() == 0 )
         std::cout << " Max-foreseeing-white-score: " << maxforeseeingwhitescore << std::endl;
   if( getNOpenconss() + getNOpenconss() == 0 )
         std::cout << " Max-foreseeing-white-aggregated-score: " << maxforeseeingwhitescoreagg << std::endl;
   if( getNOpenconss() + getNOpenconss() == 0 )
         std::cout << " PPC-max-foreseeing-white-score: " <<  setpartfwhitescore << std::endl;

   if( getNOpenconss() + getNOpenconss() == 0 )
          std::cout << " PPC-max-foreseeing-white-aggregated-score: " <<  setpartfwhitescoreagg << std::endl;

   if( getNOpenconss() + getNOpenconss() == 0 )
          std::cout << " Benderborderarea: " << benderareascore << std::endl;


   if( getNOpenconss() + getNOpenconss() == 0 )
          std::cout << " Bendersscore: " << bendersscore << std::endl;

   if( getNOpenconss() + getNOpenconss() == 0 )
          std::cout << " blockarea: " << blockareascore << std::endl;


   if( getNOpenconss() + getNOpenconss() == 0 )
          std::cout << " borderareascore: " << borderareascore << std::endl;


   std::cout << " HassetppMaster: " << hasSetppMaster() << std::endl;
   std::cout << " HassetppcMaster: " << hasSetppcMaster() << std::endl;
   std::cout << " HassetppccardMaster: " << hasSetppccardMaster() << std::endl;
   std::cout << " Seeed is for the " << ( isfromunpresolved ? "unpresolved" : "presolved" ) << " problem and "
      << ( usergiven ? "usergiven" : "not usergiven" ) << "." << std::endl;
   std::cout << " Number of constraints: " << getNConss() << std::endl;
   std::cout << " Number of variables: " << getNVars() << std::endl;

   displayAggregationInformation();
   std::cout << std::endl;



   /* detection information */
   std::cout << "-- Detection and detectors --" << std::endl;
   std::cout << " Seeed stems from the " << ( stemsFromUnpresolved ? "unpresolved" : "presolved" ) << " problem." << std::endl;
   if( isFromLegacymode() )
   {
      std::cout << " Seeed is from a detector operating in legacymode." << std::endl;
   }

   /* ancestor seeeds' ids */
   std::cout << " IDs of ancestor seeeds: ";
   if( listofancestorids.size() > 0 )
      std::cout << listofancestorids[0];
   for( int i = 1; i < (int) listofancestorids.size(); ++i )
      std::cout << ", " << listofancestorids[i];
   std::cout << std::endl;

   /* detector chain information */
   std::cout << " " << getNDetectors() << " detector" << ( getNDetectors() > 1 ? "s" : "" ) << " worked on this seeed:";
   if( getNDetectors() != 0 )
   {
      std::string detectorrepres;

      if( detectorChain[0] == NULL )
         detectorrepres = "user";
      else
      {
         /* potentially add finisher label */
         detectorrepres = (
            getNDetectors() != 1 || !isFinishedByFinisher ? DECdetectorGetName(detectorChain[0]) :
               "(finish) " + std::string(DECdetectorGetName(detectorChain[0])));
      }

      if( detailLevel > 0 )
      {
         std::cout << std::endl << " 1.: " << detectorrepres << std::endl;
         std::cout << getDetectorStatistics( 0 );
         std::cout << getDetectorClassifierInfo( 0, detailLevel > 1 && ( !stemsFromUnpresolved || isfromunpresolved ) );
      }
      else
      {
         std::cout << " " << detectorrepres;
      }

      for( int d = 1; d < getNDetectors(); ++d )
      {
         /* potentially add finisher label */
         detectorrepres = (
            getNDetectors() != d + 1 || !isFinishedByFinisher ? DECdetectorGetName(detectorChain[d]) :
               "(finish) " + std::string(DECdetectorGetName(detectorChain[d])));


         if( detailLevel > 0 )
         {
            std::cout << " " << ( d + 1 ) << ".: " << detectorrepres << std::endl;
            std::cout << getDetectorStatistics( d );
            std::cout << getDetectorClassifierInfo( d, detailLevel > 1 && ( !stemsFromUnpresolved || isfromunpresolved ) );
         }
         else
         {
            std::cout << ", " << detectorrepres;
         }
      }

      if( detailLevel <= 0 )
      {
         std::cout << std::endl;
      }
   }

   std::cout << std::endl;

   /* variable information */
   std::cout << "-- Border and unassigned --" << std::endl;
   std::cout << " Linkingvariables";
   if( detailLevel > 1 )
   {
      std::cout << " (" << getNLinkingvars() << ")";
      if( getNLinkingvars() > 0 )
         std::cout << ":  " << SCIPvarGetName( seeedpool->getVarForIndex( getLinkingvars()[0] ) );
      for( int v = 1; v < getNLinkingvars(); ++v )
      {
         std::cout << ", " << SCIPvarGetName( seeedpool->getVarForIndex( getLinkingvars()[v] ) );
      }
      std::cout << std::endl;
   }
   else
   {
      std::cout << ": " << getNLinkingvars() << std::endl;
   }
   std::cout << " Masterconstraints";
   if( detailLevel > 1 )
   {
      std::cout << " (" << getNMasterconss() << ")";
      if( getNMasterconss() > 0 )
         std::cout << ":  " << SCIPconsGetName( seeedpool->getConsForIndex( getMasterconss()[0] ) );
      for( int c = 1; c < getNMasterconss(); ++c )
      {
         std::cout << ", " << SCIPconsGetName( seeedpool->getConsForIndex( getMasterconss()[c] ) );
      }
      std::cout << std::endl;
   }
   else
   {
      std::cout << ": " << getNMasterconss() << std::endl;
   }
   std::cout << " Mastervariables";
   if( detailLevel > 1 )
   {
      std::cout << " (" << getNMastervars() << ")";
      if( getNMastervars() > 0 )
         std::cout << ":  " << SCIPvarGetName( seeedpool->getVarForIndex( getMastervars()[0] ) );
      for( int v = 1; v < getNMastervars(); ++v )
      {
         std::cout << ", " << SCIPvarGetName( seeedpool->getVarForIndex( getMastervars()[v] ) );
      }
      std::cout << std::endl;
   }
   else
   {
      std::cout << ": " << getNMastervars() << std::endl;
   }
   std::cout << " Open constraints";
   if( detailLevel > 1 )
   {
      std::cout << " (" << getNOpenconss() << ")";
      if( getNOpenconss() > 0 )
         std::cout << ":  " << SCIPconsGetName( seeedpool->getConsForIndex( getOpenconss()[0] ) );
      for( int c = 1; c < getNOpenconss(); ++c )
      {
         std::cout << ", " << SCIPconsGetName( seeedpool->getConsForIndex( getOpenconss()[c] ) );
      }
      std::cout << std::endl;
   }
   else
   {
      std::cout << ": " << getNOpenconss() << std::endl;
   }
   std::cout << " Open variables";
   if( detailLevel > 1 )
   {
      std::cout << " (" << getNOpenvars() << ")";
      if( getNOpenvars() > 0 )
         std::cout << ":  " << SCIPvarGetName( seeedpool->getVarForIndex( getOpenvars()[0] ) );
      for( int v = 1; v < getNOpenvars(); ++v )
      {
         std::cout << ", " << SCIPvarGetName( seeedpool->getVarForIndex( getOpenvars()[v] ) );
      }
      std::cout << std::endl;
   }
   else
   {
      std::cout << ": " << getNOpenvars() << std::endl;
   }

   std::cout << std::endl;

   /* block information */
   std::cout << "-- Blocks --" << std::endl;
   std::cout << " Number of blocks: " << nBlocks << std::endl;

   if( detailLevel > 0 )
   {
      for( int b = 0; b < nBlocks; ++b )
      {
         std::cout << " Block " << b << ":" << std::endl;

         std::cout << "  Constraints";
         if( detailLevel > 1 )
         {
            std::cout << " (" << getNConssForBlock( b ) << ")";
            if( getNConssForBlock( b ) > 0 )
               std::cout << ":  " << SCIPconsGetName( seeedpool->getConsForIndex( getConssForBlock( b )[0] ) );
            for( int c = 1; c < getNConssForBlock( b ); ++c )
            {
               std::cout << ", " << SCIPconsGetName( seeedpool->getConsForIndex( getConssForBlock( b )[c] ) );
            }
            std::cout << std::endl;
         }
         else
         {
            std::cout << ": " << getNConssForBlock( b ) << std::endl;
         }

         std::cout << "  Variables";
         if( detailLevel > 1 )
         {
            std::cout << " (" << getNVarsForBlock( b ) << ")";
            if( getNVarsForBlock( b ) > 0 )
               std::cout << ":  " << SCIPvarGetName( seeedpool->getVarForIndex( getVarsForBlock( b )[0] ) );
            for( int v = 1; v < getNVarsForBlock( b ); ++v )
            {
               std::cout << ", " << SCIPvarGetName( seeedpool->getVarForIndex( getVarsForBlock( b )[v] ) );
            }
            std::cout << std::endl;
         }
         else
         {
            std::cout << ": " << getNVarsForBlock( b ) << std::endl;
         }

         std::cout << "  Stairlinkingvariables";
         if( detailLevel > 1 )
         {
            std::cout << " (" << getNStairlinkingvars( b ) << ")";
            if( getNStairlinkingvars( b ) > 0 )
               std::cout << ":  " << SCIPvarGetName( seeedpool->getVarForIndex( getStairlinkingvars( b )[0] ) );
            for( int v = 1; v < getNStairlinkingvars( b ); ++v )
            {
               std::cout << ", " << SCIPvarGetName( seeedpool->getVarForIndex( getStairlinkingvars( b )[v] ) );
            }
            std::cout << std::endl;
         }
         else
         {
            std::cout << ": " << getNStairlinkingvars( b ) << std::endl;
         }
      }
   }

   std::cout << std::endl;

   return SCIP_OKAY;
}

/** displays the relevant information of the seeed */
SCIP_RETCODE Seeed::displaySeeed(
   )
{
   std::cout << "ID: " << id << std::endl;
   std::cout << "number of blocks: " << nBlocks << std::endl;
   std::cout << "hashvalue: " << hashvalue << std::endl;
   std::cout << "score: " << score << std::endl;
   if( getNOpenconss() + getNOpenconss() > 0 )
      std::cout << "maxwhitescore >= " << maxwhitescore << std::endl;
   else
      std::cout << "maxwhitescore: " << maxwhitescore << std::endl;
   std::cout << "ancestorids: ";
   for( size_t i = 0; i < listofancestorids.size(); ++ i )
      std::cout << listofancestorids[i] << "; ";
   std::cout << std::endl;

   for( int b = 0; b < nBlocks; ++ b )
   {
      std::cout << getNConssForBlock( b ) << " constraint(s) in block " << b << std::endl;
      std::cout << getNVarsForBlock( b ) << " variable(s) in block " << b << std::endl;
      std::cout << getNStairlinkingvars( b ) << " stairlinkingvariable(s) in block " << b << std::endl;
   }

   std::cout << getNLinkingvars() << " linkingvariable(s)" << std::endl;
   std::cout << getNMasterconss() << " mastercontraint(s)" << std::endl;
   std::cout << getNMastervars() << " mastervariable(s)" << std::endl;
   std::cout << getNOpenconss() << " open constraint(s)" << std::endl;
   std::cout << getNOpenvars() << " open variable(s)" << std::endl;
   std::cout << "  stems from unpresolved problem: " << stemsFromUnpresolved << std::endl;
   std::cout << getNDetectors() << " detector(s)";
   if( getNDetectors() != 0 )
   {
      std::string detectorrepres;

      if( detectorChain[0] == NULL )
         detectorrepres = "user";
      else
         detectorrepres = (
            getNDetectors() != 1 || ! isFinishedByFinisher ? DECdetectorGetName( detectorChain[0] ) :
               "(finish)" + std::string( DECdetectorGetName( detectorChain[0] ) ) );

      std::cout << ": " << detectorrepres;

      for( int d = 1; d < getNDetectors(); ++ d )
      {
         detectorrepres = (
            getNDetectors() != d + 1 || ! isFinishedByFinisher ? DECdetectorGetName( detectorChain[d] ) :
               "(finish)" + std::string( DECdetectorGetName( detectorChain[d] ) ) );

         std::cout << ", " << detectorrepres;
      }
   }
   std::cout << "\n";

   return SCIP_OKAY;
}

/** displays the assignments of the vars */
SCIP_RETCODE Seeed::displayVars(
   )
{
   for( int b = 0; b < nBlocks; ++ b )
   {
      if( getNVarsForBlock( b ) != 0 )
      {
         std::cout << "variable(s) in block " << b << ": ";
         std::cout << getVarsForBlock( b )[0] << " ("
            << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( getVarsForBlock( b )[0] ) ) ) : "" )
            << ") ";
         for( int c = 1; c < getNVarsForBlock( b ); ++ c )
            std::cout << ", " << getVarsForBlock( b )[c] << " ("
               << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( getVarsForBlock( b )[c] ) ) ) : "" )
               << ") ";
         std::cout << "\n";
      }
      else
         std::cout << "0 variables in block " << b << std::endl;
      if( getNStairlinkingvars( b ) != 0 )
      {
         std::cout << "stairlinkingvariable(s) in block " << b << ": ";
         std::cout << getStairlinkingvars( b )[0] << " ("
            << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( getStairlinkingvars( b )[0] ) ) ) : "" )
            << ") ";
         for( int c = 1; c < getNStairlinkingvars( b ); ++ c )
            std::cout << ", " << getStairlinkingvars( b )[c] << " ("
               << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( getStairlinkingvars( b )[c] ) ) ) : "" )
               << ") ";
         std::cout << "\n";
      }
      else
         std::cout << "0 stairlinkingvariables in block " << b << std::endl;
   }

   if( getNLinkingvars() != 0 )
   {
      std::cout << "linkingvariable(s): ";
      std::cout << linkingVars[0] << " ("
         << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( linkingVars[0] ) ) ) : "" ) << ") ";
      for( int c = 1; c < getNLinkingvars(); ++ c )
         std::cout << ", " << linkingVars[c] << " ("
            << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( linkingVars[c] ) ) ) : "" ) << ") ";
      std::cout << "\n";
   }
   else
      std::cout << "0 linkingvariables" << std::endl;

   if( getNMastervars() != 0 )
   {
      std::cout << "mastervariable(s): ";
      std::cout << masterVars[0] << " ("
         << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( masterVars[0] ) ) ) : "" ) << ") ";
      for( int c = 1; c < getNMastervars(); ++ c )
         std::cout << ", " << masterVars[c] << " ("
            << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( masterVars[c] ) ) ) : "" ) << ") ";
      std::cout << "\n";
   }
   else
      std::cout << "0 mastervariables" << std::endl;

   if( getNOpenvars() != 0 )
   {
      std::cout << "open variable(s): ";
      std::cout << openVars[0] << " ("
         << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( openVars[0] ) ) ) : "" ) << ") ";
      for( int c = 1; c < getNOpenvars(); ++ c )
         std::cout << ", " << openVars[c] << " ("
            << ( seeedpool != NULL ? ( SCIPvarGetName( seeedpool->getVarForIndex( openVars[c] ) ) ) : "" ) << ") ";
      std::cout << "\n";
   }
   else
      std::cout << "0 open variables" << std::endl;

   return SCIP_OKAY;
}

/** computes the score of the given seeed based on the border, the average density score and the ratio of linking variables
 *  @todo bound calculation for unfinished decompositions could be more precise */
SCIP_Real Seeed::evaluate(
   SCORETYPE sctype
   )
{
   SCIP_Real borderscore; /**< score of the border */
   SCIP_Real densityscore; /**< score of block densities */
   SCIP_Real linkingscore; /**< score related to interlinking blocks */
   SCIP_Real totalscore; /**< accumulated score */

   int matrixarea;
   int borderarea;
   int i;
   int j;
   int k;
   /*   int blockarea; */
   SCIP_Real varratio;
   int* nzblocks;
   int* nlinkvarsblocks;
   int* nvarsblocks;
   SCIP_Real* blockdensities;
   int* blocksizes;
   SCIP_Real density;

   SCIP_Real alphaborderarea;
   SCIP_Real alphalinking;
   SCIP_Real alphadensity;

   unsigned long blackarea;

   maxwhitescore = 0.;
   alphaborderarea = 0.6;
   alphalinking = 0.2;
   alphadensity = 0.2;
   blackarea = 0;

   assert( checkConsistency() );

   /* calculate bound on max white score */
   if( getNOpenconss() != 0 || getNOpenvars() != 0 )
   {
      blackarea += ( getNLinkingvars() + getNTotalStairlinkingvars() ) * getNConss();
      blackarea += (unsigned long) getNMasterconss() * (unsigned long) getNVars();
      blackarea -= (unsigned long) getNMastervars() * (unsigned long) getNMasterconss();
      for( i = 0; i < nBlocks; ++ i )
      {
         blackarea += (unsigned long) getNConssForBlock( i ) * (unsigned long) getNVarsForBlock( i );
      }


      maxwhitescore = 1. - ( (SCIP_Real) blackarea / (SCIP_Real) ( (unsigned long) getNConss() * (unsigned long) getNVars() ) );

      return maxwhitescore;

   }

   if( getNOpenconss() != 0 || getNOpenvars() != 0 )
      SCIPwarningMessage( scip, "Evaluation for seeeds is not implemented for seeeds with open conss or open vars.\n" );

 //  if ( sctype == scoretype::MAX_FORESSEEING_WHITE || sctype == scoretype::SETPART_FWHITE )

   calcAggregationInformation();

   {
      std::vector<int> nlinkingvarsforblock(getNBlocks(), 0);
      std::vector<int> nblocksforlinkingvar(getNLinkingvars() + getNTotalStairlinkingvars(), 0);

      int sumblockshittinglinkingvar;
      int sumlinkingvarshittingblock;
      int newheight;
      int newwidth;
      int newmasterarea;
      int newblockarea;
      int newblockareaagg;


      for( int lv = 0; lv < getNLinkingvars(); ++lv )
      {
         int linkingvarid = getLinkingvars()[lv];

         for( int b = 0; b < getNBlocks(); ++b )
         {
            for ( int blc = 0; blc < getNConssForBlock(b); ++blc )
            {
               int blockcons = getConssForBlock(b)[blc];
               if( !SCIPisZero( seeedpool->getScip(), seeedpool->getVal(blockcons, linkingvarid) ) )
               {
                  /** linking var hits block */
                  ++nlinkingvarsforblock[b];
                  ++nblocksforlinkingvar[lv];
                  break;
               }
            }
         }
      }

      for( int b = 0; b < getNBlocks(); ++b)
      {
         for( int slv = 0; slv < getNStairlinkingvars(b); ++slv )
         {
            ++nlinkingvarsforblock[b];
            ++nlinkingvarsforblock[b+1];
            ++nblocksforlinkingvar[getNLinkingvars() + slv];
            ++nblocksforlinkingvar[getNLinkingvars() + slv];
         }
      }

      sumblockshittinglinkingvar = 0;
      sumlinkingvarshittingblock = 0;
      for( int b = 0; b < getNBlocks(); ++b )
      {
         sumlinkingvarshittingblock += nlinkingvarsforblock[b];
      }
      for( int lv = 0; lv < getNLinkingvars(); ++lv )
      {
         sumblockshittinglinkingvar += nblocksforlinkingvar[lv];
      }

      for( int slv = 0; slv < getNTotalStairlinkingvars(); ++slv )
      {
         sumblockshittinglinkingvar += nblocksforlinkingvar[getNLinkingvars() + slv];
      }


      newheight = getNConss() + sumblockshittinglinkingvar;
      newwidth = getNVars() + sumlinkingvarshittingblock;

      newmasterarea = ( getNMasterconss() + sumblockshittinglinkingvar) * ( getNVars() + sumlinkingvarshittingblock );
      newblockarea = 0;
      newblockareaagg = 0;

      for( int b = 0; b < getNBlocks(); ++b )
      {
         newblockarea += getNConssForBlock(b) * ( getNVarsForBlock(b) + nlinkingvarsforblock[b] );
      }

      for( int br = 0; br < nrepblocks; ++br )
      {
         newblockareaagg += getNConssForBlock( reptoblocks[br][0] ) * ( getNVarsForBlock( reptoblocks[br][0] ) + nlinkingvarsforblock[reptoblocks[br][0]] );
      }

      maxforeseeingwhitescore = ((SCIP_Real ) newblockarea + (SCIP_Real) newmasterarea) / (SCIP_Real) newwidth;
      maxforeseeingwhitescore =  maxforeseeingwhitescore / (SCIP_Real) newheight ;

      maxforeseeingwhitescoreagg = ((SCIP_Real ) newblockareaagg + (SCIP_Real) newmasterarea) / (SCIP_Real) newwidth;
      maxforeseeingwhitescoreagg =  maxforeseeingwhitescoreagg / (SCIP_Real) newheight ;

      maxforeseeingwhitescore = 1. - maxforeseeingwhitescore;
      maxforeseeingwhitescoreagg = 1. - maxforeseeingwhitescoreagg;

   }

   if( hasSetppccardMaster() && !isTrivial() && getNBlocks() > 1 )
   {
      setpartfwhitescore = 0.5 * maxforeseeingwhitescore + 0.5;
      setpartfwhitescoreagg = 0.5 * maxforeseeingwhitescoreagg + 0.5;
   }
   else
   {
      setpartfwhitescore = 0.5 * maxforeseeingwhitescore;
      setpartfwhitescoreagg = 0.5 * maxforeseeingwhitescoreagg;
   }

   //if( sctype == scoretype::SETPART_FWHITE )




   SCIP_CALL( SCIPallocBufferArray( scip, & nzblocks, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & nlinkvarsblocks, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & blockdensities, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & blocksizes, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & nvarsblocks, nBlocks ) );
   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate matrix area */
   matrixarea = nVars * nConss;

   blackarea += (unsigned long) ( getNLinkingvars() + getNTotalStairlinkingvars() ) * (unsigned long) getNConss();
   blackarea += (unsigned long) getNMasterconss() * ( (unsigned long) getNVars() - ( getNLinkingvars() + getNTotalStairlinkingvars() ) ) ;

   //std::cout << " black area without blocks is " <<  "(" << getNLinkingvars() << " + " << getNTotalStairlinkingvars() << " )  * " << getNConss() <<  " + " <<  getNMasterconss() << "  * ( " << getNVars() << "  -  ( " << getNLinkingvars() << " + " <<  getNTotalStairlinkingvars() << " ) ) "
   //   <<     " = " <<   blackarea << std::endl;

   if( sctype == SCORETYPE::MAX_WHITE)
   {
      for( i = 0; i < nBlocks; ++ i )
      {
         blackarea += (unsigned long) getNConssForBlock( i ) * ( (unsigned long) getNVarsForBlock( i ) );
      }
   }

   if( sctype != SCORETYPE::MAX_WHITE)
   {
      /* calculate slave sizes, nonzeros and linkingvars */
      for( i = 0; i < nBlocks; ++ i )
      {
         int ncurconss;
         int nvarsblock;
         SCIP_Bool *ishandled;

         SCIP_CALL( SCIPallocBufferArray( scip, & ishandled, nVars ) );
         nvarsblock = 0;
         nzblocks[i] = 0;
         nlinkvarsblocks[i] = 0;

         //    std::cout << "blackarea =  " << blackarea << " +  " << getNConssForBlock( i ) << " * " << getNVarsForBlock( i ) << " = " << getNConssForBlock( i ) * ( getNVarsForBlock( i ) );

         blackarea += (unsigned long) getNConssForBlock( i ) * ( (unsigned long) getNVarsForBlock( i ) );
         //  std::cout << " =  " << blackarea  << std::endl;

         for( j = 0; j < nVars; ++ j )
         {
            ishandled[j] = FALSE;
         }
         ncurconss = getNConssForBlock( i );

         for( j = 0; j < ncurconss; ++ j )
         {
            int cons = getConssForBlock( i )[j];
            int ncurvars;
            ncurvars = seeedpool->getNVarsForCons( cons );
            for( k = 0; k < ncurvars; ++ k )
            {
               int var = seeedpool->getVarsForCons( cons )[k];
               int block = -3;
               if( isVarBlockvarOfBlock( var, i ) )
                  block = i + 1;
               else if( isVarLinkingvar( var ) || isVarStairlinkingvar( var ) )
                  block = nBlocks + 2;
               else if( isVarMastervar( var ) )
                  block = nBlocks + 1;

               ++ ( nzblocks[i] );

               if( block == nBlocks + 1 && ishandled[var] == FALSE )
               {
                  ++ ( nlinkvarsblocks[i] );
               }
               ishandled[var] = TRUE;
            }
         }

         for( j = 0; j < nVars; ++ j )
         {
            if( ishandled[j] )
            {
               ++ nvarsblock;
            }
         }

         blocksizes[i] = nvarsblock * ncurconss;
         nvarsblocks[i] = nvarsblock;
         if( blocksizes[i] > 0 )
         {
            blockdensities[i] = 1.0 * nzblocks[i] / blocksizes[i];
         }
         else
         {
            blockdensities[i] = 0.0;
         }

         assert( blockdensities[i] >= 0 && blockdensities[i] <= 1.0 );
         SCIPfreeBufferArray( scip, & ishandled );
      }
   }

   borderarea = getNMasterconss() * nVars
      + ( getNLinkingvars() + getNMastervars() + getNTotalStairlinkingvars() ) * ( nConss - getNMasterconss() );

   maxwhitescore = 1. - ( (SCIP_Real) blackarea /  (SCIP_Real) ( (unsigned long) getNConss() * (unsigned long) getNVars() ) );
//   std::cout << "black area ration =  " << blackarea << "/ ( " << getNConss() << " * " << getNVars() << " =  " << ( (unsigned long) getNConss() * (unsigned long) getNVars() ) << ")  = " << maxwhitescore << std::endl;

   //std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    this seeed has a black area ratio of " << maxwhitescore << std::endl;

   density = 1E20;
   varratio = 1.0;
   linkingscore = 1.;
   borderscore =  1.;
   densityscore = 1.;

   if( sctype != SCORETYPE::MAX_WHITE )
   {
      for( i = 0; i < nBlocks; ++ i )
      {
         density = MIN( density, blockdensities[i] );

         if( ( getNLinkingvars() + getNMastervars() + getNTotalStairlinkingvars() ) > 0 )
         {
            varratio *= 1.0 * nlinkvarsblocks[i] / ( getNLinkingvars() + getNMastervars() + getNTotalStairlinkingvars() );
         }
         else
         {
            varratio = 0;
         }
      }
      linkingscore = ( 0.5 + 0.5 * varratio );

      densityscore = ( 1 - density );
   }

   borderscore = ( 1.0 * ( borderarea ) / matrixarea );
   borderareascore = 1. - borderscore;

   DEC_DECTYPE type;
   if( getNLinkingvars() == getNTotalStairlinkingvars() && getNMasterconss() == 0 && getNLinkingvars() > 0 )
   {
      type = DEC_DECTYPE_STAIRCASE;
   }
   else if( getNLinkingvars() > 0 || getNTotalStairlinkingvars() )
   {
      type = DEC_DECTYPE_ARROWHEAD;
   }
   else if( getNMasterconss() > 0 )
   {
      type = DEC_DECTYPE_BORDERED;
   }
   else if( getNMasterconss() == 0 && getNTotalStairlinkingvars() == 0 )
   {
      type = DEC_DECTYPE_DIAGONAL;
   }
   else
   {
      type = DEC_DECTYPE_UNKNOWN;
   }

   switch( type )
   {
      case DEC_DECTYPE_ARROWHEAD:
         totalscore = 1. - (alphaborderarea * ( borderscore ) + alphalinking * ( linkingscore ) + alphadensity * ( densityscore ) );
//      score->totalscore = score->borderscore*score->linkingscore*score->densityscore;
         break;
      case DEC_DECTYPE_BORDERED:
         totalscore = 1. - ( alphaborderarea * ( borderscore ) + alphalinking * ( linkingscore ) + alphadensity * ( densityscore ) );
//      score->totalscore = score->borderscore*score->linkingscore*score->densityscore;
         break;
      case DEC_DECTYPE_DIAGONAL:
         if( nBlocks == 1 || nBlocks == 0 )
            totalscore = 0.0;
         else
            totalscore = 1.0;
         break;
      case DEC_DECTYPE_STAIRCASE:
         totalscore = 1. - ( alphaborderarea * ( borderscore ) + alphalinking * ( linkingscore ) + 0.2 * ( densityscore ) );
         break;
      case DEC_DECTYPE_UNKNOWN:
         assert (FALSE);
         totalscore = 0.0;
         break;
      default:
         SCIPerrorMessage( "No rule for this decomposition type, cannot compute score\n" );
         assert( FALSE );
         totalscore = 0.0;
         break;
   }
   if( nBlocks == 0 )
      totalscore = 0.0;
   if( nBlocks == 1 )
      totalscore *= 0.25;
   if( totalscore > 1 )
      totalscore = 1;


   SCIPfreeBufferArray( scip, & nvarsblocks );
   SCIPfreeBufferArray( scip, & blocksizes );
   SCIPfreeBufferArray( scip, & blockdensities );
   SCIPfreeBufferArray( scip, & nlinkvarsblocks );
   SCIPfreeBufferArray( scip, & nzblocks );
   score = totalscore;

   return   getScore(sctype);

}


/**
 * returns true if the master consists only setpartitioning packing, covering, or cardinality constraints
 */
SCIP_Bool Seeed::hasSetppccardMaster(
)
{
   SCIP_Bool hassetpartmaster;
   SCIP_Bool verbose;
   hassetpartmaster = TRUE;
   verbose = FALSE;



   if( getNTotalStairlinkingvars() + getNLinkingvars() > 0 )
      return FALSE;

   for( int l = 0; l < getNMasterconss(); ++l )
   {
      int consid = getMasterconss()[l];
      if( !seeedpool->isConsSetppc(consid) && !seeedpool->isConsCardinalityCons(consid) )
      {
         hassetpartmaster = FALSE;
         if( verbose )
            std::cout <<   " cons with name  " << SCIPconsGetName( seeedpool->getConsForIndex(consid) ) << " is no setppccard constraint." << std::endl;
         break;
      }
   }

   return hassetpartmaster;
}


/**
 * returns true if the master consists only setpartitioning, packing, or covering constraints
 */
SCIP_Bool Seeed::hasSetppcMaster(
)
{
   SCIP_Bool hassetpartmaster;
   hassetpartmaster = TRUE;

    if( getNTotalStairlinkingvars() + getNLinkingvars() > 0 )
      return FALSE;


   for( int l = 0; l < getNMasterconss(); ++l )
   {
      int consid = getMasterconss()[l];
      if( !seeedpool->isConsSetppc(consid)  )
      {
         hassetpartmaster = FALSE;
         break;
      }
   }
   return hassetpartmaster;
}


/**
 * returns true if the master consists only setpartitioning, or packing constraints
 */
SCIP_Bool Seeed::hasSetppMaster(
)
{
   SCIP_Bool hassetpartmaster;
   hassetpartmaster = TRUE;

   if( getNTotalStairlinkingvars() + getNLinkingvars() > 0 )
      return FALSE;

   for( int l = 0; l < getNMasterconss(); ++l )
   {
      int consid = getMasterconss()[l];
      if( !seeedpool->isConsSetpp(consid)  )
      {
         hassetpartmaster = FALSE;
         break;
      }
   }
   return hassetpartmaster;
}




/** assigns all conss to master or declares them to be open (and declares all vars to be open)
 *  according to the cons assignment information given in constoblock hashmap
 *  precondition: no cons or var is already assigned to a block */
SCIP_RETCODE Seeed::filloutBorderFromConstoblock(
   SCIP_HASHMAP* constoblock,
   int givenNBlocks
   )
{
   assert( givenNBlocks >= 0 );
   assert( nBlocks == 0 );
   assert( (int) conssForBlocks.size() == nBlocks );
   assert( (int) varsForBlocks.size() == nBlocks );
   assert( (int) stairlinkingVars.size() == nBlocks );
   assert( ! alreadyAssignedConssToBlocks() );
   nBlocks = givenNBlocks;
   nVars = seeedpool->getNVars();
   nConss = seeedpool->getNConss();
   int consnum;
   int consblock;

   changedHashvalue = true;

   for( int i = 0; i < nConss; ++ i )
   {
      consnum = i;
      consblock = ( (int) (size_t) SCIPhashmapGetImage( constoblock, (void*) (size_t) i ) ) - 1;
      assert( consblock >= 0 && consblock <= nBlocks );
      if( consblock == nBlocks )
      {
         setConsToMaster( consnum );
         deleteOpencons( consnum );
      }
   }

   nBlocks = 0;
   sort();

   assert( checkConsistency( ) );

   return SCIP_OKAY;
}

/** assigns all conss to master or a block
 *  according to the cons assignment information given in constoblock hashmap
 *  calculates implicit variable assignment through cons assignment
 *  precondition: no cons or var is already assigned to a block and constoblock contains information for every cons */
SCIP_RETCODE Seeed::filloutSeeedFromConstoblock(
   SCIP_HASHMAP* constoblock,
   int givenNBlocks
   )
{
   assert( givenNBlocks >= 0 );
   assert( nBlocks == 0 );
   assert( (int) conssForBlocks.size() == nBlocks );
   assert( (int) varsForBlocks.size() == nBlocks );
   assert( (int) stairlinkingVars.size() == nBlocks );
   assert( ! alreadyAssignedConssToBlocks() );
   nBlocks = givenNBlocks;
   nVars = seeedpool->getNVars();
   nConss = seeedpool->getNConss();
   int consnum;
   int consblock;
   int varnum;
   bool varInBlock;
   std::vector<int> varInBlocks = std::vector<int>( 0 );
   std::vector<int> emptyVector = std::vector<int>( 0 );

   changedHashvalue = true;

   for( int c = 0; c < nConss; ++ c )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "%d\n", c);
      assert( SCIPhashmapExists( constoblock, (void*) (size_t) c ) );
      assert( (int) (size_t) SCIPhashmapGetImage( constoblock, (void*) (size_t) c ) - 1 <= nBlocks );
      assert( (int) (size_t) SCIPhashmapGetImage( constoblock, (void*) (size_t) c ) - 1 >= 0 );
   }

   for( int b = (int) conssForBlocks.size(); b < nBlocks; b ++ )
      conssForBlocks.push_back( emptyVector );

   for( int b = (int) varsForBlocks.size(); b < nBlocks; b ++ )
      varsForBlocks.push_back( emptyVector );

   for( int b = (int) stairlinkingVars.size(); b < nBlocks; b ++ )
      stairlinkingVars.push_back( emptyVector );

   for( int i = 0; i < nConss; ++ i )
   {
      consnum = i;
      consblock = ( (int) (size_t) SCIPhashmapGetImage( constoblock, (void*) (size_t) i ) ) - 1;
      assert( consblock >= 0 && consblock <= nBlocks );
      if( consblock == nBlocks )
         setConsToMaster( consnum );
      else
         setConsToBlock( consnum, consblock );
   }

   for( int i = 0; i < nVars; ++ i )
   {
      varInBlocks.clear();
      varnum = i;

      /** test if the variable can be found in blocks */
      for( int b = 0; b < nBlocks; ++ b )
      {
         varInBlock = false;
         for( size_t k = 0; k < conssForBlocks[b].size() && ! varInBlock; ++ k )
         {
            for( int l = 0; l < seeedpool->getNVarsForCons( conssForBlocks[b][k] ) && ! varInBlock; ++ l )
            {
               if( varnum == ( seeedpool->getVarsForCons( conssForBlocks[b][k] ) )[l] )
               {
                  varInBlocks.push_back( b );
                  varInBlock = true;
               }
            }
         }
      }
      if( varInBlocks.size() == 1 ) /** if the var can be found in one block set the var to block var */
         setVarToBlock( varnum, varInBlocks[0] );
      else if( varInBlocks.size() == 2 ) /** if the variable can be found in two blocks check if it is a linking var or a stairlinking var*/
      {
         if( varInBlocks[0] + 1 == varInBlocks[1] )
            setVarToStairlinking( varnum, varInBlocks[0], varInBlocks[1] );
         else
            setVarToLinking( varnum );
      }
      else if( varInBlocks.size() > 2 ) /** if the variable can be found in more than two blocks it is a linking var */
         setVarToLinking( varnum );
      else
         assert( varInBlocks.size() == 0 );
      setVarToMaster( varnum );
   }
   sort();
   openVars = std::vector<int>( 0 );
   openConss = std::vector<int>( 0 );
   isvaropen = std::vector<bool>(nVars, false);
   isconsopen = std::vector<bool>(nConss, false);

   deleteEmptyBlocks(false);
   sort();
   assert( checkConsistency( ) );

   return SCIP_OKAY;
}

/** reassigns variables classified as linking to master if the variable only hits master conss */
SCIP_RETCODE Seeed::findVarsLinkingToMaster(
   )
{
   int i;
   int j;
   const int* varcons;
   bool isMasterVar;
   const int* lvars = getLinkingvars();
   std::vector<int> foundMasterVarIndices;

   changedHashvalue = true;

   // sort Master constraints for binary search
   sort();

   for( i = 0; i < getNLinkingvars(); ++ i )
   {
      isMasterVar = true;
      varcons = seeedpool->getConssForVar( lvars[i] );
      for( j = 0; j < seeedpool->getNConssForVar( lvars[i] ); ++ j )
      {
         if( ! isconsmaster[varcons[j]]  )
         {
            isMasterVar = false;
            break;
         }
      }

      if( isMasterVar )
      {
         foundMasterVarIndices.push_back( i );
      }
   }

   for( std::vector<int>::reverse_iterator it = foundMasterVarIndices.rbegin(); it != foundMasterVarIndices.rend(); ++ it )
   {
      masterVars.push_back( lvars[ * it] );
      mastervarssorted = false;
      isvarmaster[lvars[ * it]] = true;
      linkingVars.erase( linkingVars.begin() + * it );
   }

   return SCIP_OKAY;
}

/** reassigns variables classified as linking to stairlinking if the variable hits conss in exactly two consecutive blocks */
SCIP_RETCODE Seeed::findVarsLinkingToStairlinking(
   )
{
   int i;
   int j;
   int k;

   int consblock;
   int block1 = - 1;
   int block2 = - 1;

   const int* varcons;
   const int* lvars = getLinkingvars();

   std::vector<int> foundMasterVarIndices;

   sort();

   for( i = 0; i < getNLinkingvars(); ++ i )
   {
      block1 = - 1;
      block2 = - 1;
      varcons = seeedpool->getConssForVar( lvars[i] );
      for( j = 0; j < seeedpool->getNConssForVar( lvars[i] ); ++ j )
      {
         consblock = - 1;
         for( k = 0; k < nBlocks; ++ k )
         {
            if( std::binary_search( conssForBlocks[k].begin(), conssForBlocks[k].end(), varcons[j] ) )
            {
               consblock = k;
               break;
            }
         }

         if( consblock == - 1 )
         {
            block1 = - 1;
            block2 = - 1;
            break;
         }
         else if( block1 == consblock || block2 == consblock )
         {
            continue;
         }
         else if( block1 == - 1 )
         {
            block1 = consblock;
            continue;
         }
         else if( block2 == - 1 )
         {
            block2 = consblock;
            continue;
         }
         else
         {
            block1 = - 1;
            block2 = - 1;
            break;
         }
      }

      if( block1 != - 1 && block2 != - 1 && ( block1 == block2 + 1 || block1 + 1 == block2 ) )
      {
//    	 std::cout << "Var " << lvars[i] << " hits block " << block1 << " and " << block2 << "\n";

         setVarToStairlinking( lvars[i], block1, block2 );
         foundMasterVarIndices.push_back( i );
      }
   }

   for( std::vector<int>::reverse_iterator it = foundMasterVarIndices.rbegin(); it != foundMasterVarIndices.rend(); ++ it )
   {
      linkingVars.erase( linkingVars.begin() + * it );
   }

   return SCIP_OKAY;
}

/** returns a vector of pairs of var indices and vectors of (two) block indices
 *  the related linking variable hits exactly the two blocks given in the related vector */
std::vector< std::pair< int, std::vector< int > > > Seeed::findLinkingVarsPotentiallyStairlinking(
   )
{
	std::vector< std::pair< int, std::vector< int > > > blocksOfVars( 0 );
	const int* varcons;
	const int* lvars;
	int blockcounter;

   /* if there is at least one linking variable, then the blocks of vars must be created. */
   if( getNLinkingvars() > 0 )
   {
      lvars = getLinkingvars();
      sort();

      /* check every linking var */
      for ( int v = 0; v < getNLinkingvars(); ++v )
      {
         std::vector< int > blocksOfVar( 0 );
         blockcounter = 0;

         varcons = seeedpool->getConssForVar( lvars[v] );

         /* find all blocks that are hit by this linking var */
         for ( int c = 0; c < seeedpool->getNConssForVar( lvars[v] ) && blockcounter <= 2; ++c )
         {
            for ( int b = 0; b < nBlocks && blockcounter <= 2; ++b )
            {
               if ( std::binary_search( conssForBlocks[b].begin(),
                     conssForBlocks[b].end(), varcons[c] ) )
               {
                  /* if the hit block is new, add it to blockOfVar vector */
                  if ( std::find( blocksOfVar.begin(), blocksOfVar.end(), b ) == blocksOfVar.end() )
                  {
                     //std::cout << "Var " << lvars[v] << " hits block " << b << "\n" ;
                     ++blockcounter;
                     blocksOfVar.push_back( b );
                  }
               }
            }
         }

         /* if the var hits exactly two blocks, it is potentially stairlinking */
         if ( blockcounter == 2 )
         {
            std::pair< int, std::vector< int > > pair( v, blocksOfVar );
            blocksOfVars.push_back( pair );
         }
      }
   }

	return blocksOfVars;
}

/** assigns all booked constraints and variables and deletes them from list of open cons and open vars */
SCIP_RETCODE Seeed::flushBooked()
{
   std::vector<int>::const_iterator bookedIter;
   std::vector<int>::const_iterator bookedIterEnd;
   std::vector<std::pair<int, int>>::iterator bookedIter2;
   std::vector<std::pair<int, int>>::iterator bookedIterEnd2;

   std::vector < SCIP_Bool > varislinking( getNVars(), FALSE );

   changedHashvalue = true;

   bookedIter = bookedAsMasterConss.begin();
   bookedIterEnd = bookedAsMasterConss.end();
   for( ; bookedIter != bookedIterEnd; ++ bookedIter )
   {
      setConsToMaster( * bookedIter );
      deleteOpencons( * bookedIter );
   }
   bookedAsMasterConss.clear();

   bookedIter2 = bookedAsBlockConss.begin();
   bookedIterEnd2 = bookedAsBlockConss.end();
   for( ; bookedIter2 != bookedIterEnd2; ++ bookedIter2 )
   {
      setConsToBlock( ( * bookedIter2 ).first, ( * bookedIter2 ).second );
      deleteOpencons( ( * bookedIter2 ).first );
   }
   bookedAsBlockConss.clear();

   bookedIter = bookedAsLinkingVars.begin();
   bookedIterEnd = bookedAsLinkingVars.end();
   for( ; bookedIter != bookedIterEnd; ++ bookedIter )
   {
      varislinking[ * bookedIter] = TRUE;
      setVarToLinking( * bookedIter );
      deleteOpenvar( * bookedIter );
   }
   bookedAsLinkingVars.clear();

   bookedIter = bookedAsMasterVars.begin();
   bookedIterEnd = bookedAsMasterVars.end();

   for( ; bookedIter != bookedIterEnd; ++ bookedIter )
   {
      setVarToMaster( * bookedIter );
      deleteOpenvar( * bookedIter );
   }
   bookedAsMasterVars.clear();

   bookedIter2 = bookedAsBlockVars.begin();
   bookedIterEnd2 = bookedAsBlockVars.end();
   for( ; bookedIter2 != bookedIterEnd2; ++ bookedIter2 )
   {
      if( varislinking[( * bookedIter2 ).first] )
         continue;
      setVarToBlock( ( * bookedIter2 ).first, ( * bookedIter2 ).second );
      deleteOpenvar( ( * bookedIter2 ).first );
   }
   bookedAsBlockVars.clear();

   bookedIter2 = bookedAsStairlinkingVars.begin();
   bookedIterEnd2 = bookedAsStairlinkingVars.end();
   for( ; bookedIter2 != bookedIterEnd2; ++ bookedIter2 )
   {
      setVarToStairlinking( ( * bookedIter2 ).first, ( * bookedIter2 ).second, ( * bookedIter2 ).second + 1 );
      deleteOpenvar( ( * bookedIter2 ).first );
   }
   bookedAsStairlinkingVars.clear();

   sort();

   return SCIP_OKAY;
}

/** returns ancestor id of given ancestor */
int Seeed::getAncestorID(
   int ancestorindex
   )
{
   assert( 0 <= ancestorindex && ancestorindex < (int) listofancestorids.size() );

   return listofancestorids[ancestorindex];
}


/** returns ancestor id of given ancestor */
std::vector<int> Seeed::getAncestorList(
   )
{
   return listofancestorids;
}

const std::vector<int> & Seeed::getBlocksForRep(int repid)
{
   return reptoblocks[repid];
}


/** returns ancestor id of given ancestor */
void Seeed::setAncestorList(
   std::vector<int> newlist
   )
{
   listofancestorids = newlist ;
}



/** returns ancestor id of given ancestor */
void Seeed::addAncestorID(
   int ancestor
   )
{
   assert( 0 <= ancestor );
   listofancestorids.push_back(ancestor);
}


/** returns detectorchainstring */
char* Seeed::getDetectorChainString()
{
   return detectorchainstring;
}

/** returns detectorchain info of detector related to given detectorchain index */
std::string Seeed::getDetectorchainInfo(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorchaininfo.size() );

   return detectorchaininfo[detectorchainindex];
}

/** returns the time that the detector related to the given detectorchainindex needed for detecting */
SCIP_Real Seeed::getDetectorClockTime(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorClockTimes.size() );

   return detectorClockTimes[ detectorchainindex ];
}

/** returns the time that the detectors needed for detecting */
std::vector<SCIP_Real> Seeed::getDetectorClockTimes()
{
   return detectorClockTimes;
}


std::string Seeed::getComponentInformation(
   ){

   std::stringstream buf;

   int minrow = getNConss();
   int maxrow = 0;
   SCIP_Real medianrow = 0.;
   SCIP_Real meanrow = 0.;
   std::vector<int> nrows(getNBlocks(), 0);

   int ntotalcols = getNVars();
   int ntotalrows = seeedpool->getNTotalConss();

   int minvar = getNVars();
   int maxvar = 0;
   SCIP_Real medianvar = 0;
   SCIP_Real meanvar = 0.;
   std::vector<int> ncols(getNBlocks(), 0);

   for( int b = 0; b  < getNBlocks(); ++b )
   {
       for( int c = 0; c < getNConssForBlock(b); ++c  )
       {
          int cons = getConssForBlock(b)[c];
          SCIP_Cons* scipcons = seeedpool->getScipCons(cons);

          nrows[b] += GCGconsIsRanged(scip, scipcons) ? 2 : 1;
          meanrow += GCGconsIsRanged(scip, scipcons) ? 2 : 1;
       }

       for( int v = 0; v < getNVarsForBlock(b); ++v )
       {
          ncols[b]++;
          meanvar += 1;
       }
   }

   std::sort(nrows.begin(), nrows.end());
   std::sort(ncols.begin(), ncols.end());

   minrow = nrows[0];
   minvar = ncols[0];

   maxrow = nrows[getNBlocks()-1];
   maxvar = ncols[getNBlocks()-1];

   meanrow = meanrow / getNBlocks();
   meanvar = meanvar / getNBlocks();

   if( getNBlocks() % 2 == 0)
   {
      medianrow =  ( nrows[ getNBlocks() / 2  ] + nrows[ getNBlocks() / 2  - 1 ] ) / 2.;
       medianvar = ( ncols[ getNBlocks() / 2  ] + ncols[ getNBlocks() / 2  - 1 ] ) / 2.;
   }
   else
   {
      medianrow = nrows[ getNBlocks() / 2  ];
      medianvar = ncols[ getNBlocks() / 2  ];
   }

   buf << getNBlocks() << ", " << ( (SCIP_Real) minrow ) / ntotalrows << ", " << ( (SCIP_Real) maxrow ) / ntotalrows << ", ";
   buf << ( (SCIP_Real) medianrow ) / ntotalrows << ", " << ( (SCIP_Real) meanrow ) / ntotalrows << ", ";

   buf << ( (SCIP_Real) minvar ) / ntotalcols << ", " << ( (SCIP_Real) maxvar ) / ntotalcols << ", ";
   buf << ( (SCIP_Real) medianvar ) / ntotalcols << ", " << ( (SCIP_Real) meanvar ) / ntotalcols << " ";


   return buf.str();
}

/** returns the data of the consclassifier that the given detector made use of */
SCIP_RETCODE Seeed::getConsClassifierData(
   int detectorchainindex,
   ConsClassifier** classifier,
   std::vector<int>& consclassesmaster
   )
{
   assert( consClassifierUsed( detectorchainindex ) );

   *classifier = dynamic_cast<ConsClassifier*>( usedClassifier[detectorchainindex] );
   consclassesmaster = classesToMaster[detectorchainindex];

   return SCIP_OKAY;
}

/** returns the time that the detectors needed for detecting */
void Seeed::setDetectorClockTimes(
   std::vector<SCIP_Real> newvector)
{
   detectorClockTimes = newvector;
}

/** returns true if the given detector used a varclassifier */
bool Seeed::varClassifierUsed(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) usedClassifier.size() );

   return ( usedClassifier[detectorchainindex] != NULL )
      && ( dynamic_cast<VarClassifier*>( usedClassifier[detectorchainindex] ) != NULL );
}


/** returns array containing constraints assigned to a block */
const int* Seeed::getConssForBlock(
   int block
   )
{
   assert( block >= 0 && block < nBlocks );
   return & conssForBlocks[block][0];
}

/** returns the detectorchain */
DEC_DETECTOR** Seeed::getDetectorchain()
{
   return & detectorChain[0];
}

/** returns the detectorchain as a vector */
std::vector<DEC_DETECTOR*> Seeed::getDetectorchainVector()
{
   return detectorChain;
}

/** returns a string displaying detector-related information, i.e. clock times and assignment data */
std::string Seeed::getDetectorStatistics(
   int detectorchainindex
   )
{
   std::stringstream output;

   if( (int) getDetectorClockTimes().size() > detectorchainindex )
      output << "  Detection time: " << getDetectorClockTime( detectorchainindex ) << std::endl;
   if( (int) getPctConssFromFreeVector().size() > detectorchainindex )
      output << "  % newly assigned constraints: " << getPctConssFromFree( detectorchainindex ) << std::endl;
   if( (int) getPctConssToBorderVector().size() > detectorchainindex )
      output << "  % constraints the detector assigned to border: " << getPctConssToBorder( detectorchainindex ) << std::endl;
   if( (int) getPctConssToBlockVector().size() > detectorchainindex )
      output << "  % constraints the detector assigned to blocks: " << getPctConssToBlock( detectorchainindex ) << std::endl;
   if( (int) getPctVarsFromFreeVector().size() > detectorchainindex )
      output << "  % newly assigned variables: " << getPctVarsFromFree( detectorchainindex ) << std::endl;
   if( (int) getPctVarsToBorderVector().size() > detectorchainindex )
      output << "  % variables the detector assigned to border: " << getPctVarsToBorder( detectorchainindex ) << std::endl;
   if( (int) getPctVarsToBlockVector().size() > detectorchainindex )
      output << "  % variables the detector assigned to blocks: " << getPctVarsToBlock( detectorchainindex ) << std::endl;
   if( (int) getNNewBlocksVector().size() > detectorchainindex )
         output << "  New blocks: " << getNNewBlocks( detectorchainindex ) << std::endl;

   return output.str();
}

/** returns a string displaying classifier information if such a classifier was used */
std::string Seeed::getDetectorClassifierInfo(
   int detectorchainindex,
   bool displayConssVars
   )
{
   std::stringstream output;

   if( consClassifierUsed( detectorchainindex ) )
   {
      ConsClassifier* classifier;
      std::vector<int> constomaster;

      getConsClassifierData( detectorchainindex, &classifier, constomaster );

      output << "  Used consclassifier: " << classifier->getName() << std::endl;
      output << "   Pushed to master:";

      if( constomaster.size() > 0 )
      {
         if( displayConssVars )
         {
            output << std::endl << "    " << classifier->getClassName( constomaster[0] ) << " ("
               << classifier->getClassDescription( constomaster[0] ) << "): ";
            bool first = true;
            for( int c = 0; c < classifier->getNConss(); ++c )
            {
               if( classifier->getClassOfCons( c ) == constomaster[0] )
               {
                  if( first )
                  {
                     output << SCIPconsGetName( seeedpool->getConsForIndex( c ) );
                     first = false;
                  }
                  else
                  {
                     output << ", " << SCIPconsGetName( seeedpool->getConsForIndex( c ) );
                  }
               }
            }
            output << std::endl;
         }
         else
         {
            output << " " << classifier->getClassName( constomaster[0] );
         }
      }

      for( size_t i = 1; i < constomaster.size(); ++i )
      {
         if( displayConssVars )
         {
            output << "    " << classifier->getClassName( constomaster[i] ) << " ("
               << classifier->getClassDescription( constomaster[i] ) << "): ";
            bool first = true;
            for( int c = 0; c < classifier->getNConss(); ++c )
            {
               if( classifier->getClassOfCons( c ) == constomaster[i] )
               {
                  if( first )
                  {
                     output << SCIPconsGetName( seeedpool->getConsForIndex( c ) );
                     first = false;
                  }
                  else
                  {
                     output << ", " << SCIPconsGetName( seeedpool->getConsForIndex( c ) );
                  }
               }
            }
            output << std::endl;
         }
         else
         {
            output << ", " << classifier->getClassName( constomaster[i] );
         }
      }

      if ( !displayConssVars || constomaster.size() == 0 )
      {
         output << std::endl;
      }
   }

   if( varClassifierUsed( detectorchainindex ) )
   {
      VarClassifier* classifier;
      std::vector<int> vartolinking;
      std::vector<int> vartomaster;

      getVarClassifierData( detectorchainindex, &classifier, vartolinking, vartomaster );

      output << "  Used varclassifier: " << classifier->getName() << std::endl;
      output << "   Pushed to linking:";

      if( vartolinking.size() > 0 )
      {
         if( displayConssVars )
         {
            output << std::endl << "    " << classifier->getClassName( vartolinking[0] ) << " ("
               << classifier->getClassDescription( vartolinking[0] ) << "): ";
            bool first = true;
            for( int v = 0; v < classifier->getNVars(); ++v )
            {
               if( classifier->getClassOfVar( v ) == vartolinking[0] )
               {
                  if( first )
                  {
                     output << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                     first = false;
                  }
                  else
                  {
                     output << ", " << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                  }
               }
            }
            output << std::endl;
         }
         else
         {
            output << " " << classifier->getClassName( vartolinking[0] );
         }
      }

      for( size_t i = 1; i < vartolinking.size(); ++i )
      {
         if( displayConssVars )
         {
            output << "    " << classifier->getClassName( vartolinking[i] ) << " ("
               << classifier->getClassDescription( vartolinking[i] ) << "): ";
            bool first = true;
            for( int v = 0; v < classifier->getNVars(); ++v )
            {
               if( classifier->getClassOfVar( v ) == vartolinking[i] )
               {
                  if( first )
                  {
                     output << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                     first = false;
                  }
                  else
                  {
                     output << ", " << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                  }
               }
            }
            output << std::endl;
         }
         else
         {
            output << ", " << classifier->getClassName( vartolinking[i] );
         }
      }

      if ( !displayConssVars || vartolinking.size() == 0 )
      {
         output << std::endl;
      }

      output << "   Pushed to master:";

      if( vartomaster.size() > 0 )
      {
         if( displayConssVars )
         {
            output << std::endl << "    " << classifier->getClassName( vartomaster[0] ) << " ("
               << classifier->getClassDescription( vartomaster[0] ) << "): ";
            bool first = true;
            for( int v = 0; v < classifier->getNVars(); ++v )
            {
               if( classifier->getClassOfVar( v ) == vartomaster[0] )
               {
                  if( first )
                  {
                     output << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                     first = false;
                  }
                  else
                  {
                     output << ", " << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                  }
               }
            }
            output << std::endl;
         }
         else
         {
            output << " " << classifier->getClassName( vartomaster[0] );
         }
      }

      for( size_t i = 1; i < vartomaster.size(); ++i )
      {
         if( displayConssVars )
         {
            output << "    " << classifier->getClassName( vartomaster[i] ) << " ("
               << classifier->getClassDescription( vartomaster[i] ) << "): ";
            bool first = true;
            for( int v = 0; v < classifier->getNVars(); ++v )
            {
               if( classifier->getClassOfVar( v ) == vartolinking[i] )
               {
                  if( first )
                  {
                     output << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                     first = false;
                  }
                  else
                  {
                     output << ", " << SCIPvarGetName( seeedpool->getVarForIndex( v ) );
                  }
               }
            }
            output << std::endl;
         }
         else
         {
            output << ", " << classifier->getClassName( vartomaster[i] );
         }
      }

      if ( !displayConssVars || vartomaster.size() == 0 )
      {
         output << std::endl;
      }
   }

   return output.str();
}

/** returns true if this seeed was finished by finishSeeed() method of a detector */
bool Seeed::getFinishedByFinisher()
{
   return isFinishedByFinisher;
}

/** returns true if the seeed is finished by a finisher in the unpresolved problem */
bool Seeed::getFinishedByFinisherUnpresolved()
{
   return isFinishedByFinisherUnpresolved;
}

/** returns the detector that finished this seeed in the unpresolved problem if there exists one, NULL otherwise */
DEC_DETECTOR* Seeed::getFinishedUnpresolvedBy()
{
   return finishedUnpresolvedBy;
}

/** returns the calculated hash value of this seeed */
long Seeed::getHashValue()
{
   if( changedHashvalue )
      calcHashvalue();
   changedHashvalue = false;
   return this->hashvalue;
}

/** returns the id of the seeed */
int Seeed::getID()
{
   return id;
}

/** returns array containing all linking vars */
const int* Seeed::getLinkingvars()
{
   return & linkingVars[0];
}

/** returns array containing all master conss */
const int* Seeed::getMasterconss()
{
   return & masterConss[0];
}

/** returns array containing all master vars (every constraint containing a master var is in master) */
const int* Seeed::getMastervars()
{
   return & masterVars[0];
}

/** returns the "maximum white score" */
SCIP_Real Seeed::getMaxWhiteScore()
{

   return getScore(SCORETYPE::MAX_WHITE);
}


/** returns the "maximum white score" with adaptions for benders*/
SCIP_Real Seeed::getBendersScore()
{

   return getScore(SCORETYPE::BENDERS);
}


/** returns the number of nonzero coeffs in a certain block */
int  Seeed::getNCoeffsForBlock(
   int blockid
   ){

   if( !calculatedncoeffsforblock )
      calcNCoeffsForBlocks();

   return ncoeffsforblock[blockid];
}


/** returns the number of nonzero coeffs in master */
int  Seeed::getNCoeffsForMaster(
   ){

   if( !calculatedncoeffsforblock )
      calcNCoeffsForBlocks();

   return ncoeffsformaster;
}


/** returns the score of the seeed (depending on used scoretype) */
SCIP_Real Seeed::getScore(
   SCORETYPE type
   )
{
   /** if there are indicator constraints in the master we want to reject this decomposition */
   for( int mc = 0; mc < getNMasterconss(); ++mc )
   {
      SCIP_CONS* cons;
      cons = getSeeedpool()->getScipCons(getMasterconss()[mc]);
      if( GCGconsGetType(cons) == consType::indicator )
         return 0.;
   }

   /** calculate maximum white score anyway */
   if( maxwhitescore == -1. )
      calcmaxwhitescore();

   if( type == scoretype::MAX_WHITE )
      return maxwhitescore;

   if( type == scoretype::CLASSIC )
   {
      if ( score == -1. )
         SCIP_CALL_ABORT(calcclassicscore() );
      return score;
   }

   if( type == scoretype::BORDER_AREA )
   {
      if( borderareascore == -1. )
         calcborderareascore();
      return borderareascore;
   }

   if( type == scoretype::MAX_FORESSEEING_WHITE )
   {
      if( maxforeseeingwhitescore == -1. )
         calcmaxforeseeingwhitescore();
      return maxforeseeingwhitescore;
   }

   if( type == scoretype::MAX_FORESEEING_AGG_WHITE )
   {
      if( maxforeseeingwhitescoreagg == -1. )
         calcmaxforeseeingwhitescoreagg();
      return maxforeseeingwhitescoreagg;
   }

   if( type == scoretype::SETPART_FWHITE )
   {
      if( setpartfwhitescore == -1. )
         calcsetpartfwhitescore();
      return setpartfwhitescore;
   }

   if( type == scoretype::SETPART_AGG_FWHITE )
   {
      if( setpartfwhitescoreagg == -1. )
         calcsetpartfwhitescoreagg();
      return setpartfwhitescoreagg;
   }

   if( type == scoretype::BENDERS )
      {
         if( bendersscore == -1. )
            calcbendersscore();
         return bendersscore;
      }


   return 0;
}

/** returns whether this seeed is usergiven */
USERGIVEN Seeed::getUsergiven()
{
   return usergiven;
}

/** returns number of ancestor seeeds */
int Seeed::getNAncestors()
{
   return listofancestorids.size();
}

/** returns number of blocks */
int Seeed::getNBlocks()
{
   return nBlocks;
}

/** returns number of conss */
int Seeed::getNConss()
{
   return nConss;
}

/** returns size of the vector containing conss assigned to a block */
int Seeed::getNConssForBlock(
   int block
   )
{
   assert( block >= 0 && block < nBlocks );
   return (int) conssForBlocks[block].size();
}

/** returns size of the detectorchain info vector */
int Seeed::getNDetectorchainInfo()
{
   return detectorchaininfo.size();
}

/** returns the number of detectors the seeed is propagated by */
int Seeed::getNDetectors()
{
   if ( usergiven == USERGIVEN::NOT )
      return (int) detectorChain.size();
   else
      return 0;
}

/** returns the number used classifiers */
int Seeed::getNUsedClassifier()
{
   return (int) usedClassifier.size();
}


/** returns size of the vector containing linking vars */
int Seeed::getNLinkingvars()
{
   return (int) linkingVars.size();
}

/** returns number of blocks a detector added */
int Seeed::getNNewBlocks(
      int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorChain.size() );

   return nNewBlocks[detectorchainindex];
}

/** returns number of blocks the detectors in the detectorchain added */
std::vector<int> Seeed::getNNewBlocksVector()
{
   return nNewBlocks;
}

/** returns number of blocks the detectors in the detectorchain added */
void Seeed::setNNewBlocksVector(
   std::vector<int>  newvector
   )
{
   nNewBlocks = newvector;
}


/** returns size of the vector containing master conss */
int Seeed::getNMasterconss()
{
   return (int) masterConss.size();
}

/** returns size of the vector containing master vars (hitting only constraints in the master) */
int Seeed::getNMastervars()
{
   return (int) masterVars.size();
}


/** returns total number of stairlinking vars */
int Seeed::getNTotalStairlinkingvars()
{
   int nstairlinkingvars = 0;
   for( int b = 0; b < getNBlocks(); ++ b )
      nstairlinkingvars += getNStairlinkingvars( b );

   return nstairlinkingvars;
}

/** returns size of vector containing constraints not assigned yet */
int Seeed::getNOpenconss()
{
   return (int) openConss.size();
}

/** returns size of vector containing variables not assigned yet */
int Seeed::getNOpenvars()
{
   return (int) openVars.size();
}

/** returns the number of blockrepresentatives */
int Seeed::getNReps(){

   return nrepblocks;
}

/** returns size of the vector containing stairlinking vars */
int Seeed::getNStairlinkingvars(
   int block
   )
{
   assert( block >= 0 && block < nBlocks );
   return (int) stairlinkingVars[block].size();
}

/** returns number of vars */
int Seeed::getNVars()
{
   return nVars;
}

/** returns size of the vector containing vars assigned to a block */
int Seeed::getNVarsForBlock(
   int block
   )
{
   assert( block >= 0 && block < nBlocks );
   return (int) varsForBlocks[block].size();
}

/** returns array containing constraints not assigned yet */
const int* Seeed::getOpenconss()
{
   return & openConss[0];
}

/** returns array containing constraints not assigned yet*/
std::vector<int> Seeed::getOpenconssVec()
{
   return openConss;
}


/** returns array containing variables not assigned yet*/
const int* Seeed::getOpenvars()
{
   return & openVars[0];
}

/** returns array containing variables not assigned yet*/
std::vector<int> Seeed::getOpenvarsVec()
{
   return openVars;
}

/** returns fraction of variables assigned to the border for a detector */
SCIP_Real Seeed::getPctVarsToBorder(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorChain.size() );

   return pctVarsToBorder[detectorchainindex];
}

/** returns fraction of variables assigned to the border for detectors in detectorchain */
std::vector<SCIP_Real> Seeed::getPctVarsToBorderVector()
{
   return pctVarsToBorder;
}

/** returns fraction of variables assigned to the border for detectors in detectorchain */
void Seeed::setPctVarsToBorderVector(
   std::vector<SCIP_Real> newvector
   )
{
   pctVarsToBorder = newvector;
}



/** returns fraction of variables assigned to a block for a detector */
SCIP_Real Seeed::getPctVarsToBlock(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorChain.size() );

   return pctVarsToBlock[detectorchainindex];
}

/** returns fraction of variables assigned to a block for detectors in detectorchain */
std::vector<SCIP_Real> Seeed::getPctVarsToBlockVector()
{
   return pctVarsToBlock;
}


/** returns fraction of variables assigned to a block for detectors in detectorchain */
void Seeed::setPctVarsToBlockVector(
   std::vector<SCIP_Real> newvector
)
{
   pctVarsToBlock = newvector;
}



/** returns fraction of variables that are not longer open for a detector */
SCIP_Real Seeed::getPctVarsFromFree(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorChain.size() );

   return pctVarsFromFree[detectorchainindex];
}

/** returns fraction of variables that are not longer open for detectors in detectorchain */
std::vector<SCIP_Real> Seeed::getPctVarsFromFreeVector()
{
   return pctVarsFromFree;
}

/** returns fraction of variables that are not longer open for detectors in detectorchain */
void Seeed::setPctVarsFromFreeVector(
   std::vector<SCIP_Real> newvector
   )
{
   pctVarsFromFree = newvector;
}


/** returns fraction of constraints assigned to the border for a detector */
SCIP_Real Seeed::getPctConssToBorder(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorChain.size() );

   return pctConssToBorder[detectorchainindex];
}

/** returns fraction of constraints assigned to the border for detectors in detectorchain */
std::vector<SCIP_Real> Seeed::getPctConssToBorderVector()
{
   return pctConssToBorder;
}

/** returns fraction of constraints assigned to the border for detectors in detectorchain */
void Seeed::setPctConssToBorderVector(
   std::vector<SCIP_Real> newvector
   )
{
   pctConssToBorder = newvector;
}


/** returns fraction of constraints assigned to a block for a detector */
SCIP_Real Seeed::getPctConssToBlock(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorChain.size() );

   return pctConssToBlock[detectorchainindex];
}

/** returns fraction of constraints assigned to a block for detectors in detectorchain */
std::vector<SCIP_Real> Seeed::getPctConssToBlockVector()
{
   return pctConssToBlock;
}

/** returns fraction of constraints assigned to a block for detectors in detectorchain */
void Seeed::setPctConssToBlockVector(
   std::vector<SCIP_Real> newvector  )
{
   pctConssToBlock = newvector;
}


/** returns fraction of constraints that are not longer open for a detector */
SCIP_Real Seeed::getPctConssFromFree(
   int detectorchainindex
   )
{
   assert( 0 <= detectorchainindex && detectorchainindex < (int) detectorChain.size() );

   return pctConssFromFree[detectorchainindex];
}

/** returns fraction of constraints that are not longer open for detectors in detectorchain */
std::vector<SCIP_Real> Seeed::getPctConssFromFreeVector()
{
   return pctConssFromFree;
}

/** returns index of the representative block */
int Seeed::getRepForBlock(
   int blockid
   ){
     return blockstorep[blockid];
}

std::vector<int> & Seeed::getRepVarmap(
      int repid,
      int blockrepid
      )
{
   return pidtopidvarmaptofirst[repid][blockrepid];
}


/** returns the corresponding seeedpool */
Seeedpool* Seeed::getSeeedpool()
{
   return seeedpool;
}


/** returns fraction of constraints that are not longer open for detectors in detectorchain */
void Seeed::setPctConssFromFreeVector(
   std::vector<SCIP_Real> newvector)
{
   pctConssFromFree = newvector;
}


/** returns array containing stairlinking vars */
const int* Seeed::getStairlinkingvars(
   int block
   )
{
   assert( block >= 0 && block < nBlocks );
   return & stairlinkingVars[block][0];
}

/** returns true if this seeed stems from an unpresolved problem seeed */
bool Seeed::getStemsFromUnpresolved()
{
   return stemsFromUnpresolved;
}

/** returns the data of the varclassifier that the given detector made use of */
SCIP_RETCODE Seeed::getVarClassifierData(
   int detectorchainindex,
   VarClassifier** classifier,
   std::vector<int>& varclasseslinking,
   std::vector<int>& varclassesmaster
   )
{
   assert( varClassifierUsed( detectorchainindex ) );

   *classifier = dynamic_cast<VarClassifier*>( usedClassifier[detectorchainindex] );
   varclasseslinking = classesToLinking[detectorchainindex];
   varclassesmaster = classesToMaster[detectorchainindex];

   return SCIP_OKAY;
}

/** returns array containing vars of a block */
const int* Seeed::getVarsForBlock(
   int block
   )
{
   assert( block >= 0 && block < nBlocks );
   return & varsForBlocks[block][0];
}

/** returns array containing vars of a block */
int Seeed::getVarProbindexForBlock(
   int varid,
   int block
){
   std::vector<int>::iterator lb = lower_bound( varsForBlocks[block].begin(), varsForBlocks[block].end(), varid );

   if( lb != varsForBlocks[block].end() )
      return (int) ( lb - varsForBlocks[block].begin() );
   else
      return -1;

}



/** returns true if this seeed is complete,
 *  i.e. it has at no more open constraints and variables */
bool Seeed::isComplete()
{
   return ( 0 == getNOpenconss() && 0 == getNOpenvars() );
}

/** returns true if the cons is a cons of the block */
bool Seeed::isConsBlockconsOfBlock(
   int cons,
   int block
   )
{
   assert( cons >= 0 && cons < nConss );
   assert( block >= 0 && block < nBlocks );
   std::vector<int>::iterator lb = lower_bound( conssForBlocks[block].begin(), conssForBlocks[block].end(), cons );
   if( lb != conssForBlocks[block].end() &&  *lb == cons )
      return true;
   else
      return false;
}

/** returns true if the cons is a master cons */
bool Seeed::isConsMastercons(
   int cons
   )
{
   assert( cons >= 0 && cons < nConss );
  return isconsmaster[cons];
}

/** returns true if the cons is an open cons */
bool Seeed::isConsOpencons(
   int cons
   )
{
   assert( cons >= 0 && cons < nConss );
   return isconsopen[cons];
}

/** returns true if the seeed is from a detector operating in legacymode */
bool Seeed::isFromLegacymode()
{
   return isfromlegacymode;
}

/** returns true if the seeed is from the unpresolved problem */
bool Seeed::isFromUnpresolved()
{
   return isfromunpresolved;
}


/* method to check whether this seeed is equal to given other seeed (calls isEqual(Seeed*)) */
SCIP_RETCODE Seeed::isEqual(
   Seeed* otherseeed,
   SCIP_Bool* isequal,
   bool sortseeeds
   )
{
   if( sortseeeds )
   {
      sort();
      otherseeed->sort();
   }

   * isequal = isEqual( otherseeed );

   return SCIP_OKAY;
}


/* method to check whether this seeed is equal to given other seeed */
bool Seeed::isEqual(
   Seeed* other
   )
{
   if( getNMasterconss() != other->getNMasterconss() || getNMastervars() != other->getNMastervars()
      || getNBlocks() != other->getNBlocks() || getNLinkingvars() != other->getNLinkingvars() )
      return false;

   if( getHashValue() != other->getHashValue() )
      return false;

   std::vector<std::pair<int, int>> blockorderthis = std::vector < std::pair<int, int> > ( 0 );
   std::vector<std::pair<int, int>> blockorderother = std::vector < std::pair<int, int> > ( 0 );

   /** find sorting for blocks (non decreasing according smallest row index) */
   for( int i = 0; i < this->nBlocks; ++ i )
   {
      blockorderthis.push_back( std::pair<int, int>( i, conssForBlocks[i][0] ) );
      blockorderother.push_back( std::pair<int, int>( i, other->conssForBlocks[i][0] ) );
   }

   std::sort( blockorderthis.begin(), blockorderthis.end(), compare_blocks );
   std::sort( blockorderother.begin(), blockorderother.end(), compare_blocks );

   /** compares the number of stairlinking vars */
   for( int b = 0; b < getNBlocks(); ++ b )
   {
      int blockthis = blockorderthis[b].first;
      int blockother = blockorderother[b].first;

      if( getNStairlinkingvars( blockthis ) != other->getNStairlinkingvars( blockother ) )
         return false;
   }

   /** compares the number of constraints and variables in the blocks*/
   for( int b = 0; b < getNBlocks(); ++ b )
   {
      int blockthis = blockorderthis[b].first;
      int blockother = blockorderother[b].first;

      if( ( getNVarsForBlock( blockthis ) != other->getNVarsForBlock( blockother ) )
         || ( getNConssForBlock( blockthis ) != other->getNConssForBlock( blockother ) ) )
         return false;
   }

   /** compares the master cons */
   for( int j = 0; j < getNMasterconss(); ++ j )
   {
      if( getMasterconss()[j] != other->getMasterconss()[j] )
         return false;
   }

   /** compares the master vars */
   for( int j = 0; j < getNMastervars(); ++ j )
   {
      if( getMastervars()[j] != other->getMastervars()[j] )
         return false;
   }

   /** compares the constrains and variables in the blocks */
   for( int b = 0; b < getNBlocks(); ++ b )
   {
      int blockthis = blockorderthis[b].first;
      int blockother = blockorderother[b].first;

      for( int j = 0; j < getNConssForBlock( blockthis ); ++ j )
      {
         if( getConssForBlock( blockthis )[j] != other->getConssForBlock( blockother )[j] )
            return false;
      }

      for( int j = 0; j < getNVarsForBlock( blockthis ); ++ j )
      {
         if( getVarsForBlock( blockthis )[j] != other->getVarsForBlock( blockother )[j] )
            return false;
      }

      for( int j = 0; j < getNStairlinkingvars( blockthis ); ++ j )
      {
         if( getStairlinkingvars( blockthis )[j] != other->getStairlinkingvars( blockother )[j] )
            return false;
      }
   }

   /** compares the linking vars */
   for( int j = 0; j < getNLinkingvars(); ++ j )
   {
      if( getLinkingvars()[j] != other->getLinkingvars()[j] )
         return false;
   }

   return true;
}


/** returns true if this seeed was propagated by a detector */
bool Seeed::isPropagatedBy(
   DEC_DETECTOR* detectorID
   )
{
   std::vector<DEC_DETECTOR*>::const_iterator iter = std::find( detectorChain.begin(), detectorChain.end(), detectorID );

   return iter != detectorChain.end();
}

/** returns true if this seeed is trivial,
 *  i.e. all conss are in one block, all conss are in border, all variables linking or mastervars */
bool Seeed::isTrivial()
{
   if( getNBlocks() == 1 && (SCIP_Real) getNConssForBlock( 0 ) >= 0.95 * getNConss() )
      return true;

   if( getNConss() == getNMasterconss() )
      return true;

   if( getNConss() == getNOpenconss() && getNVars() == getNOpenvars() )
      return true;

   if( getNVars() == getNMastervars() + getNLinkingvars() )
      return true;

   return false;
}


/** returns true if the seeed is selected */
bool Seeed::isSelected()
{
   return isselected;
}


/** returns true if the var is assigned to the block */
bool Seeed::isVarBlockvarOfBlock(
   int var,
   int block
   )
{
   assert( var >= 0 && var < nVars );
   assert( block >= 0 && block < nConss );

   std::vector<int>::iterator lb = lower_bound( varsForBlocks[block].begin(), varsForBlocks[block].end(), var );
   if( lb != varsForBlocks[block].end() &&  *lb == var )
      return true;
   else
      return false;
}

/** returns true if the var is a master var */
bool Seeed::isVarMastervar(
   int var
   )
{
   assert( var >= 0 && var < nVars );
  return isvarmaster[var];
}

/** returns true if the var is a linking var */
bool Seeed::isVarLinkingvar(
   int var
   )
{
   assert( var >= 0 && var < nVars );
   std::vector<int>::iterator lb = lower_bound( linkingVars.begin(), linkingVars.end(), var );
   if( lb != linkingVars.end() &&  *lb == var )
      return true;
   else
      return false;
}

/** returns true if the var is an open var */
bool Seeed::isVarOpenvar(
   int var
   )
{
   assert( var >= 0 && var < nVars );
   return isvaropen[var];
}

/** returns true if the var is a stairlinking var */
bool Seeed::isVarStairlinkingvar(
   int var
   )
{
   for( int b = 0; b < nBlocks; ++ b )
   {
      std::vector<int>::iterator lb = lower_bound( stairlinkingVars[b].begin(), stairlinkingVars[b].end(), var );
      if( lb != stairlinkingVars[b].end() &&  *lb == var )
         return true;
   }
   return false;
}

/** returns true if the var is a stairlinkingvar of the block */
bool Seeed::isVarStairlinkingvarOfBlock(
   int var,
   int block
   )
{
   assert( var >= 0 && var < nVars );
   assert( block >= 0 && block < nBlocks );
   std::vector<int>::iterator lb = lower_bound( stairlinkingVars[block].begin(), stairlinkingVars[block].end(), var );
   if( lb != stairlinkingVars[block].end() &&  *lb == var )
      return true;
   else
   {
      if( block == 0 )
         return false;
      else
      {
         lb = lower_bound( stairlinkingVars[block - 1].begin(), stairlinkingVars[block - 1].end(), var );
         return ( lb != stairlinkingVars[block-1].end() &&  *lb == var );
      }
   }
}


SCIP_RETCODE Seeed::printClassifierInformation(
   SCIP*                givenscip,
   FILE*                file
   )
{

   int nusedclassifier = (int) getNUsedClassifier();
   int nconsclassifier = 0;
   int nvarclassifier = 0;

   //SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n",  nusedclassifier );

   for( int classif = 0; classif < nusedclassifier; ++classif)
   {
      if( usedClassifier[classif] == NULL )
         continue;

      if( dynamic_cast<ConsClassifier*>( usedClassifier[classif] ) != NULL )
      {
         /** classifier is cons classifier */
         ++nconsclassifier;
      }
      else
      {
         /** classifier is var classifier */
         ++nvarclassifier;
      }
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n",  nconsclassifier );

   for( int classif = 0; classif < nusedclassifier; ++classif)
   {
      if( dynamic_cast<ConsClassifier*>( usedClassifier[classif] ) != NULL )
      {
         /** classifier is cons classifier */
         int nmasterclasses = (int) classesToMaster[classif].size();
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%s\n", usedClassifier[classif]->getName() );
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", nmasterclasses  );
         for ( int mclass = 0; mclass < (int) classesToMaster[classif].size(); ++mclass )
         {
            int classid = classesToMaster[classif][mclass];
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%s\n", usedClassifier[classif]->getClassName(classid), usedClassifier[classif]->getClassDescription(classid)  );
         }
      }
   }

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n",  nvarclassifier );

   for( int classif = 0; classif < nusedclassifier; ++classif)
   {
      if( dynamic_cast<VarClassifier*>( usedClassifier[classif] ) != NULL )
      {
         /** classifier is var classifier */
         int nmasterclasses = (int) classesToMaster[classif].size();
         int nlinkingclasses = (int) classesToLinking[classif].size();
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%s\n", usedClassifier[classif]->getName() );
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", nmasterclasses  );
         for ( int mclass = 0; mclass < (int) classesToMaster[classif].size();   ++mclass )
         {
            int classid = classesToMaster[classif][mclass];
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%s : %s\n", usedClassifier[classif]->getClassName(classid), usedClassifier[classif]->getClassDescription(classid)  );
         }

         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", nlinkingclasses  );
         for ( int linkingclass = 0; linkingclass < nlinkingclasses;   ++linkingclass )
         {
            int classid = classesToLinking[classif][linkingclass];
            SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%s : %s\n", usedClassifier[classif]->getClassName(classid),  usedClassifier[classif]->getClassDescription(classid) );
         }

      }
   }


   return SCIP_OKAY;
}


/** refine seeed with focus on blocks: assigns open conss and vars if they can be
 *  found in blocks (assignHittingOpenconss(), assignHittingOpenvars()) */
SCIP_RETCODE Seeed::refineToBlocks(
   )
{
   bool success = true;

   changedHashvalue = true;

   while( success )
      success = assignHittingOpenconss( ) || assignHittingOpenvars( );
   sort();
   return SCIP_OKAY;
}

/** refine seeed with focus on master: do obvious (considerImplicits()) assignments and
 *  assign other conss and vars to master if possible (assignOpenPartialHittingToMaster()) */
SCIP_RETCODE Seeed::refineToMaster(
    )
{
   changedHashvalue = true;

   SCIP_CALL( considerImplicits( ) );
   SCIP_CALL( assignOpenPartialHittingToMaster( ) );

   return SCIP_OKAY;
}

/** registers statistics for a used consclassifier */
void Seeed::setConsClassifierStatistics(
   int detectorchainindex,
   ConsClassifier* classifier,
   std::vector<int> consclassesmaster
   )
{
   assert( 0 <= detectorchainindex  );

   if( detectorchainindex >= (int) usedClassifier.size() )
   {
      usedClassifier.resize(detectorchainindex + 1);
      classesToMaster.resize(detectorchainindex + 1);
   }

   usedClassifier[detectorchainindex] = classifier;
   classesToMaster[detectorchainindex] = consclassesmaster;
}

/** directly adds a constraint to a block
 *  does not delete this cons from list of open conss */
SCIP_RETCODE Seeed::setConsToBlock(
   int consToBlock,
   int block
   )
{
   assert( consToBlock >= 0 && consToBlock < nConss );
   assert( block >= 0 && block < nBlocks );
   assert( (int) conssForBlocks.size() > block );

   changedHashvalue = true;

   conssForBlocks[block].push_back( consToBlock );
   conssforblocksorted = false;

   return SCIP_OKAY;
}

/** directly adds a constraint to the master constraints
 *  does not delete this cons from list of open conss */
SCIP_RETCODE Seeed::setConsToMaster(
   int consToMaster
   )
{
   assert( consToMaster >= 0 && consToMaster < nConss );
   masterConss.push_back( consToMaster );
   isconsmaster[consToMaster] = true;
   masterconsssorted = false;
   changedHashvalue = true;

   return SCIP_OKAY;
}

/** sets the whole detectorchain */
void Seeed::setDetectorchain(
   std::vector<DEC_DETECTOR*> givenDetectorChain
   )
{
   detectorChain = givenDetectorChain;
}

/** sets seeed to be propagated by a detector */
SCIP_RETCODE Seeed::setDetectorPropagated(
   DEC_DETECTOR* detectorID
   )
{
   detectorChain.push_back( detectorID );
   detectorChainFinishingUsed.push_back( FALSE );
   addEmptyClassifierStatistics();

   return SCIP_OKAY;
}

/** sets seeed to be propagated by a finishing detector */
SCIP_RETCODE Seeed::setFinishingDetectorPropagated(
   DEC_DETECTOR* detectorID
   )
{
   isFinishedByFinisher = true;
   detectorChain.push_back( detectorID );
   detectorChainFinishingUsed.push_back( TRUE );
   addEmptyClassifierStatistics();

   return SCIP_OKAY;
}

/** sets whether this seeed was finished by a detector */
void Seeed::setFinishedByFinisher(
   bool finished
   )
{
   isFinishedByFinisher = finished;
}


/** sets whether this seeed is finished by a finisher in the unpresolved problem */
void Seeed::setFinishedByFinisherUnpresolved(
   bool finishedByFinisherUnpresolved
   )
{
   isFinishedByFinisherUnpresolved = finishedByFinisherUnpresolved;
}

/** sets the detector that finished the seeed in the unpresolved problem */
void Seeed::setFinishedUnpresolvedBy(
   DEC_DETECTOR* detector
   )
{
   finishedUnpresolvedBy = detector;
}

/** sets whether this seeed stems from a detector operating in legacymode */
void Seeed::setLegacymode(
   bool legacymode
   )
{
   isfromlegacymode = legacymode;
}

/** sets number of blocks, only increasing number allowed */
SCIP_RETCODE Seeed::setNBlocks(
   int newNBlocks
   )
{
   assert( newNBlocks >= nBlocks );

   assert( (int) conssForBlocks.size() == nBlocks );
   assert( (int) varsForBlocks.size() == nBlocks );
   assert( (int) stairlinkingVars.size() == nBlocks );
   /** increase number of blocks in conssForBlocks and varsForBlocks */

   changedHashvalue = true;

   for( int b = nBlocks; b < newNBlocks; ++ b )
   {
      conssForBlocks.push_back( std::vector<int>( 0 ) );
      varsForBlocks.push_back( std::vector<int>( 0 ) );
      stairlinkingVars.push_back( std::vector<int>( 0 ) );
   }

   nBlocks = newNBlocks;

   return SCIP_OKAY;
}

/** sets the id of this seeed */
SCIP_RETCODE Seeed::setID(
   int newid
   )
{
   this->id = newid;
   return SCIP_OKAY;
}


/** sets whether this seeed is from the unpresolved problem */
void Seeed::setIsFromUnpresolved(
   bool unpresolved
   )
{
   isfromunpresolved = unpresolved;
}

/** set selection information about this seeed */
void Seeed::setSelected(
   bool selected
   )
{
   isselected = selected;
}


/** set the corresponding seeedpool */
void Seeed::setSeeedpool(
   Seeedpool* givenseeedpool
){
   this->seeedpool = givenseeedpool;
}


/** sets whether this seeed stems from an unpresolved problem seeed */
void Seeed::setStemsFromUnpresolved(
   bool stemsfromunpresolved
   )
{
   stemsFromUnpresolved = stemsfromunpresolved;
}

/** sets whether this seeed is usergiven */
void Seeed::setUsergiven(
   USERGIVEN givenusergiven
   )
{
   usergiven = givenusergiven;
}

/** registers statistics for a used varclassifier */
void Seeed::setVarClassifierStatistics(
   int detectorchainindex,
   VarClassifier* classifier,
   std::vector<int> varclasseslinking,
   std::vector<int> varclassesmaster
   )
{
   assert( 0 <= detectorchainindex );

   if( detectorchainindex >= (int) usedClassifier.size() )
    {
       usedClassifier.resize(detectorchainindex + 1);
       classesToMaster.resize(detectorchainindex + 1);
       classesToLinking.resize(detectorchainindex + 1);
    }


   usedClassifier[detectorchainindex] = classifier;
   classesToLinking[detectorchainindex] = varclasseslinking;
   classesToMaster[detectorchainindex] = varclassesmaster;
}

/** directly adds a variable to the linking variables
 *  does not delete this var from list of open vars */
SCIP_RETCODE Seeed::setVarToBlock(
   int varToBlock,
   int block
   )
{
   assert( varToBlock >= 0 && varToBlock < nVars );
   assert( block >= 0 && block < nBlocks );
   assert( (int) varsForBlocks.size() > block );

   changedHashvalue = true;

   varsForBlocks[block].push_back( varToBlock );
   varsforblocksorted = false;

   return SCIP_OKAY;
}

/** directly adds a variable to the linking variables
 *  does not delete this var from list of open vars */
SCIP_RETCODE Seeed::setVarToLinking(
   int varToLinking
   )
{
   assert( varToLinking >= 0 && varToLinking < nVars );
   linkingVars.push_back( varToLinking );
   changedHashvalue = true;
   linkingvarssorted = false;

   return SCIP_OKAY;
}

/** directly adds a variable to the master variables (hitting only constraints in the master)
 *  does not delete this var from list of open vars */
SCIP_RETCODE Seeed::setVarToMaster(
   int varToMaster
   )
{
   assert( varToMaster >= 0 && varToMaster < nVars );
   masterVars.push_back( varToMaster );
   isvarmaster[varToMaster] = true;
   mastervarssorted = false;
   changedHashvalue = true;

   return SCIP_OKAY;
}

/** directly adds a variable to the stairlinking variables
 *  does not delete this var from list of open vars */
SCIP_RETCODE Seeed::setVarToStairlinking(
   int varToStairlinking,
   int block1,
   int block2
   )
{
   assert( varToStairlinking >= 0 && varToStairlinking < nVars );
   assert( block1 >= 0 && block1 <= nBlocks );
   assert( block2 >= 0 && block2 <= nBlocks );
   assert( ( block1 + 1 == block2 ) || ( block2 + 1 == block1 ) );

   changedHashvalue = true;

   if( block1 > block2 )
      stairlinkingVars[block2].push_back( varToStairlinking );
   else
      stairlinkingVars[block1].push_back( varToStairlinking );

   stairlinkingvarsforblocksorted = false;

   return SCIP_OKAY;
}

/** generates and opens a gp visualization of the seeed */
void Seeed::showVisualisation()
{
   int returnvalue;

   MiscVisualization* miscvisu = new MiscVisualization();

   /* get names for gp file and output file */
   char filename[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   miscvisu->GCGgetVisualizationFilename(scip, this, ".gp", filename);
   miscvisu->GCGgetVisualizationFilename(scip, this, ".pdf", outname);

   /* generate gp file */
   GCGwriteGpVisualization( scip, filename, outname, getID() );

   /* compile gp file */
   char command[SCIP_MAXSTRLEN];
   strcpy(command, "gnuplot ");
   strcat(command, filename);
   SCIPinfoMessage(seeedpool->getScip(), NULL, "%s\n", command);
   returnvalue = system(command);
   if( returnvalue == -1 )
      SCIPwarningMessage(scip, "Unable to write gnuplot file\n");

   /* open outputfile */
   strcpy(command, GCGVisuGetPdfReader());
   strcat(command, " ");
   strcat(command, outname);
   strcat(command, " &");
   SCIPinfoMessage(seeedpool->getScip(), NULL, "%s\n", command);
   returnvalue = system(command);
   if( returnvalue == -1 )
      SCIPwarningMessage(scip, "Unable to open gnuplot file\n");

   return;
}

/** returns true if this seeed is a userseeed that should be completed by setting unspecified constraints to master */
SCIP_Bool Seeed::shouldCompletedByConsToMaster()
{
   return usergiven == USERGIVEN::COMPLETED_CONSTOMASTER;
}

/** sorts the vars and conss by their indices */
void Seeed::sort()
{
   for( int b = 0; b < nBlocks; ++ b )
   {
      if( ! varsforblocksorted )
         std::sort( varsForBlocks[b].begin(), varsForBlocks[b].end() );
      if( ! stairlinkingvarsforblocksorted )
         std::sort( stairlinkingVars[b].begin(), stairlinkingVars[b].end() );
      if( ! conssforblocksorted )
         std::sort( conssForBlocks[b].begin(), conssForBlocks[b].end() );
   }
   if( ! linkingvarssorted )
      std::sort( linkingVars.begin(), linkingVars.end() );
   if( !mastervarssorted )
      std::sort( masterVars.begin(), masterVars.end() );
   if( !masterconsssorted )
      std::sort( masterConss.begin(), masterConss.end() );

   varsforblocksorted = true;
   stairlinkingvarsforblocksorted = true;
   conssforblocksorted = true;
   linkingvarssorted = true;
   mastervarssorted = true;
   masterconsssorted = true;

}


/** returns a short caption for this seeed */
const char* Seeed::getShortCaption()
{
   static char shortcaption[SCIP_MAXSTRLEN];
   if( getNOpenconss() + getNOpenvars() > 0 )
      sprintf( shortcaption, "id %d; nB %d; maxW$\\geq$ %.2f ", getID(), getNBlocks(), maxwhitescore );
   else
      sprintf( shortcaption, "id %d; nB %d; maxW %.2f ", getID(), getNBlocks(), maxwhitescore );

   return shortcaption;
}


/** sets the detector chain short string */
SCIP_RETCODE Seeed::setDetectorChainString(
   char* givenDetectorchainstring
   )
{
   if ( this->detectorchainstring != NULL )
      SCIPfreeBlockMemoryArray(scip, & detectorchainstring, SCIP_MAXSTRLEN);
   SCIP_CALL( SCIPduplicateBlockMemoryArray( scip, & this->detectorchainstring, givenDetectorchainstring, SCIP_MAXSTRLEN ) );
   return SCIP_OKAY;
}

/**
 * finds  a translation from orig to transformed seeedpool
 * @param consindex already allocated
 * @param varindex already allocated
 */
SCIP_RETCODE findTranslationForDec(
   Seeedpool*        origseeedpool,
   Seeedpool*        transseeedpool,
   std::vector<int>* consindex,
   std::vector<int>* varindex,
   SCIP_Bool         fromunpresolvedtopresolved,
   SCIP_Bool*        success
)
{
   SCIP* scip;
   int norigconss;
   int norigvars;
   SCIP_CONS** origconss;
   SCIP_VAR**  origvars;

   scip = origseeedpool->getScip();
   norigconss = SCIPgetNOrigConss(scip);
   origconss = SCIPgetOrigConss(scip);
   norigvars = SCIPgetNOrigVars(scip);
   origvars = SCIPgetOrigVars(scip);

   for( size_t oc = 0; oc < consindex->size(); ++oc )
   {
      consindex->at(oc) = -1;
   }

   for( size_t ov = 0; ov < varindex->size(); ++ov )
   {
      varindex->at(ov) = -1;
   }

   for( int oc = 0; oc < norigconss; ++oc )
   {
      SCIP_CONS* origcons;
      SCIP_CONS* transcons;
      int origconsid;
      int transconsid;

      origconsid = -1;
      transconsid = -1;

      origcons = origconss[oc];
      origconsid = origseeedpool->getIndexForCons(origcons);

      SCIPgetTransformedCons(scip, origcons, &transcons);
      if( transcons == NULL )
      {
 //       std::cout << "consname: " << SCIPconsGetName(origcons) << " ; oc:" << oc << " has no transformed constraint "  << std::endl;
        continue;
      }
      transconsid = transseeedpool->getIndexForCons(transcons);

 //     std::cout << "consname: " << SCIPconsGetName(origcons) << " ; oc:" << oc << " ;transconsid: " << transconsid << " ; origconsid: " << origconsid << " transformed: " << origseeedpool->getTransformedInfo() << std::endl;

      if( fromunpresolvedtopresolved  )
      {
         consindex->at(origconsid) = transconsid;
      }
      else
      {
         consindex->at(transconsid) = origconsid;
      }
   }

   for( int ov = 0; ov < norigvars; ++ov )
   {
      SCIP_VAR* origvar;
      SCIP_VAR* transvar;
      int origvarid;
      int transvarid;

      origvarid = -1;
      transvarid = -1;

      origvar = origvars[ov];
      origvarid = origseeedpool->getIndexForVar(origvar);
      SCIPgetTransformedVar(scip, origvar, &transvar);

      if( transvar == NULL)
      {
         continue;
      }


      transvarid = transseeedpool->getIndexForVar( SCIPvarGetProbvar(transvar) );

      if( fromunpresolvedtopresolved  )
      {
         varindex->at(origvarid) = transvarid;
      }
      else
      {
         varindex->at(transvarid) = origvarid;
      }
   }

   *success = TRUE;

   return SCIP_OKAY;

}

SCIP_RETCODE Seeed::writeAsDec(
   FILE* file,
   Seeedpool*   seeedpooltowriteto,
   SCIP_RESULT*  result
   )
{

   Seeed*      helpseeed;
   static const char commentchars[] = "\\";

   int nconss;
   int nvars;
   std::vector<int> consindex(0);
   std::vector<int> varindex(0);

   assert(seeedpooltowriteto != NULL);

   helpseeed = this;
   nconss = seeedpooltowriteto->getNConss();
   nvars = seeedpooltowriteto->getNVars();

   consindex = std::vector<int>(nconss);
   varindex = std::vector<int>(nvars);

   /* is there no translation needed ? */
   if( getSeeedpool() == seeedpooltowriteto )
   {
      for( int i = 0; i < nconss; ++i )
         consindex[i] = i;
      for( int i = 0; i < nvars; ++i )
         varindex[i] = i;
   }
   else /** translation is neeeded */
   {
      SCIP_Bool success;

      Seeed* transseeed;
      success = FALSE;

      if( isFromUnpresolved() )
         findTranslationForDec(getSeeedpool(), seeedpooltowriteto, &consindex, &varindex, isFromUnpresolved() , &success );
      else
         findTranslationForDec( seeedpooltowriteto, getSeeedpool(), &consindex, &varindex, isFromUnpresolved() , &success );

      transseeed = new Seeed(seeedpooltowriteto->getScip(), -1, seeedpooltowriteto );

      /** translate seeed */
      transseeed->setNBlocks(getNBlocks() );

      for( int b = 0; b < getNBlocks(); ++b )
      {
         for( int c = 0; c < getNConssForBlock(b); ++c )
         {
            int cons = consindex[getConssForBlock(b)[c] ];
            if( cons == -1 )
               continue;
            transseeed->bookAsBlockCons(cons, b);
         }

         for( int v = 0; v < getNVarsForBlock(b); ++v )
         {
            int var = varindex[getVarsForBlock(b)[v] ];
            if( var == -1 )
               continue;
            transseeed->bookAsBlockVar(var, b);
         }

         for( int v = 0; v < getNStairlinkingvars(b); ++v )
         {
            int var = varindex[getStairlinkingvars(b)[v] ];
            if( var == -1 )
               continue;
            transseeed->bookAsStairlinkingVar(var, b);
         }
      }

      for( int c = 0; c < getNMasterconss(); ++c )
      {
         int cons = consindex[getMasterconss()[c] ];
         if( cons == -1 )
            continue;
         transseeed->bookAsMasterCons(cons);
      }

      for( int v = 0; v < getNLinkingvars(); ++v )
      {
         int var = varindex[getLinkingvars()[v] ];
         if( var == -1 )
            continue;
         transseeed->bookAsLinkingVar(var);
      }

      for( int v = 0; v < getNMastervars(); ++v )
      {
         int var = varindex[getMastervars()[v] ];
         if( var == -1 )
            continue;
         transseeed->bookAsMasterVar(var);
      }

      transseeed->flushBooked();

      transseeed->considerImplicits();

      transseeed->deleteEmptyBlocks(false);

      //displayInfo(getSeeedpool(), 0 );

      //transseeed->displayInfo(seeedpooltowriteto, 0 );

      if( transseeed->isComplete() != isComplete() )
         success = FALSE;

      if( ! success )
      {
         if( isFromUnpresolved() )
            SCIPwarningMessage(seeedpooltowriteto->getScip(), "Writing dec-file is not possible since translation to presolved (transformed) problem failed. Please consider writing for original problem. Ignore next message about written problem if it is there.\n" );
         else
            SCIPwarningMessage(seeedpooltowriteto->getScip(), "Writing dec-file is not possible since translation to unpresolved (non-transformed) problem failed. Please consider writing for transformed problem. Ignore next message about written problem if it is there.\n" );
         /** unforunately there is no appropiate result type */
         *result = SCIP_SUCCESS;
         return SCIP_OKAY;
      }
      else
         *result = SCIP_SUCCESS;

      helpseeed = transseeed;

   }

   /** @TODO: statistical stuff  */
   /* at first: write meta data of decomposition as comment */
   SCIPinfoMessage(scip, file, "%s%s ndetectors \n", commentchars, commentchars );
   SCIPinfoMessage(scip, file, "%s%s %d \n", commentchars, commentchars, getNDetectorchainInfo() );

   SCIPinfoMessage(scip, file, "%s%s name time nnewblocks %%ofnewborderconss %%ofnewblockconss %%ofnewlinkingvars %%ofnewblockvars  \n", commentchars, commentchars );

   for ( int i = 0; i < getNDetectorchainInfo() ; ++i )
   {
      SCIPinfoMessage(scip, file, "%s%s %s %f %d %f %f %f %f \n", commentchars, commentchars, DECdetectorGetName(getDetectorchain()[i] ), detectorClockTimes[i],
         nNewBlocks[i], pctConssToBorder[i], pctConssToBlock[i], pctVarsToBorder[i],
         pctVarsToBlock[i]) ;
   }


   if( !helpseeed->isComplete() )
         SCIPinfoMessage(scip, file, "INCOMPLETE\n1\n" );

   if( ( isFromUnpresolved() && seeedpooltowriteto == getSeeedpool() ) || ( !isFromUnpresolved() && seeedpooltowriteto != getSeeedpool() )  )
      SCIPinfoMessage(scip, file, "PRESOLVED\n0\n" );
   else
      SCIPinfoMessage(scip, file, "PRESOLVED\n1\n" );

   SCIPinfoMessage(scip, file, "NBLOCKS\n%d\n", getNBlocks() );


   for( int b = 0; b < helpseeed->getNBlocks(); ++b )
   {
      SCIPinfoMessage(scip, file, "BLOCK %d\n", b+1 );
      for( size_t c = 0; c < helpseeed->conssForBlocks[b].size(); ++c )
      {
         SCIPinfoMessage(scip, file, "%s\n", SCIPconsGetName(seeedpooltowriteto->getConsForIndex(  helpseeed->conssForBlocks[b][c]  )) );
      }
   }

    SCIPinfoMessage(scip, file, "MASTERCONSS\n" );
   for( int mc = 0; mc < helpseeed->getNMasterconss(); ++mc )
   {
      SCIPinfoMessage(scip, file, "%s\n", SCIPconsGetName(seeedpooltowriteto->getConsForIndex( helpseeed->masterConss[mc])) );
   }

   if( helpseeed->isComplete() )
   {
      if( this != helpseeed )
         delete helpseeed;
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }
   for( int b = 0; b < helpseeed->getNBlocks(); ++b )
   {
      SCIPinfoMessage(scip, file, "BLOCKVARS %d\n", b+1 );
      for( size_t v = 0; v < helpseeed->varsForBlocks[b].size(); ++v )
      {
         SCIPinfoMessage(scip, file, "%s\n", SCIPvarGetName(seeedpooltowriteto->getVarForIndex( helpseeed->varsForBlocks[b][v])) );
      }
   }

   SCIPinfoMessage(scip, file, "LINKINGVARS\n" );
   for( int lv = 0; lv < helpseeed->getNLinkingvars(); ++lv )
   {
      SCIPinfoMessage(scip, file, "%s\n", SCIPvarGetName(seeedpooltowriteto->getVarForIndex( helpseeed->linkingVars[lv])) );
   }

   SCIPinfoMessage(scip, file, "MASTERVARS\n" );
   for( int mv = 0; mv < helpseeed->getNMastervars(); ++mv )
   {
      SCIPinfoMessage(scip, file, "%s\n", SCIPvarGetName(seeedpooltowriteto->getVarForIndex( helpseeed->masterVars[mv])) );
   }

   if( this != helpseeed )
      delete helpseeed;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;

}


/** creates and sets a detector chain short string for this seeed */
SCIP_RETCODE Seeed::buildDecChainString()
{
   char decchaininfo[SCIP_MAXSTRLEN];
   /** set detector chain info string */
   SCIPsnprintf( decchaininfo, SCIP_MAXSTRLEN, "" );
   if( this->usergiven == USERGIVEN::PARTIAL || this->usergiven == USERGIVEN::COMPLETE
      || this->usergiven == USERGIVEN::COMPLETED_CONSTOMASTER || this->getDetectorchain() == NULL
      || this->getDetectorchain()[0] == NULL )
   {
      char str1[2] = "\0"; /* gives {\0, \0} */
      str1[0] = 'U';
      (void) strncat( decchaininfo, str1, 1 );
   }
   for( int d = 0; d < this->getNDetectors(); ++ d )
   {
      if( d == 0 && this->getDetectorchain()[d] == NULL )
         continue;
      char str[2] = "\0"; /* gives {\0, \0} */
      str[0] = DECdetectorGetChar( this->getDetectorchain()[d] );
      (void) strncat( decchaininfo, str, 1 );
   }

   SCIP_CALL( this->setDetectorChainString( decchaininfo ) );

   return SCIP_OKAY;
}

void Seeed::calcmaxwhitescore(){

   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );

   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );

   if( blockareascore == -1. )
      calcblockareascore();

   if( borderareascore == -1. )
      calcborderareascore();

   /** maxwhitescore = 1 - ( 1 - blackareascore + (1 - borderareascore ) ) */
   maxwhitescore = blockareascore + borderareascore - 1.;



   if( maxwhitescore < 0. )
     maxwhitescore = 0.;

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return;
}

void Seeed::calcbenderareascore()
{
   unsigned long nrelevantconss;
   unsigned long nrelevantvars;

   unsigned long nrelevantconss2;
   unsigned long nrelevantvars2;

   unsigned long badblockvararea;


   long benderborderarea;
   unsigned long totalarea;

   nrelevantconss = 0;
   nrelevantvars = 0;

   nrelevantconss2 = 0;
   nrelevantvars2 = 0;

   badblockvararea = 0;

   for( int  c = 0; c < getNMasterconss(); ++c )
   {
      bool relevant = true;
      int cons = getMasterconss()[c];
      for( int v = 0; v < seeedpool->getNVarsForCons(cons); ++v )
      {
         int var = seeedpool->getVarsForCons(cons)[v];
         if ( isVarOpenvar(var) || isVarMastervar(var) || isVarLinkingvar(var) )
         {
            relevant = false;
            break;
         }

      }
      if( relevant )
         ++nrelevantconss;
   }

   for( int b = 0; b < getNBlocks(); ++b )
   {
      for(int v = 0; v < getNVarsForBlock(b); ++v )
      {
         bool relevant = true;
         int var = getVarsForBlock(b)[v];
         for( int c = 0; c < seeedpool->getNConssForVar(var); ++c )
         {
            int cons = seeedpool->getConssForVar(var)[c];
            if( isConsMastercons(cons) || isConsOpencons(cons)  )
            {
               relevant  = false;
               for( int b2 = 0; b2 < getNBlocks(); ++b2 )
               {
                  if( b2 != b )
                     badblockvararea += getNConssForBlock(b2);
               }
               break;
            }
         }
         if( relevant )
            ++nrelevantvars;
      }
   }

   for( int  v = 0; v < getNLinkingvars(); ++v )
   {
      bool relevant = true;
      int var = getLinkingvars()[v];
      for( int c = 0; c < seeedpool->getNConssForVar(var); ++c )
      {
         int cons = seeedpool->getConssForVar(var)[c];
         if ( isConsOpencons(cons) || isConsMastercons(cons) )
         {
            relevant = false;
            break;
         }

      }
      if( relevant )
         ++nrelevantvars2;
   }


   for( int b = 0; b < getNBlocks(); ++b )
   {
      for(int c = 0; c < getNConssForBlock(b); ++c )
      {
         bool relevant = true;
         int cons = getConssForBlock(b)[c];
         for( int v = 0; v < seeedpool->getNVarsForCons(cons); ++v )
         {
            int var = seeedpool->getVarsForCons(cons)[v];
            if( isVarLinkingvar(var) || isVarOpenvar(var)  )
            {
               relevant  = false;
               break;
            }
         }
         if( relevant )
            ++nrelevantconss2;
      }
   }



   benderborderarea = ( nrelevantconss * nrelevantvars  ) + ( nrelevantconss2 * nrelevantvars2  ) - badblockvararea;
   totalarea = ( (unsigned long) getNConss() * (unsigned long) getNVars() );
   benderareascore =  ( SCIP_Real) benderborderarea / totalarea;

}

void Seeed::calcbendersscore(){

   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );

   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );

   if( blockareascore == -1. )
      calcblockareascore();

   if( benderareascore == -1. )
      calcbenderareascore();

   if( borderareascore == -1. )
      calcborderareascore();

   /** bendersscore = 1 - ( 1 - blockareascore + (1 - borderareascore - benderborderscore ) ) */
   bendersscore = blockareascore + benderareascore + borderareascore - 1.;

    if( bendersscore < 0. )
     bendersscore = 0.;

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return;
}


SCIP_RETCODE Seeed::calcclassicscore()
{
   int i;
   int j;
   int k;

   unsigned long matrixarea;
   unsigned long borderarea;
   SCIP_Real borderscore; /**< score of the border */
   SCIP_Real densityscore; /**< score of block densities */
   SCIP_Real linkingscore; /**< score related to interlinking blocks */
   SCIP_Real totalscore; /**< accumulated score */

   SCIP_Real varratio;
   int* nzblocks;
   int* nlinkvarsblocks;
   int* nvarsblocks;
   SCIP_Real* blockdensities;
   int* blocksizes;
   SCIP_Real density;

   SCIP_Real alphaborderarea;
   SCIP_Real alphalinking;
   SCIP_Real alphadensity;

   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );



   SCIP_CALL( SCIPallocBufferArray( scip, & nzblocks, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & nlinkvarsblocks, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & blockdensities, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & blocksizes, nBlocks ) );
   SCIP_CALL( SCIPallocBufferArray( scip, & nvarsblocks, nBlocks ) );

   alphaborderarea = 0.6;
   alphalinking = 0.2;
   alphadensity = 0.2;


   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate slave sizes, nonzeros and linkingvars */
   for( i = 0; i < nBlocks; ++ i )
   {
      int ncurconss;
      int nvarsblock;
      SCIP_Bool *ishandled;

      SCIP_CALL( SCIPallocBufferArray( scip, & ishandled, nVars ) );
      nvarsblock = 0;
      nzblocks[i] = 0;
      nlinkvarsblocks[i] = 0;

      //    std::cout << "blackarea =  " << blackarea << " +  " << getNConssForBlock( i ) << " * " << getNVarsForBlock( i ) << " = " << getNConssForBlock( i ) * ( getNVarsForBlock( i ) );

   //   blackarea += (unsigned long) getNConssForBlock( i ) * ( (unsigned long) getNVarsForBlock( i ) );
      //  std::cout << " =  " << blackarea  << std::endl;

      for( j = 0; j < nVars; ++ j )
      {
         ishandled[j] = FALSE;
      }
      ncurconss = getNConssForBlock( i );

      for( j = 0; j < ncurconss; ++ j )
      {
         int cons = getConssForBlock( i )[j];
         int ncurvars;
         ncurvars = seeedpool->getNVarsForCons( cons );
         for( k = 0; k < ncurvars; ++ k )
         {
            int var = seeedpool->getVarsForCons( cons )[k];
            int block = -3;
            if( isVarBlockvarOfBlock( var, i ) )
               block = i + 1;
            else if( isVarLinkingvar( var ) || isVarStairlinkingvar( var ) )
               block = nBlocks + 2;
            else if( isVarMastervar( var ) )
               block = nBlocks + 1;

            ++ ( nzblocks[i] );

            if( block == nBlocks + 1 && ishandled[var] == FALSE )
            {
               ++ ( nlinkvarsblocks[i] );
            }
            ishandled[var] = TRUE;
         }
      }

      for( j = 0; j < nVars; ++ j )
      {
         if( ishandled[j] )
         {
            ++ nvarsblock;
         }
      }

      blocksizes[i] = nvarsblock * ncurconss;
      nvarsblocks[i] = nvarsblock;
      if( blocksizes[i] > 0 )
      {
         blockdensities[i] = 1.0 * nzblocks[i] / blocksizes[i];
      }
      else
      {
         blockdensities[i] = 0.0;
      }

      assert( blockdensities[i] >= 0 && blockdensities[i] <= 1.0 );
      SCIPfreeBufferArray( scip, & ishandled );
   }


   borderarea = ((unsigned long) getNMasterconss() * nVars )  + ( ((unsigned long) getNLinkingvars() + getNMastervars() + getNTotalStairlinkingvars() ) ) * ( nConss - getNMasterconss() );

   matrixarea = ((unsigned long) nVars ) * ((unsigned long) nConss );

//   std::cout << "black area ration =  " << blackarea << "/ ( " << getNConss() << " * " << getNVars() << " =  " << ( (unsigned long) getNConss() * (unsigned long) getNVars() ) << ")  = " << maxwhitescore << std::endl;

//std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    this seeed has a black area ratio of " << maxwhitescore << std::endl;

   density = 1E20;
   varratio = 1.0;
   linkingscore = 1.;
   borderscore =  1.;
   densityscore = 1.;

   for( i = 0; i < nBlocks; ++ i )
   {
      density = MIN( density, blockdensities[i] );

      if( ( getNLinkingvars() + getNMastervars() + getNTotalStairlinkingvars() ) > 0 )
      {
         varratio *= 1.0 * nlinkvarsblocks[i] / ( getNLinkingvars() + getNMastervars() + getNTotalStairlinkingvars() );
      }
      else
      {
         varratio = 0.;
      }
   }
   linkingscore = ( 0.5 + 0.5 * varratio );

   densityscore = ( 1. - density );

   borderscore = ( 1.0 * ( borderarea ) / matrixarea );

   totalscore = 1. - (alphaborderarea * ( borderscore ) + alphalinking * ( linkingscore ) + alphadensity * ( densityscore ) );

   score = totalscore;

   SCIPfreeBufferArray( scip, & nzblocks );
   SCIPfreeBufferArray(  scip, & nlinkvarsblocks) ;
   SCIPfreeBufferArray(  scip, & blockdensities);
   SCIPfreeBufferArray(  scip, & blocksizes);
   SCIPfreeBufferArray(  scip, & nvarsblocks);

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return SCIP_OKAY;
}

void Seeed::calcborderareascore(){

   unsigned long matrixarea;
   unsigned long borderarea;


   matrixarea = (unsigned long) getNVars() * (unsigned long)getNConss();
   borderarea = 0;

   borderarea += (unsigned long) ( getNLinkingvars() + getNTotalStairlinkingvars() ) * (unsigned long) getNConss();
   borderarea += (unsigned long) getNMasterconss() * ( (unsigned long) getNVars() - ( getNLinkingvars() + getNTotalStairlinkingvars() ) ) ;

   borderareascore = 1. - ( (SCIP_Real) borderarea / (SCIP_Real) matrixarea );

   return;
}

void Seeed::calcmaxforeseeingwhitescore(){

   std::vector<int> nlinkingvarsforblock(getNBlocks(), 0);
   std::vector<int> nblocksforlinkingvar(getNLinkingvars() + getNTotalStairlinkingvars(), 0);

   unsigned long sumblockshittinglinkingvar;
   unsigned long sumlinkingvarshittingblock;
   unsigned long newheight;
   unsigned long newwidth;
   unsigned long newmasterarea;
   unsigned long newblockarea;

   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );

   for( int lv = 0; lv < getNLinkingvars(); ++lv )
   {
      int linkingvarid = getLinkingvars()[lv];

      for( int b = 0; b < getNBlocks(); ++b )
      {
         for ( int blc = 0; blc < getNConssForBlock(b); ++blc )
         {
            int blockcons = getConssForBlock(b)[blc];
            if( !SCIPisZero( seeedpool->getScip(), seeedpool->getVal(blockcons, linkingvarid) ) )
            {
               /** linking var hits block */
               ++nlinkingvarsforblock[b];
               ++nblocksforlinkingvar[lv];
               break;
            }
         }
      }
   }

   for( int b = 0; b < getNBlocks(); ++b)
   {
      for( int slv = 0; slv < getNStairlinkingvars(b); ++slv )
      {
         ++nlinkingvarsforblock[b];
         ++nlinkingvarsforblock[b+1];
         ++nblocksforlinkingvar[getNLinkingvars() + slv];
         ++nblocksforlinkingvar[getNLinkingvars() + slv];
      }
   }

   sumblockshittinglinkingvar = 0;
   sumlinkingvarshittingblock = 0;
   for( int b = 0; b < getNBlocks(); ++b )
   {
      sumlinkingvarshittingblock += nlinkingvarsforblock[b];
   }
   for( int lv = 0; lv < getNLinkingvars(); ++lv )
   {
      sumblockshittinglinkingvar += nblocksforlinkingvar[lv];
   }

   for( int slv = 0; slv < getNTotalStairlinkingvars(); ++slv )
   {
      sumblockshittinglinkingvar += nblocksforlinkingvar[getNLinkingvars() + slv];
   }


   newheight = getNConss() + sumblockshittinglinkingvar;
   newwidth = getNVars() + sumlinkingvarshittingblock;

   newmasterarea = ( getNMasterconss() + sumblockshittinglinkingvar) * ( getNVars() + sumlinkingvarshittingblock );
   newblockarea = 0;

   for( int b = 0; b < getNBlocks(); ++b )
   {
      newblockarea += getNConssForBlock(b) * ( getNVarsForBlock(b) + nlinkingvarsforblock[b] );
   }

   maxforeseeingwhitescore = ((SCIP_Real ) newblockarea + (SCIP_Real) newmasterarea) / (SCIP_Real) newwidth;
   maxforeseeingwhitescore =  maxforeseeingwhitescore / (SCIP_Real) newheight ;

   maxforeseeingwhitescore = 1. - maxforeseeingwhitescore;

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return;
}

void Seeed::calcmaxforeseeingwhitescoreagg(){

   std::vector<int> nlinkingvarsforblock(getNBlocks(), 0);
   std::vector<int> nblocksforlinkingvar(getNLinkingvars() + getNTotalStairlinkingvars(), 0);

   unsigned long sumblockshittinglinkingvar;
   unsigned long sumlinkingvarshittingblock;
   unsigned long newheight;
   unsigned long newwidth;
   unsigned long newmasterarea;
   unsigned long newblockareaagg;

   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );

   calcAggregationInformation();

   for( int lv = 0; lv < getNLinkingvars(); ++lv )
   {
      int linkingvarid = getLinkingvars()[lv];

      for( int b = 0; b < getNBlocks(); ++b )
      {
         for ( int blc = 0; blc < getNConssForBlock(b); ++blc )
         {
            int blockcons = getConssForBlock(b)[blc];
            if( !SCIPisZero( seeedpool->getScip(), seeedpool->getVal(blockcons, linkingvarid) ) )
            {
               /** linking var hits block */
               ++nlinkingvarsforblock[b];
               ++nblocksforlinkingvar[lv];
               break;
            }
         }
      }
   }

   for( int b = 0; b < getNBlocks(); ++b)
   {
      for( int slv = 0; slv < getNStairlinkingvars(b); ++slv )
      {
         ++nlinkingvarsforblock[b];
         ++nlinkingvarsforblock[b+1];
         ++nblocksforlinkingvar[getNLinkingvars() + slv];
         ++nblocksforlinkingvar[getNLinkingvars() + slv];
      }
   }

   sumblockshittinglinkingvar = 0;
   sumlinkingvarshittingblock = 0;
   for( int b = 0; b < getNBlocks(); ++b )
   {
      sumlinkingvarshittingblock += nlinkingvarsforblock[b];
   }
   for( int lv = 0; lv < getNLinkingvars(); ++lv )
   {
      sumblockshittinglinkingvar += nblocksforlinkingvar[lv];
   }

   for( int slv = 0; slv < getNTotalStairlinkingvars(); ++slv )
   {
      sumblockshittinglinkingvar += nblocksforlinkingvar[getNLinkingvars() + slv];
   }


   newheight = getNConss() + sumblockshittinglinkingvar;
   newwidth = getNVars() + sumlinkingvarshittingblock;

   newmasterarea = ( getNMasterconss() + sumblockshittinglinkingvar) * ( getNVars() + sumlinkingvarshittingblock );
   newblockareaagg = 0;

   assert(nrepblocks > 0 );
   for( int br = 0; br < nrepblocks; ++br )
   {
      newblockareaagg += getNConssForBlock( reptoblocks[br][0] ) * ( getNVarsForBlock( reptoblocks[br][0] ) + nlinkingvarsforblock[reptoblocks[br][0]] );
   }

   maxforeseeingwhitescoreagg = ((SCIP_Real ) newblockareaagg + (SCIP_Real) newmasterarea) / (SCIP_Real) newwidth;
   maxforeseeingwhitescoreagg =  maxforeseeingwhitescoreagg / (SCIP_Real) newheight ;

   maxforeseeingwhitescoreagg = 1. - maxforeseeingwhitescoreagg;

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return;
}

void Seeed::calcsetpartfwhitescore(){

   if( maxforeseeingwhitescore == -1. )
      calcmaxforeseeingwhitescore();
   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );

   if( hasSetppccardMaster() && !isTrivial() && getNBlocks() > 1 )
   {
      setpartfwhitescore = 0.5 * maxforeseeingwhitescore + 0.5;
   }
   else
   {
      setpartfwhitescore = 0.5 * maxforeseeingwhitescore;
   }

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return;
}

void Seeed::calcsetpartfwhitescoreagg(){

   if( maxforeseeingwhitescoreagg == -1. )
      calcmaxforeseeingwhitescoreagg();

   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );

   if( hasSetppccardMaster() && !isTrivial() && getNBlocks() > 1 )
   {
      setpartfwhitescoreagg = 0.5 * maxforeseeingwhitescoreagg + 0.5;
   }
   else
   {
      setpartfwhitescoreagg = 0.5 * maxforeseeingwhitescoreagg;
   }

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return;
}

void Seeed::calcblockareascore(){

   unsigned long matrixarea;
   unsigned long blockarea;


   matrixarea = (unsigned long) getNVars()  * (unsigned long) getNConss() ;
   blockarea = 0;

   for( int i = 0; i < getNBlocks(); ++ i )
   {
      blockarea += (unsigned long) getNConssForBlock( i ) * ( (unsigned long) getNVarsForBlock( i ) );
   }

   blockareascore = 1. - ( (SCIP_Real) blockarea / (SCIP_Real) matrixarea );

   return;
}

void Seeed::calcblockareascoreagg(){

   unsigned long matrixarea;
   unsigned long blockarea;

   SCIP_CLOCK* clock;

   SCIP_CALL_ABORT( SCIPcreateClock( seeedpool->getScip(), &clock) );
   SCIP_CALL_ABORT( SCIPstartClock( seeedpool->getScip(), clock) );

   matrixarea = (unsigned long)getNVars() * (unsigned long)getNConss();
   blockarea = 0;

   for( int i = 0; i < nrepblocks; ++ i )
   {
      blockarea += (unsigned long) getNConssForBlock( reptoblocks[i][0] ) * ( (unsigned long) getNVarsForBlock( reptoblocks[i][0] ) );
   }

   blockareascoreagg = 1. - ( (SCIP_Real) blockarea / (SCIP_Real) matrixarea );

   SCIP_CALL_ABORT(SCIPstopClock( seeedpool->getScip(), clock) );
   seeedpool->scorecalculatingtime += SCIPgetClockTime( seeedpool->getScip(), clock);
   SCIP_CALL_ABORT(SCIPfreeClock( seeedpool->getScip(), &clock) );

   return;
}

} /* namespace gcg */
