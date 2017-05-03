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

/**@file   class_Seeed.h
 * @brief  class with functions for seeed
 * @author Michael Bastubbe
 * @author Hannah Hechenrieder
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_SEEED_H__
#define GCG_CLASS_SEEED_H__


#include "objscip/objscip.h"
#include <vector>
#include "struct_detector.h"
#include <string>



namespace gcg {

enum USERGIVEN
{
   NOT = 0,
   PARTIAL = -1,
   COMPLETE = -2,
   COMPLETED_CONSTOMASTER = -3
};

class Seeedpool;


class Seeed
{

private:
	SCIP*							         scip;
   int								      id;						         /**< id of the seeed */
   int 								      nBlocks;				            /**< number of blocks the decomposition currently has */
   int 								      nVars;                        /**< number of variables */
   int								      nConss;                       /**< number of constraints */
   std::vector<int>					   masterConss;			         /**< vector containing indices of master constraints */
   std::vector<int>					   masterVars;				         /**< vector containing indices of master variables */
   std::vector<std::vector<int>>    conssForBlocks; 		         /**< conssForBlocks[k] contains a vector of indices of all constraints assigned to block k */
   std::vector<std::vector<int>>    varsForBlocks; 			      /**< varsForBlocks[k] contains a vector of indices of all variables assigned to block k */
   std::vector<int> 				      linkingVars;			         /**< vector containing indices of linking variables */
   std::vector<std::vector<int>>    stairlinkingVars;		         /**< vector containing indices of staircase linking variables of the blocks (stairlinking variables can be found only in the vector of their first block) */
   std::vector<int> 				      openVars;				         /**< vector containing indices of variables that are not assigned yet*/
   std::vector<int> 				      openConss;				         /**< vector containing indices of constraints that are not assigned yet*/
   std::vector<int>                 bookedAsMasterConss;          /**< vector containing indices of constraints that are not assigned yet but booked as master conss */
   std::vector<std::pair<int, int>> bookedAsBlockConss;           /**< vector containing indices of constraints that are not assigned yet but booked as block conss and the block */
   std::vector<int>                 bookedAsLinkingVars;          /**< vector containing indices of variables that are not assigned yet but booked as linking vars */
   std::vector<int>                 bookedAsMasterVars;           /**< vector containing indices of variables that are not assigned yet but booked as master vars */
   std::vector<std::pair<int, int>> bookedAsBlockVars;            /**< vector containing indices of variables that are not assigned yet but booked as block vars and the block */
   std::vector<std::pair<int, int>> bookedAsStairlinkingVars;     /**< vector containing indices of variables that are not assigned yet but booked as stairlinking vars and the first block of the stairlinking var */
   std::vector<bool> 				   propagatedByDetector;	      /**< propagatedByDetector[i] is this seeed propagated by detector i */
   bool 							         openVarsAndConssCalculated;   /**< are the openVars and openCons calculated */
   long                             hashvalue;
   SCIP_Real                        score;                        /**< score to evaluate the seeeds */
   SCIP_Real                        maxwhitescore;                /**< score corresponding to the max white measure */

   bool                             changedHashvalue;             /**< are there any changes concerning the hash value since it was calculated last time */

   const static int              primes[];
   const static int              nPrimes;

public:

   bool                             isFinishedByFinisher;         /**< was this seeed finished by the finishseeed() method of a detector */
   /** statistic information */
   std::vector<DEC_DETECTOR*>       detectorChain;              /**< vector containing detectors that worked on that seeed */
   std::vector<std::string>         detectorchaininfo;          /**< vector containing information about the detector call */
   std::vector<SCIP_Bool>           detectorChainFinishingUsed; /**< vector containing whether the finishing method of the corresponding detector was used on that seeed */
   std::vector<SCIP_Real>           detectorClockTimes;         /**< vector containing detector times in seconds  */
   std::vector<SCIP_Real>           pctVarsToBorder;            /**< vector containing the fraction of variables assigned to the border for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctVarsToBlock;             /**< vector containing the fraction of variables assigned to a block for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctVarsFromFree;            /**< vector containing the fraction of variables that are not longer open for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctConssToBorder;           /**< vector containing the fraction of constraints assigned to the border for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctConssToBlock;            /**< vector containing the fraction of constraints assigned to a block for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctConssFromFree;           /**< vector containing the fraction of constraints that are not longer open for each detector working on that seeed*/
   std::vector<int>                 nNewBlocks;                 /**< vector containing detector indices that worked on that seeed */

   std::vector<int>                 listofancestorids;          /**< vector containing detector indices that worked on that seeed */

   USERGIVEN                        usergiven;                  /**< is this seeed partially or completely given by user */

   /** datastructure to store information if this seeed stems from a seeed concerning the unpresolved problem */
   bool                             stemsFromUnpresolved;
   bool                             isFinishedByFinisherUnpresolved; /**< was the ancestor seeed for the unpresolved problem finished by the finishseeed() method of a detector */
   DEC_DETECTOR*                    finishedUnpresolvedBy;           /**< index of dinishing detector of unpresolved ancestor seeed */


   /** constructor(s) */
   Seeed(
      SCIP*          scip,
	  int             id,      		   	/**< id that is given to this seeed */
	  int             nDetectors,          /**< number of detectors */
	  int				   nConss,				   /**number of constraints */
	  int 				nVars				      /**number of variables */
      );


   Seeed(const Seeed *seeedToCopy, Seeedpool* seeedpool);

   ~Seeed();


   /** adds a block, returns the number of the new block */
   int addBlock(
   );

   /** incorporates the the needed time of a certain detector in the detector chain */
   void addClockTime(
      SCIP_Real clocktime
   );

   /** incorporates the changes from ancestor seeed */
   void addDecChangesFromAncestor(
      Seeed* ancestor
   );

   /** adds a detector chain info */
   void addDetectorChainInfo(
      const char* decinfo
   );

   /** returns true if at least one constraint is assigned to a block */
   bool alreadyAssignedConssToBlocks();

   /** returns true if the open vars and conss are calculated */
   bool areOpenVarsAndConssCalculated(
   );

   /** assigns open conss and vars if they can be found in blocks */
   SCIP_RETCODE assignAllDependent(
      Seeedpool*       seeedpool
   );

   /** fills out the border of the seeed with the hashmap constoblock if there are still assigned conss and vars */
   SCIP_RETCODE assignBorderFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool
   );

   /** assigns open vars to stairlinking if they can be found in two consecutive blocks, returns true if stairlinkingvars are assigned */
   bool assignCurrentStairlinking(
      Seeedpool*       seeedpool
   );

   /** assigns open conss if they include blockvars, returns true if open conss are assigned */
   bool assignHittingOpenconss(
      Seeedpool* seeedpool
   );

   /** assigns open vars if they can be found in one block, returns true if open vars are assigned */
   bool assignHittingOpenvars(
      Seeedpool* seeedpool
   );

   /** fills out the seeed with the hashmap constoblock if there are still assigned conss and vars */
   SCIP_RETCODE assignSeeedFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool
   );

   /** fills out the seeed with the vector constoblock  */
   SCIP_RETCODE assignSeeedFromConstoblockVector(
      std::vector<int> constoblock,
      int additionalNBlocks,
      Seeedpool* seeedpool
   );

   /** books a constraint to be added to the block constraints of the given block (after calling flushBooked) */
   SCIP_RETCODE bookAsBlockCons(
      int consToBlock,
      int block
   );

   /** books a variable to be added to the block constraints of the given block (after calling flushBooked) */
   SCIP_RETCODE bookAsBlockVar(
      int varToBlock,
      int block
   );

   /** books a constraint to be added to the master constraints (after calling flushBooked) */
   SCIP_RETCODE bookAsMasterCons(
      int consToMaster
   );

   /** books a variable to be added to the master variables (after calling flushBooked) */
   SCIP_RETCODE bookAsMasterVar(
      int varToMaster
   );

   /** books a variable to be added to the linking variables (after calling flushBooked) */
   SCIP_RETCODE bookAsLinkingVar(
      int varToLinking
   );

   /** books a variable to be added to the stairlinking variables of the given block and the following block (after calling flushBooked) */
   SCIP_RETCODE bookAsStairlinkingVar(
      int varToStairlinking,
      int firstBlock
   );

   /** calculates the hashvalue of the seeed for comparing */
   void calcHashvalue(
   );

   /** calculates vector containing constraints not assigned yet */
   void  calcOpenconss(
   );

   /** constructs vector containing variables not assigned yet */
   void  calcOpenvars(
   );

   /** returns whether all constraints are assigned and deletes the vector open cons if all are assigned */
   bool checkAllConsAssigned(
   );

   /** checks the consistency of this seeed */
   bool checkConsistency(
      Seeedpool* seeedpool);

   /** checks the consistency of this seeed */
   bool checkConsistency(
   );

   bool checkVarsAndConssConsistency(
         Seeedpool* seeedpool
   );

   /** assigns the open cons and open vars */
   SCIP_RETCODE completeByConnected(
         Seeedpool* seeedpool
   );

   /** assigns the open cons and open vars */
   SCIP_RETCODE completeGreedily(
         Seeedpool* seeedpool
   );

   /** assigns the open cons and open vars which are implicit assigned */
   SCIP_RETCODE considerImplicits(
         Seeedpool* seeedpool
   );

   /** deletes empty blocks */
   SCIP_RETCODE deleteEmptyBlocks(
   );

   /** deletes an open conss */
   SCIP_RETCODE deleteOpencons(
         int opencons
   );

   /** deletes an open var */
   SCIP_RETCODE deleteOpenvar(
         int openvar
   );

   /** displays the assignments of the conss */
   SCIP_RETCODE displayConss(
   );

   /** displays the relevant information of the seeed */
   SCIP_RETCODE displaySeeed(Seeedpool* seeedpool = NULL
   );

   /** displays the assignments of the vars */
   SCIP_RETCODE displayVars( Seeedpool* seeedpool = NULL
   );

   /** computes the score of the given seeed based on the border, the average density score and the ratio of linking variables*/
   SCIP_Real evaluate(
      Seeedpool* seeedpool
   );

   /** fills out the border of the seeed with the hashmap constoblock */
   SCIP_RETCODE filloutBorderFromConstoblock(
         SCIP_HASHMAP* constoblock,
         int givenNBlocks,
         Seeedpool* seeedpool
   );

   /** fills out a seeed with the hashmap constoblock */
   SCIP_RETCODE filloutSeeedFromConstoblock(
         SCIP_HASHMAP* constoblock,
         int givenNBlocks,
         Seeedpool* seeedpool
   );

   /** finds linking-variables that are actually master-variables. I.e. the variable is adjacent to only master-constraints. */
   SCIP_RETCODE findVarsLinkingToMaster(
       Seeedpool* seeedpool
   );

   /** finds linking-variables that are actually stairlinking-ariables. I.e. the variable is adjacent to constraints in exactly two block. */
   SCIP_RETCODE findVarsLinkingToStairlinking(
       Seeedpool* seeedpool
   );

   /** add all booked constraints to master and delete them from openConss */
   SCIP_RETCODE flushBooked(
   );

   /** returns vector containing master conss */
   const int* getConssForBlock(
         int block
   );

   /** returns the detectorchain */
   DEC_DETECTOR** getDetectorchain(
   );

   /** returns if theis seeed was finished by finishSeeed() method of a detector */
   bool getFinishedByFinisher(
   );



   /** returns the calculated has value of this seeed */
   long getHashValue(
   );

   /** returns the id of the seeed */
   int getID(
   );

   /** returns vector containing linking vars */
   const int* getLinkingvars(
   );

   /** returns array containing master conss */
   const int* getMasterconss(
   );

   /** returns vector containing master vars (every constraint containing a master var is in master )*/
   const int* getMastervars(
   );

   /** return the "maximum white score" (the smaller the better) */
   SCIP_Real getMaxWhiteScore();

   /** returns number of blocks */
   int getNBlocks(
   );

   /** get number of conss */
   int getNConss(
   );

   /** returns size of vector containing master conss */
   int getNConssForBlock(
         int block
   );

   /** returns the number of detectors the seeed is propagated by */
   int getNDetectors(
   );

   /** returns size ofvector containing linking vars */
   int getNLinkingvars(
   );

   /** returns size of array containing master conss */
   int getNMasterconss(
   );

   /** returns size of array containing master vars */
   int getNMastervars(
   );

   /** returns total number of stairlinking vars */
   int getNTotalStairlinkingvars(
      );

   /** returns size of vector containing constraints not assigned yet */
   int getNOpenconss(
   );

   /** returns size of vector containing variables not assigned yet */
   int getNOpenvars(
   );

   /** returns size of vector containing stairlinking vars */
   int getNStairlinkingvars(
         int block
   );

   /** get number of vars */
   int getNVars(
   );

   /** returns size of vector containing vars of a certain block */
   int getNVarsForBlock(
         int block
   );

   /**  returns vector containing constraints not assigned yet */
   const int* getOpenconss(
   );

   /** returns vector containing variables not assigned yet */
   const int* getOpenvars(
   );

   /** returns vector containing stairlinking vars */
   const int* getStairlinkingvars(
         int block
   );

   /** returns vector containing vars of a certain block */
   const int* getVarsForBlock(
         int block
   );

   /** returns whether the cons is a cons of the block */
   bool isConsBlockconsOfBlock(
         int cons, int block
   );

   /** returns whether the cons is a master cons*/
   bool isConsMastercons(
         int cons
   );

   /** return whether the cons is an open conss */
   bool isConsOpencons(
         int cons
   );

   /** returns whether this seeed was propagated by certain detector */
   bool isPropagatedBy(
         DEC_DETECTOR* detectorID
   );

   /** is this seeed trivial (i.e. all constraints in one block, or all conss in border, or all variables linking or mastervars  ) */
   bool isTrivial(
   );

   /** is this seeed complete (i.e. has it least one open constraint or one open variable ) */
   bool isComplete(
   );


/* method to check whether seeed is equal to given other seeed */
   SCIP_RETCODE isEqual(
      Seeed*               otherseeed,          /**< other seeed */
      SCIP_Bool*           isequal,             /**< pointer to store whether seeeds are identical */
      bool                 sortseeeds           /**< should conss and vars be sorted before comparing the seeeds? */
      );


   bool isEqual(
      Seeed* other                              /**< other seeed */
      );


   /** returns whether the var is a var of the block */
   bool isVarBlockvarOfBlock(
         int var,
         int block
   );

   /** returns whether the var is a linking var */
   bool isVarLinkingvar(
         int var
   );

   /** returns whether the var is a master var */
   bool isVarMastervar(
         int var
   );

   /** returns whether the var is an open var */
   bool isVarOpenvar(
         int var
   );

   /** returns whether the var is a stairlinking var */
   bool isVarStairlinkingvar(
      int var
   );

   /** returns whether the var is a stairlinkingvar of the block */
   bool isVarStairlinkingvarOfBlock(
         int var,
         int block
   );

   /** refine seeed: do obvious (considerImplicits()) and some non-obvious assignments assignOpenPartialHittingToMaster() */
   SCIP_RETCODE refineToMaster(
      Seeedpool*       seeedpool
   );

   /** add a constraint to a block */
   SCIP_RETCODE setConsToBlock(
         int consToBlock,
         int block
   );

   /** add a constraint to the master constraints */
   SCIP_RETCODE setConsToMaster(
         int consToMaster
   );

   SCIP_RETCODE setDetectorPropagated(
      DEC_DETECTOR* detectorID
   );

   SCIP_RETCODE setFinishingDetectorPropagated(
      DEC_DETECTOR* detectorID
   );


   /** set if this seeed was finished by finishSeeed() method of a detector */
   void setFinishedByFinisher(
         bool finished
      );

   /** set number of blocks, atm only increasing number of blocks  */
   SCIP_RETCODE setNBlocks(
         int nBlocks
   );

   /** set the id */

   SCIP_RETCODE setID(
            int id
      );

   SCIP_RETCODE setOpenVarsAndConssCalculated(
         bool value
   );

   /** add a variable to a block */
   SCIP_RETCODE setVarToBlock(
         int varToBlock,
         int block
   );


   /** add a variable to the linking variables */
   SCIP_RETCODE setVarToLinking(
         int varToLinking
   );

   /** add a variable to the master variables (every constraint consisting it is in master) */
   SCIP_RETCODE setVarToMaster(
         int varToMaster
   );

   /** add a variable to the stair linking variables */
   SCIP_RETCODE setVarToStairlinking(
         int varToStairLinking, int block1, int block2
   );

   void showScatterPlot(  Seeedpool* seeedpool, SCIP_Bool writeonly = FALSE, const char* filename = NULL, SCIP_Bool draft = FALSE, SCIP_Bool colored = TRUE );

   /** is this seeed a userseeed that should be completed by setting unspecified constraints to master */
   SCIP_Bool shouldCompletedByConsToMaster();

   /** sorts the vars and conss according their numbers */
   void sort(
   );

   const char* getShortCaption();

   /** displays the assignments of the vars */
   SCIP_RETCODE writeScatterPlot(
         Seeedpool* seeedpool,
         const char* filename
   );


private:

   /** assigns open conss (and vars) that hit a block and other open vars (or conss) that are open to border */
   SCIP_RETCODE assignOpenPartialHittingToMaster(
      Seeedpool*       seeedpool
   );

   /** assigns open conss that hit a block and other open vars that are open to border */
   SCIP_RETCODE assignOpenPartialHittingConsToMaster(
      Seeedpool*       seeedpool
   );

   /** assigns open vars that hit a block and other open conss that are open to border */
   SCIP_RETCODE assignOpenPartialHittingVarsToMaster(
      Seeedpool*       seeedpool
   );


 //  bool compare_blocks(int a, int b);

};

} /* namespace gcg */
#endif /* GCG_CLASS_Seeed_H__ */
