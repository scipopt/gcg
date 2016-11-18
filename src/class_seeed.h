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





namespace gcg {

class Seeedpool;


class Seeed
{

private:
	SCIP*							         scip;
   int								      id;						         /**< id of the seeed */
   int 								      nBlocks;				            /**< number of blocks the decomposition currently has */
   int 								      nVars;                        /**< number of variables */
   int								      nConss;                       /**< numver of constraints */
   std::vector<int>					   masterConss;			         /**< vector containing indices of master constraints */
   std::vector<int>					   masterVars;				         /**< vector containing indices of master variables */
   std::vector<std::vector<int>>    conssForBlocks; 		         /**< conssForBlocks[k] contains a vector of indices of all constraints assigned to block k */
   std::vector<std::vector<int>>    varsForBlocks; 			      /**< varsForBlocks[k] contains a vector of indices of all variables assigned to block k */
   std::vector<int> 				      linkingVars;			         /**< vector containing indices of linking variables */
   std::vector<std::vector<int>>    stairlinkingVars;		         /**< vector containing indices of staircase linking variables of the blocks */
   std::vector<int> 				      openVars;				         /**< vector containing indices of  variables that are not assigned yet*/
   std::vector<int> 				      openConss;				         /**< vector containing indices of  constraints that are not assigned yet*/
   std::vector<int>                 bookedAsMasterConss;          /**< vector containing indices of  constraints that are not assigned yet but booked as master conss */
   std::vector<std::pair<int, int>> bookedAsBlockConss;           /**< vector containing indices of constraints that are not assigned yet but booked as block conss and the refering block */
   std::vector<bool> 				   propagatedByDetector;	      /**< propagatedByDetector[i] is this seeed propagated by detector i */
   bool 							         openVarsAndConssCalculated;   /**< are the openVars and openCons calculated */
   long                             hashvalue;

   const static int              primes[];
   const static int              nPrimes;

public:
   std::vector<int>                 detectorChain;             /**< vector containing detectors */

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

   /**check-methods */

   /** check the consistency of this seeed */
   bool checkConsistency(
   );

   /** set-methods */

   /** set number of blocks, atm only increasing number of blocks  */
      SCIP_RETCODE setNBlocks(
	   int nBlocks
      );

   /** add a constraint to the master constraints */
   SCIP_RETCODE setConsToMaster(
		   int consToMaster
   );


   /** book a constraint to be added to the master constraints (after calling flushBooked) */
   SCIP_RETCODE bookAsMasterCons(
           int consToMaster
   );

   /** book a constraint to be added to the block constraints of the given block (after calling flushBookes) */
   SCIP_RETCODE bookAsBlockCons(
           int consToBlock,
           int block
   );

   /** add all booked constraints to master and delete them from openConss */
   SCIP_RETCODE flushBooked(
   );


   /** add a variable to the master variables (every constraint consisting it is in master) */
   SCIP_RETCODE setVarToMaster(
		   int varToMaster
   );

   /** add a constraint to a block */
   SCIP_RETCODE setConsToBlock(
		   int consToBlock,
		   int block
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

   /** add a variable to the stair linking variables */
   SCIP_RETCODE setVarToStairlinking(
		   int varToStairLinking, int block1, int block2
   );

   /** add a block, returns the number of the new block */
   int addBlock();

   /** finds linking-variables that are actually master-variables. I.e. the variable is adjacent to only master-constraints. */
   SCIP_RETCODE findVarsLinkingToMaster(
       Seeedpool* seeedpool
   );

   /** finds linking-variables that are actually stairlinking-ariables. I.e. the variable is adjacent to constraints in exactly two block. */
   SCIP_RETCODE findVarsLinkingToStairlinking(
       Seeedpool* seeedpool
   );

   SCIP_RETCODE setDetectorPropagated(
   		   int detectorID
     );

   SCIP_RETCODE setOpenVarsAndConssCalculated(
            bool value
   );


   /** get-methods */

   /** returns array containing master conss */
   const int* getMasterconss(
   );

   /** returns size of array containing master conss */
   int getNMasterconss(
   );

   /** returns vector containing master vars (every constraint containing a master var is in master )*/
   const int* getMastervars(
   );

   /** returns size of array containing master vars */
   int getNMastervars(
   );

   /** returns number of blocks */
   int getNBlocks(
   );

   /** returns vector containing master conss */
   const int* getConssForBlock(
		   int block
   );

   /** returns vector containing vars of a certain block */
   const int* getVarsForBlock(
		   int block
   );

   /** returns vector containing linking vars */
   const int* getLinkingvars(
   );

   /** returns vector containing stairlinking vars */
   const int* getStairlinkingvars(
         int block
   );

   /** returns vector containing variables not assigned yet */
   const int* getOpenvars(
   );

   /**  returns vector containing constraints not assigned yet */
   const int* getOpenconss(
   );

   bool alreadyAssignedConssToBlocks();

   /** returns size of vector containing master conss */
   int getNConssForBlock(
		   int block
   );

   /** returns size of vector containing vars of a certain block */
   int getNVarsForBlock(
		   int block
   );

   /** returns size ofvector containing linking vars */
   int getNLinkingvars(
   );

   /** returns size of vector containing stairlinking vars */
   int getNStairlinkingvars(
      int block
   );

   /** returns the calculated has value of this seeed */
   long getHashValue(
   );

   /** calculates vector containing constraints not assigned yet */
   void  calcOpenconss();

   /** constructs vector containing variables not assigned yet */
   void  calcOpenvars();

   /** sorts the vars and conss according their numbers */
   void sort(
   );

   /** returns size of vector containing variables not assigned yet */
   int getNOpenvars(
   );

   /** returns size of vector containing constraints not assigned yet */
   int getNOpenconss(
   );

   /** returns whether this seeed was propagated by certain detector */
   bool isPropagatedBy(
		   int detectorID
   );

   /** returns if the open vars and conss are calculated */
   bool areOpenVarsAndConssCalculated(
   );

   /** assigns the open cons and open vars */
   SCIP_RETCODE completeGreedily(
         Seeedpool* seeedpool
   );

   /** assigns the open cons and open vars which are implicit assigned */
   SCIP_RETCODE considerImplicits(
         Seeedpool* seeedpool
   );

   /** assigns open conss if they includes blockvars, returns true if open conss are assigned */
   bool assignHittingOpenconss(
         Seeedpool* seeedpool
   );

   /** assigns open vars if they can be found in one block, returns true if open vars are assigned */
   bool assignHittingOpenvars(
         Seeedpool* seeedpool
   );

   /** assigns open vars to stairlinking if they can be found in two consecutive  blocks, returns true if stairlinkingvars are assigned */
   bool assignCurrentStairlinking(
         Seeedpool*       seeedpool
   );

   /** assigns open conss and vars if they can be found in blocks */
   SCIP_RETCODE assignAllDependent(
         Seeedpool*       seeedpool
   );

   /** returns whether the var is a linking var */
   bool isVarLinkingvar(
         int var
   );

   /** returns whether the var is a var of the block */
   bool isVarBlockvarOfBlock(
         int var, int block
   );

   /** returns whether the var is a stairlinkingvar of the block */
   bool isVarStairlinkingvarOfBlock(
         int var, int block
   );

   /** returns whether the var is a master var */
   bool isVarMastervar(
         int var
   );

   /** returns whether the var is an open var */
   bool isVarOpenvar(
         int var
   );

   /** returns whether the cons is a master cons*/
   bool isConsMastercons(
         int cons
   );

   /** return whether the cons is an open conss */
   bool isConsOpencons(
         int cons
   );

   /** returns whether the cons is a cons of the block */
   bool isConsBlockconsOfBlock(
         int cons, int block
   );

   /** returns index of the Openvar in the vector */
   int getIndexOfOpenvar(
         int var
   );

   /** returns whether all cons are assigned and deletes the vector open cons if all are assigned */
   bool checkAllConsAssigned(
   );

   /** returns the detectorchain */
   int* getDetectorchain(
   );

   /** returns the number of detectors the seeed is propagated by */
   int getNDetectors(
   );

   /** returns the id of the seeed */
   int getID(
   );

   /** get number of conss */
   int getNConss(
   );

   /** get number of vars */
   int getNVars(
   );

   /** fills out a seeed with the hashmap constoblock */
   SCIP_RETCODE filloutSeeedFromConstoblock(
         SCIP_HASHMAP* constoblock,
         int givenNBlocks,
         Seeedpool* seeedpool
   );

   /** fills out the border of the seeed with the hashmap constoblock */
   SCIP_RETCODE filloutBorderFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool
   );

   /** fills out the seeed with the hashmap constoblock if there are still assigned conss and vars */
   SCIP_RETCODE assignSeeedFromConstoblock(
         SCIP_HASHMAP* constoblock,
         int givenNBlocks,
         Seeedpool* seeedpool
   );

   /** fills out the vorder of the seeed with the hashmap constoblock if there are still assigned conss and vars */
   SCIP_RETCODE assignBorderFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool
   );

   /** deletes an open var */
   SCIP_RETCODE deleteOpenvar(
         int openvar
   );

   /** deletes an open conss */
   SCIP_RETCODE deleteOpencons(
         int opencons
   );

   /** deletes empty blocks */
   SCIP_RETCODE deleteEmptyBlocks(
   );

   /** fills out the vectors for free conss and free vars */
   SCIP_RETCODE identifyFreeConssAndVars(
         Seeedpool* seeedpool
   );

   /** calculates the hashvalue of the seeed for comparing */
   void calcHashvalue(
   );

   bool checkVarsAndConssConsistency(
         Seeedpool* seeedpool
   );

   /** displays the relevant information of the seeed */
   SCIP_RETCODE displaySeeed(
   );

   /** displays the assignments of the conss */
   SCIP_RETCODE displayConss(
   );

   /** displays the assignments of the vars */
   SCIP_RETCODE displayVars(
   );




private:

 //  bool compare_blocks(int a, int b);

};

} /* namespace gcg */
#endif /* GCG_CLASS_Seeed_H__ */
