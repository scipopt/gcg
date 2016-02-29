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

/**@file   class_seeedpool.h
 * @brief  class with functions for seeed pool where a seeed is a (potentially incomplete) description of a decomposition (not to confuse with the band from German capital)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_SEEEDPOOL_H__
#define GCG_CLASS_SEEEDPOOL_H__

#include "objscip/objscip.h"
#include <vector>
#include <tr1/unordered_map> //c++ hashmap

#include "gcg.h"

#include "class_seeed.h"


namespace gcg {

//typedef boost::shared_ptr<Seeed> SeeedPtr;
typedef Seeed* SeeedPtr;


class Seeedpool
{   /*lint -esym(1712,Seeedpool)*/

private:
   SCIP*                 						scip;              	/**< SCIP data structure */
   std::vector<SeeedPtr> 						seeeds;				/**< vector of current (open) seeeds */

   std::vector<std::vector<int>> 				varsForConss; 		/** stores for every constraint the indices of variables that are contained in the constraint */
   std::vector<std::vector<int>> 				conssForVars; 		/** stores for every variable the indices of constraints containing this variable */

   std::vector<SCIP_CONS*> 						consToScipCons;	/** stores the corresponding scip constraints pointer */
   std::vector<SCIP_VAR*> 						varToScipVar;		/** stores the corresponding scip variable pointer */
   std::vector<DEC_DETECTOR*> 					detectorToScipDetector; /** stores the corresponding SCIP detector pinter */
   std::tr1::unordered_map<SCIP_CONS*, int> 	scipConsToIndex;	/** maps SCIP_CONS* to the corresponding index */
   std::tr1::unordered_map<SCIP_VAR*, int>  	scipVarToIndex;		/** maps SCIP_VAR* to the corresponding index */
   std::tr1::unordered_map<DEC_DETECTOR*, int>  scipDetectorToIndex;		/** maps SCIP_VAR* to the corresponding index */

   int 											nVars;
   int 											nConss;
   int											nDetectors;


public:

   /** constructor */
   Seeedpool(
      SCIP*             scip, /**< SCIP data structure */
	  const char*	conshdlrName
      );

   ~Seeedpool();

   /** finds decompositions  */
   DEC_DECOMP**       	findDecompostions(
   );

   /** access coefficient matrix constraint-wise */
   std::vector<int> const & getVarsForCons(int cons);

   /** access coefficient matrix variable-wise */
   std::vector<int> const & getConssForVar(int varid);


};

} /* namespace gcg */
#endif /* GCG_CLASS_SEEEDPOOL_H__ */
