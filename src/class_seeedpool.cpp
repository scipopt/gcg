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

/**@file   class_colpool.cpp
 * @brief  class with functions for colpool
 * @author Michael Bastubbe
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg.h"
#include "objscip/objscip.h"
#include "class_seeedpool.h"


#include <exception>

#define SCIP_CALL_EXC(x)   do                                                                                  \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) !=  SCIP_OKAY )                                                \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             throw std::exception();                                                          \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

namespace gcg {

/** constructor */
 Seeedpool::Seeedpool(
    SCIP*             givenScip               /**< SCIP data structure */
    ):scip(givenScip), seeeds(0)
 {

 }


 /** finds decompositions  */
 DEC_DECOMP**       	Seeedpool::findDecompostions(
 ){
	 return NULL;
 }

} /* namespace gcg */
