/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file    type_nauty.h
 * @ingroup PUBLICCOREAPI
 * @brief   types for nauty automorphism detection
 *
 * @author  Erik Muehmer
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_TYPE_NAUTY_H__
#define GCG_TYPE_NAUTY_H__

#include "nauty/nausparse.h"

#ifdef __cplusplus
struct struct_graph
{
   void add_vertex(int id);
   void add_edge(int v1, int v2);
   unsigned int get_nof_vertices();
   void find_automorphisms(
      void* ptrhook,
      void (*fhook)(void*, unsigned int, const unsigned int*),
      unsigned int searchnodelimit,                                /**< search node limit (requires patched bliss version) */
      unsigned int generatorlimit                                  /**< generator limit (requires patched bliss version or version >=0.76) */
      );

   statsblk stats;
   sparsegraph graph;
   DEFAULTOPTIONS_SPARSEGRAPH(options);
   int* lab;
   int* ptn;
   int* orbits;
};
#endif

#endif /* GCG_TYPE_NAUTY_H__ */
