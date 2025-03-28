/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   dialog_graph.h
 * @ingroup DIALOGS
 * @brief  A dialog to write graph representations of the matrix and read partitions as decompositions.
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef DIALOG_GRAPH_H_
#define DIALOG_GRAPH_H_

#include "gcg/objdialog.h"
#include "graph/graph.h"

namespace gcg
{

class DialogGraph: public ObjDialog
{
public:
   DialogGraph(
      GCG*               gcg                 /**< GCG data structure */
   );
   virtual ~DialogGraph() {}
   virtual SCIP_DECL_DIALOGEXEC(scip_exec);
};

class DialogWriteGraph: public ObjDialog
{
public:
   DialogWriteGraph(
      GCG*               gcgstruct           /**< GCG data structure */
   );
   virtual ~DialogWriteGraph() {}
   virtual SCIP_DECL_DIALOGEXEC(scip_exec);
};

class DialogReadPartition: public ObjDialog
{
public:
   DialogReadPartition(
      GCG*               gcg                 /**< GCG data structure */
   );
   virtual ~DialogReadPartition() {}
   virtual SCIP_DECL_DIALOGEXEC(scip_exec);
};

template<class T, template <class T1> class G>
class DialogReadGraphs: public ObjDialog
{
private:
   typedef G<T> GRAPH_TYPE;
public:
   DialogReadGraphs(
      GCG*               gcgstruct           /**< GCG data structure */
   );
   virtual ~DialogReadGraphs() {}
   virtual SCIP_DECL_DIALOGEXEC(scip_exec);
};

template<class T, template <class T1> class G>
class DialogWriteGraphs: public ObjDialog
{
private:
   typedef G<T> GRAPH_TYPE;
public:
   DialogWriteGraphs(
      GCG*               gcgstruct           /**< GCG data structure */
   );
   virtual ~DialogWriteGraphs() {}
   virtual SCIP_DECL_DIALOGEXEC(scip_exec);
};
} /* namespace gcg */



#endif /* DIALOG_GRAPH_H_ */
