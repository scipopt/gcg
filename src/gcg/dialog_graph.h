/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
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
