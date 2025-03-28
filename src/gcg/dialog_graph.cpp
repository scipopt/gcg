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

/**@file   dialog_graph.cpp
 * @brief  A dialog to write graph representations of the matrix and read partitions as decompositions.
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/dialog_graph.h"
#include "scip/dialog_default.h"
#include "graph/bipartitegraph.h"
#include "graph/hyperrowcolgraph.h"
#include "graph/hyperrowgraph.h"
#include "graph/hypercolgraph.h"
#include "graph/columngraph.h"
#include "graph/rowgraph.h"
#include "gcg/scip_misc.h"
#include "gcg/cons_decomp.h"
#include "graph/graph_tclique.h"

namespace gcg
{

DialogWriteGraph::DialogWriteGraph(
   GCG*               gcgstruct                 /**< GCG data structure */
) : ObjDialog(gcgstruct, "write", "write graph to file", TRUE)
{

}

SCIP_DECL_DIALOGEXEC(DialogWriteGraph::scip_exec) {
   SCIP_CALL(SCIPdialogExecMenu(scip, dialog, dialoghdlr, nextdialog));
   return SCIP_OKAY;
}

DialogGraph::DialogGraph(
   GCG*               gcgstruct                 /**< GCG data structure */
) : ObjDialog(gcgstruct, "graph", "graph submenu to read and write graph", TRUE)
{

}

SCIP_DECL_DIALOGEXEC(DialogGraph::scip_exec) {
   SCIP_CALL(SCIPdialogExecMenu(scip, dialog, dialoghdlr, nextdialog));
   return SCIP_OKAY;
}
DialogReadPartition::DialogReadPartition(
   GCG*               gcgstruct                 /**< GCG data structure */
) : ObjDialog(gcgstruct, "read", "read partition from file", TRUE)
{

}

SCIP_DECL_DIALOGEXEC(DialogReadPartition::scip_exec) {

   SCIP_CALL(SCIPdialogExecMenu(scip, dialog, dialoghdlr, nextdialog));
   return SCIP_OKAY;
}

template<class T, template <class T1> class G>
DialogWriteGraphs<T,G>::DialogWriteGraphs(
   GCG*               gcgstruct                 /**< GCG data structure */
):  ObjDialog(gcgstruct, G<T>(gcgstruct, Weights()).name.c_str(), "writes graph of given type", FALSE)
{
   (void)static_cast<MatrixGraph<T>*>((G<T>*)0); /* assure we only get descendants of type Graph */
}


template<class T,template <class T1> class G>
SCIP_RETCODE DialogWriteGraphs<T, G>::scip_exec(SCIP* scip, SCIP_DIALOG* dialog, SCIP_DIALOGHDLR* dialoghdlr, SCIP_DIALOG** nextdialog)
{
   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      SCIPdialogMessage(scip, NULL, "No problem exists, read in a problem first.\n");
      return SCIP_OKAY;
   }

   char* filename;
   SCIP_Bool endoffile;

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {

      char* extension;
      int fd;
      FILE* file;

      extension = filename;

      file = fopen(filename, "wx");
      if( file == NULL )
         return SCIP_FILECREATEERROR;

      fd = fileno(file);
      if( fd == -1 )
               return SCIP_FILECREATEERROR;

      MatrixGraph<T>* graph = new G<T>(gcg, Weights());
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, extension, TRUE) );
      SCIP_CALL( graph->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
      SCIP_CALL( graph->writeToFile(fd, FALSE) );
      delete graph;
      SCIPdialogMessage(scip, NULL, "graph written to <%s>\n", extension);
   }
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
   return SCIP_OKAY;
}

template<class T, template <class T1> class G>
DialogReadGraphs<T,G>::DialogReadGraphs(
   GCG*               gcgstruct                /**< GCG data structure */
): ObjDialog(gcgstruct, G<T>(gcgstruct, Weights()).name.c_str(), "reads graph of given type", FALSE)
{
   (void)static_cast<MatrixGraph<T>*>((G<T>*)0); /* assure we only get descendants of type Graph */
}

template<class T, template <class T1> class G>
SCIP_RETCODE DialogReadGraphs<T, G>::scip_exec(SCIP* scip, SCIP_DIALOG* dialog, SCIP_DIALOGHDLR* dialoghdlr, SCIP_DIALOG** nextdialog)
{
   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
      SCIPdialogMessage(scip, NULL, "No problem exists, read in a problem first.\n");
      return SCIP_OKAY;
   }

   char* filename;
   SCIP_Bool endoffile;

   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ", &filename, &endoffile) );
   if( endoffile )
   {
      *nextdialog = NULL;
      return SCIP_OKAY;
   }
   if( filename[0] != '\0' )
   {
      MatrixGraph<T>* graph = new G<T>(gcg, Weights());
      char* extension;
      extension = filename;
      GCG_DECOMP* decomp;
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, extension, TRUE) );
      SCIP_CALL( graph->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
      SCIP_CALL( graph->readPartition(extension) );
      SCIP_CALL( graph->createDecompFromPartition(&decomp) );
      delete graph;

      SCIP_CALL( GCGconshdlrDecompAddPreexistingDecomp(gcg, decomp) );
      GCGdecompFree(gcg, &decomp);
      SCIPdialogMessage(scip, NULL, "decomposition read from <%s>\n", extension);
   }
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
   return SCIP_OKAY;
}

} /* namespace gcg */

/** include the graph entries for both writing the graph and reading in the partition */
template<class T, template <class T1> class G>
SCIP_RETCODE GCGincludeGraphEntries(
   GCG*               gcg                 /**< GCG data structure */
)
{
   SCIP_DIALOG* graphdialog;
   SCIP_DIALOG* subdialog;
   SCIP* origprob = GCGgetOrigprob(gcg);

   (void)static_cast<gcg::MatrixGraph<T>*>((G<T>*)0); /* assure we only get descendants of type Graph */

   (void) SCIPdialogFindEntry(SCIPgetRootDialog(origprob), "graph", &graphdialog);
   assert(graphdialog != NULL);

   (void) SCIPdialogFindEntry(graphdialog, "write", &subdialog);
   assert(subdialog != NULL);
   SCIP_CALL( GCGincludeObjDialog(gcg, subdialog, new gcg::DialogWriteGraphs<T,G>(gcg), true) );

   (void) SCIPdialogFindEntry(graphdialog, "read", &subdialog);
   assert(subdialog != NULL);
   SCIP_CALL( GCGincludeObjDialog(gcg, subdialog, new gcg::DialogReadGraphs<T,G >(gcg), true) );

   return SCIP_OKAY;
}

/** inludes all graph submenu entries */
extern "C"
SCIP_RETCODE GCGincludeDialogsGraph(
   GCG*               gcg                 /**< GCG data structure */
   )
{
   SCIP_DIALOG* dialog;
   SCIP_DIALOG* subdialog;
   SCIP* origprob = GCGgetOrigprob(gcg);
   dialog = SCIPgetRootDialog(origprob);
   SCIP_CALL( GCGincludeObjDialog(gcg, dialog, new gcg::DialogGraph(gcg), TRUE) );
   (void) SCIPdialogFindEntry(dialog, "graph", &subdialog);
   assert(subdialog != NULL);
   SCIP_CALL( GCGincludeObjDialog(gcg, subdialog, new gcg::DialogWriteGraph(gcg), TRUE) );
   SCIP_CALL( GCGincludeObjDialog(gcg, subdialog, new gcg::DialogReadPartition(gcg), TRUE) );

   SCIP_CALL( (GCGincludeGraphEntries<gcg::GraphTclique,gcg::RowGraph>(gcg)) );
#ifdef SCIP_DISABLED_CODE
   /*SCIP_CALL*/( GCGincludeGraphEntries<gcg::GraphTclique,gcg::BipartiteGraph>(origprob) );
   /*SCIP_CALL*/( GCGincludeGraphEntries<gcg::GraphTclique,gcg::ColumnGraph>(origprob) );
   /*SCIP_CALL*/( GCGincludeGraphEntries<gcg::GraphTclique,gcg::HyperrowcolGraph>(origprob) );
   /*SCIP_CALL*/( GCGincludeGraphEntries<gcg::GraphTclique,gcg::HyperrowGraph>(origprob) );
   /*SCIP_CALL*/( GCGincludeGraphEntries<gcg::GraphTclique,gcg::HypercolGraph>(origprob) );
#endif
   return SCIP_OKAY;
}
