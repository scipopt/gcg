/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2013 Operations Research, RWTH Aachen University       */
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

/**@file   dialog_graph.cpp
 * @brief  A dialog to write graph representations of the matrix and read partitions as decompositions.
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "dialog_graph.h"
#include "scip/dialog_default.h"
#include "graph/bipartitegraph.h"
#include "graph/hyperrowcolgraph.h"
#include "graph/hyperrowgraph.h"
#include "graph/hypercolgraph.h"
#include "graph/columngraph.h"
#include "graph/rowgraph.h"
#include "scip_misc.h"
namespace gcg
{

DialogWriteGraph::DialogWriteGraph(
   SCIP*              scip                /**< SCIP data structure */
) : ObjDialog(scip, "write", "write graph to file", TRUE)
{

}

SCIP_DECL_DIALOGEXEC(DialogWriteGraph::scip_exec) {
   SCIP_CALL(SCIPdialogExecMenu(scip, dialog, dialoghdlr, nextdialog));
   return SCIP_OKAY;
}

DialogGraph::DialogGraph(
   SCIP*              scip                /**< SCIP data structure */
) : ObjDialog(scip, "graph", "graph submenu to read and write graph", TRUE)
{

}

SCIP_DECL_DIALOGEXEC(DialogGraph::scip_exec) {
   SCIP_CALL(SCIPdialogExecMenu(scip, dialog, dialoghdlr, nextdialog));
   return SCIP_OKAY;
}
DialogReadPartition::DialogReadPartition(
   SCIP*              scip                /**< SCIP data structure */
) : ObjDialog(scip, "read", "read partition from file", TRUE)
{

}

SCIP_DECL_DIALOGEXEC(DialogReadPartition::scip_exec) {

   SCIP_CALL(SCIPdialogExecMenu(scip, dialog, dialoghdlr, nextdialog));
   return SCIP_OKAY;
}

template<class G>
DialogWriteGraphs<G>::DialogWriteGraphs(
   SCIP*              scip                /**< SCIP data structure */
):  ObjDialog(scip, G(scip, Weights()).name.c_str(), "writes graph of given type", FALSE)
{
   (void)static_cast<Graph*>((G*)0); /* assure we only get descendants of type Graph */
}

template<class G>
SCIP_DECL_DIALOGEXEC(DialogWriteGraphs<G>::scip_exec)
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
      Graph* graph = new G(scip, Weights());
      char* extension;
      extension = filename;
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, extension, TRUE) );
      SCIP_CALL( graph->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
      SCIP_CALL( graph->writeToFile(extension, FALSE) );
      delete graph;
      SCIPdialogMessage(scip, NULL, "graph written to <%s>\n", extension);
   }
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
   return SCIP_OKAY;
}

template<class G>
DialogReadGraphs<G>::DialogReadGraphs(
   SCIP*              scip               /**< SCIP data structure */
): ObjDialog(scip, G(scip, Weights()).name.c_str(), "reads graph of given type", FALSE)
{
   (void)static_cast<Graph*>((G*)0); /* assure we only get descendants of type Graph */
}

template<class G>
SCIP_DECL_DIALOGEXEC(DialogReadGraphs<G>::scip_exec)
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
      Graph* graph = new G(scip, Weights());
      char* extension;
      extension = filename;
      SCIP_CALL( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, extension, TRUE) );
      SCIP_CALL( graph->createFromMatrix(SCIPgetConss(scip), SCIPgetVars(scip), SCIPgetNConss(scip), SCIPgetNVars(scip)) );
      SCIP_CALL( graph->readPartition(extension) );
      //SCIP_CALL( graph->createDecompFromPartition(scip, &decomp) );
      delete graph;
      SCIPdialogMessage(scip, NULL, "decomposition read from <%s>\n", extension);
   }
   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
   return SCIP_OKAY;
}

} /* namespace gcg */

/** include the graph entries for both writing the graph and reading in the partition */
template<class G>
SCIP_RETCODE GCGincludeGraphEntries(
   SCIP*              scip                /**< SCIP data structure */
)
{
   SCIP_DIALOG* graphdialog;
   SCIP_DIALOG* subdialog;

   (void)static_cast<gcg::Graph*>((G*)0); /* assure we only get descendants of type Graph */

   SCIPdialogFindEntry(SCIPgetRootDialog(scip), "graph", &graphdialog);
   assert(graphdialog != NULL);

   SCIPdialogFindEntry(graphdialog, "write", &subdialog);
   assert(subdialog != NULL);
   SCIP_CALL( SCIPincludeObjDialog(scip, subdialog, new gcg::DialogWriteGraphs<G>(scip), true) );

   SCIPdialogFindEntry(graphdialog, "read", &subdialog);
   assert(subdialog != NULL);
   SCIP_CALL( SCIPincludeObjDialog(scip, subdialog, new gcg::DialogReadGraphs<G>(scip), true) );

   return SCIP_OKAY;
}

/** inludes all graph submenu entries */
SCIP_RETCODE GCGincludeDialogsGraph(
   SCIP*              scip                /**< SCIP data structure */
   )
{
   SCIP_DIALOG* dialog;
   SCIP_DIALOG* subdialog;
   dialog = SCIPgetRootDialog(scip);
   SCIP_CALL( SCIPincludeObjDialog(scip, dialog, new gcg::DialogGraph(scip), TRUE) );
   SCIPdialogFindEntry(dialog, "graph", &subdialog);
   assert(subdialog != NULL);
   SCIP_CALL( SCIPincludeObjDialog(scip, subdialog, new gcg::DialogWriteGraph(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjDialog(scip, subdialog, new gcg::DialogReadPartition(scip), TRUE) );

   SCIP_CALL( GCGincludeGraphEntries<gcg::BipartiteGraph>(scip) );
   SCIP_CALL( GCGincludeGraphEntries<gcg::RowGraph>(scip) );
   SCIP_CALL( GCGincludeGraphEntries<gcg::ColumnGraph>(scip) );
   SCIP_CALL( GCGincludeGraphEntries<gcg::HyperrowcolGraph>(scip) );
   SCIP_CALL( GCGincludeGraphEntries<gcg::HyperrowGraph>(scip) );
   SCIP_CALL( GCGincludeGraphEntries<gcg::HypercolGraph>(scip) );
   return SCIP_OKAY;
}
