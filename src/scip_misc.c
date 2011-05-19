/*
 * scip_misc.cpp
 *
 *  Created on: Apr 8, 2010
 *      Author: mbergner
 */

#include "scip_misc.h"
#include "scip/scipdefplugins.h"
#include <string.h>

/**
 * Returns the rhs of an arbitrary SCIP constraint
 * @param scip the scip instance
 * @param cons the constraint
 * @return
 */

consType SCIPconsGetType( SCIP* scip, SCIP_CONS *cons )
{
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return linear;
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch ( SCIPgetTypeSetppc(scip, cons) ) {
      case SCIP_SETPPCTYPE_COVERING:
         return setcovering;
      case SCIP_SETPPCTYPE_PACKING:
         return setpacking;
      case SCIP_SETPPCTYPE_PARTITIONING:
         return setpartitioning;
      default:
         return unknown;
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return logicor;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return knapsack;
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return varbound;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      return sos1;
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      return sos2;
   }
   return unknown;
}

SCIP_Real SCIPgetRhsXXX( SCIP * scip, SCIP_CONS * cons )
{

   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetRhsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch ( SCIPgetTypeSetppc(scip, cons) ) {
      case SCIP_SETPPCTYPE_PARTITIONING: // fall through desired
      case SCIP_SETPPCTYPE_PACKING:
         return 1.0;
      case SCIP_SETPPCTYPE_COVERING:
         return SCIPinfinity(scip);
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPinfinity(scip);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return SCIPgetCapacityKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return SCIPgetRhsVarbound(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPerrorMessage
      ("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPerrorMessage
      ("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else
   {
      SCIPerrorMessage
      ("WARNING: NOT IMPLEMENTED");
   }
   return -SCIPinfinity(scip);
}

/**
 * Returns the lhs of an arbitrary SCIP constraint
 * @param scip the scip instance
 * @param cons the constraint
 * @return
 */
SCIP_Real SCIPgetLhsXXX( SCIP * scip, SCIP_CONS * cons )
{

   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetLhsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      switch ( SCIPgetTypeSetppc(scip, cons) ) {
      case SCIP_SETPPCTYPE_PARTITIONING: // fall through desired
      case SCIP_SETPPCTYPE_COVERING:
         return 1.0;
      case SCIP_SETPPCTYPE_PACKING:
         return -SCIPinfinity(scip);
      }
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return 1.0;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return -SCIPinfinity(scip);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return SCIPgetLhsVarbound(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPerrorMessage
      ("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPerrorMessage
      ("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else
   {
      SCIPerrorMessage
      ("WARNING: NOT IMPLEMENTED");
   }
   return SCIPinfinity(scip);
}

/**
 * Returns the number of variables in an arbitrary SCIP constraint
 * @param scip the scip instance
 * @param cons the constraint
 * @return
 */
int SCIPgetNVarsXXX( SCIP * scip, SCIP_CONS * cons )
{

   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetNVarsLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      return SCIPgetNVarsSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPgetNVarsLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return SCIPgetNVarsKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return 2;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      return SCIPgetNVarsSOS1(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      return SCIPgetNVarsSOS2(scip, cons);
   }
   else
   {
      SCIPerrorMessage("WARNING: NOT IMPLEMENTED <%s>\n", conshdlrname);
      return 0;
   }
}

/**
 * Returns the variable array of an arbitrary SCIP constraint
 * @param scip the scip instance
 * @param cons the constraint
 * @return
 */
SCIP_VAR ** SCIPgetVarsXXX( SCIP * scip, SCIP_CONS * cons )
{

   SCIP_RETCODE retcode;
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   SCIP_VAR ** vars;
   int nvars;
   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      nvars = SCIPgetNVarsLinear(scip, cons);
      retcode = SCIPduplicateMemoryArray(scip, &vars, SCIPgetVarsLinear(scip, cons), nvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      return vars;
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      nvars = SCIPgetNVarsSetppc(scip, cons);
      retcode = SCIPduplicateMemoryArray(scip, &vars, SCIPgetVarsSetppc(scip, cons), nvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      return vars;
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      nvars = SCIPgetNVarsLogicor(scip, cons);
      retcode = SCIPduplicateMemoryArray(scip, &vars, SCIPgetVarsLogicor(scip, cons), nvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      return vars;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      nvars = SCIPgetNVarsKnapsack(scip, cons);
      retcode = SCIPduplicateMemoryArray(scip, &vars, SCIPgetVarsKnapsack(scip, cons), nvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      return vars;
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      retcode = SCIPallocMemoryArray(scip, &vars, 2);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      vars[0] = SCIPgetVarVarbound(scip, cons);
      vars[1] = SCIPgetVbdvarVarbound(scip, cons);
      return vars;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      nvars = SCIPgetNVarsSOS1(scip, cons);
      retcode = SCIPduplicateMemoryArray(scip, &vars, SCIPgetVarsSOS1(scip, cons), nvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      return vars;
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      nvars = SCIPgetNVarsSOS2(scip, cons);
      retcode = SCIPduplicateMemoryArray(scip, &vars, SCIPgetVarsSOS2(scip, cons), nvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      return vars;
   }
   else
   {
      SCIPerrorMessage("WARNING: NOT IMPLEMENTED <%s>\n", conshdlrname);
   }
   return NULL;
}

/**
 * Returns the dual solution value of an arbitrary SCIP constraint
 * @param scip
 * @param cons
 * @return
 */
SCIP_Real SCIPgetDualsolXXX( SCIP * scip, SCIP_CONS * cons )
{

   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      return SCIPgetDualsolLinear(scip, cons);
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      return SCIPgetDualsolSetppc(scip, cons);
   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      return SCIPgetDualsolLogicor(scip, cons);
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {
      return SCIPgetDualsolKnapsack(scip, cons);
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      return SCIPgetDualsolVarbound(scip, cons);
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      SCIPerrorMessage
      ("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      SCIPerrorMessage
      ("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else
   {
      SCIPerrorMessage("WARNING: NOT IMPLEMENTED: ");
      SCIPerrorMessage
      (conshdlrname);
   }
   return 0;
}

/**
 * Returns the value array of an arbitrary SCIP constraint
 * @todo varbound, SOS1 & SOS2 not implemented yet
 * @param scip the scip instance
 * @param cons the constraint
 * @return
 */
SCIP_Real * SCIPgetValsXXX( SCIP * scip, SCIP_CONS * cons )
{

   SCIP_RETCODE retcode;
   SCIP_CONSHDLR* conshdlr;
   const char * conshdlrname;
   int i;
   SCIP_Real * vals;
   assert(scip != NULL);
   assert(cons != NULL);
   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrname = SCIPconshdlrGetName(conshdlr);

   if( strcmp(conshdlrname, "linear") == 0 )
   {
      retcode = SCIPduplicateMemoryArray(scip, &vals, SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons));
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      return vals;
   }
   else if( strcmp(conshdlrname, "setppc") == 0 )
   {
      int nconsvars = SCIPgetNVarsSetppc(scip, cons);
      retcode = SCIPallocMemoryArray(scip, &vals, nconsvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      for( i = 0; i < nconsvars; i++ )
         vals[i] = 1.0;
      return vals;

   }
   else if( strcmp(conshdlrname, "logicor") == 0 )
   {
      int nconsvars = SCIPgetNVarsLogicor(scip, cons);
      retcode = SCIPallocMemoryArray(scip, &vals, nconsvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      for( i = 0; i < nconsvars; i++ )
         vals[i] = 1.0;
      return vals;
   }
   else if( strcmp(conshdlrname, "knapsack") == 0 )
   {

      /* copy Longint array to SCIP_Real array */
      long long int * w = SCIPgetWeightsKnapsack(scip, cons);
      int nconsvars = SCIPgetNVarsKnapsack(scip, cons);
      retcode = SCIPallocMemoryArray(scip, &vals, nconsvars);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      for( i = 0; i < nconsvars; i++ )
         vals[i] = w[i];
      return vals;
   }
   else if( strcmp(conshdlrname, "varbound") == 0 )
   {
      retcode = SCIPallocMemoryArray(scip, &vals, 2);
      if(retcode != SCIP_OKAY)
      {
         SCIPerrorMessage("Memory for array could not be allocated!");
         return NULL;
      }
      vals[0] = 1.0;
      vals[1] = SCIPgetVbdcoefVarbound(scip, cons);
      return vals;
   }
   else if( strcmp(conshdlrname, "SOS1") == 0 )
   {
      /* store constraint */
      SCIPerrorMessage
      ("WARNING: SOS1 NOT IMPLEMENTED\n");
   }
   else if( strcmp(conshdlrname, "SOS2") == 0 )
   {
      /* store constraint */
      SCIPerrorMessage
      ("WARNING: SOS2 NOT IMPLEMENTED\n");
   }
   else
   {
      SCIPerrorMessage
      ("WARNING: UNKNOWN NOT IMPLEMENTED\n");
   }
   return NULL;
}
