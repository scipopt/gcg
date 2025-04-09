/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
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

/**@file    gcgvarhistory.c
 * @brief   methods for managing variable history
 * @author  Til Mohr
 */

#include "gcg/gcgvarhistory.h"
#include "gcg/struct_gcgvarhistory.h"

/** free a history buffer */
static
SCIP_RETCODE historybufferFree(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORYBUFFER** buffer              /**< pointer to the history buffer */
   )
{
   GCG_VARHISTORYBUFFER* next;

   assert(scip != NULL);
   assert(buffer != NULL);
   assert(*buffer != NULL);
   assert((*buffer)->nuses == 0);
   assert(((*buffer)->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*buffer)->next != NULL));

   SCIPdebugMessage("Freeing history buffer with %d variables\n", (*buffer)->nvars);

   for( int i = 0; i < (*buffer)->nvars; i++ )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*buffer)->vars[i]) );
   }

   next = (*buffer)->next;

   SCIPfreeBlockMemory(scip, buffer);

   if( next != NULL )
   {
      SCIP_CALL( GCGvarhistoryReleaseBuffer(scip, &next) );
   }

   return SCIP_OKAY;
}

/** get the variable behinde the pointer */
SCIP_RETCODE GCGvarhistoryGetVar(
   GCG_VARHISTORY*        pointer,           /**< pointer to the history */
   SCIP_VAR**             var                /**< pointer to store the variable */
   )
{
   assert(pointer != NULL);
   assert(0 <= pointer->pos);
   assert(var != NULL);
   assert(pointer->buffer != NULL);
   assert(pointer->pos < pointer->buffer->nvars);
   assert((pointer->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= (pointer->buffer->next != NULL));

   if( pointer->pos < 0 )
   {
      *var = NULL;
      return SCIP_OKAY;
   }

   *var = pointer->buffer->vars[pointer->pos];

   return SCIP_OKAY;
}

/** check if there is a next history event */
SCIP_Bool GCGvarhistoryHasNext(
   GCG_VARHISTORY*        pointer            /**< pointer to the history */
   )
{
   assert(pointer != NULL);
   assert(pointer->buffer != NULL);

   assert(pointer->buffer != NULL);
   assert(pointer->pos < pointer->buffer->nvars);
   assert(pointer->buffer->nvars <= GCG_VARHISTORYBUFFER_SIZE);
   assert((pointer->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= (pointer->buffer->next != NULL));

   if( pointer->pos < pointer->buffer->nvars - 1 )
   {
      return TRUE;
   }

   assert(pointer->pos == pointer->buffer->nvars - 1);

   if( pointer->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE )
   {
      return pointer->buffer->next != NULL;
   }

   assert(pointer->buffer->next == NULL);
   return FALSE;
}

/** get the next history event */
SCIP_RETCODE GCGvarhistoryNext(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORY**       pointer            /**< pointer to the history */
   )
{
   GCG_VARHISTORYBUFFER* next;

   assert(scip != NULL);
   assert(pointer != NULL);
   assert(*pointer != NULL);

   assert((*pointer)->buffer != NULL);
   assert((*pointer)->pos < (*pointer)->buffer->nvars);
   assert((*pointer)->buffer->nvars <= GCG_VARHISTORYBUFFER_SIZE);
   assert(((*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*pointer)->buffer->next != NULL));

   if( (*pointer)->pos < (*pointer)->buffer->nvars - 1 && (*pointer)->buffer->nvars > 0)
   {
      SCIPdebugMessage("Advancing history pointer\n");
      if( (*pointer)->pos < 0 )
      {
         (*pointer)->pos = 0;
      }
      else
      {
         (*pointer)->pos += 1;
      }
      return SCIP_OKAY;
   }

   assert((*pointer)->pos == (*pointer)->buffer->nvars - 1);

   next = (*pointer)->buffer->next;
   if( (*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE && next != NULL )
   {
      SCIPdebugMessage("Advancing history pointer to next buffer\n");
      SCIP_CALL( GCGvarhistoryCaptureBuffer(next) );
      SCIP_CALL( GCGvarhistoryReleaseBuffer(scip, &(*pointer)->buffer) );

      (*pointer)->buffer = next;
      (*pointer)->pos = 0;

      assert((*pointer)->buffer->nvars > 0);
      assert((*pointer)->buffer->vars[0] != NULL);

      return SCIP_OKAY;
   }

   assert((*pointer)->buffer->next == NULL);

   (*pointer) = NULL;
   return SCIP_OKAY;
}

/** jump to the latest history event */
SCIP_RETCODE GCGvarhistoryJumpToLatest(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORY**       pointer            /**< pointer to the history */
   )
{
   GCG_VARHISTORYBUFFER* next;
   assert(scip != NULL);
   assert(pointer != NULL);
   assert(*pointer != NULL);
   assert((*pointer)->buffer != NULL);
   assert((*pointer)->pos < (*pointer)->buffer->nvars);
   assert((*pointer)->buffer->nvars <= GCG_VARHISTORYBUFFER_SIZE);
   assert(((*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*pointer)->buffer->next != NULL));

   next = (*pointer)->buffer->next;
   while( next != NULL )
   {
      assert(((*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*pointer)->buffer->next != NULL));
      SCIPdebugMessage("Jumping history pointer to next buffer\n");
      SCIP_CALL( GCGvarhistoryCaptureBuffer(next) );
      SCIP_CALL( GCGvarhistoryReleaseBuffer(scip, &(*pointer)->buffer) );

      (*pointer)->buffer = next;
      next = (*pointer)->buffer->next;

      assert((*pointer)->buffer->nvars > 0);
      assert((*pointer)->buffer->vars[0] != NULL);
   }

   (*pointer)->pos = (*pointer)->buffer->nvars - 1;
   return SCIP_OKAY;
}

/** jump to the latest history event and retrieve all new variables */
SCIP_RETCODE GCGvarhistoryJumpAndRetrieveVars(
   SCIP*                  scip,              /**< SCIP data structure */
   GCG_VARHISTORY**       pointer,           /**< pointer to the history */
   SCIP_VAR***            vars,              /**< pointer to store the variables */
   int*                   nvars              /**< pointer to store the number of variables */
   )
{
   GCG_VARHISTORY weakcopy;
   int curridx;
   int i;
   GCG_VARHISTORYBUFFER* next;

   assert(scip != NULL);
   assert(pointer != NULL);
   assert(*pointer != NULL);
   assert((*pointer)->buffer != NULL);
   assert((*pointer)->pos < (*pointer)->buffer->nvars);
   assert((*pointer)->buffer->nvars <= GCG_VARHISTORYBUFFER_SIZE);
   assert(((*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*pointer)->buffer->next != NULL));
   assert(vars != NULL);
   assert(*vars == NULL);
   assert(nvars != NULL);
   assert(*nvars == 0);

   *vars = NULL;
   *nvars = 0;

   if( (*pointer)->buffer->nvars == 0 )
   {
      return SCIP_OKAY;
   }

   if( (*pointer)->pos < 0 )
   {
      (*pointer)->pos = -1;
   }

   weakcopy.buffer = (*pointer)->buffer;
   weakcopy.pos = (*pointer)->pos;
   do
   {
      *nvars += weakcopy.buffer->nvars - weakcopy.pos - 1;

      if( weakcopy.buffer->next != NULL )
      {
         weakcopy.buffer = weakcopy.buffer->next;
         weakcopy.pos = -1;
      }
      else
      {
         break;
      }
   }
   while( TRUE );

   if( *nvars == 0 )
   {
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, vars, *nvars) );

   curridx = 0;

   do
   {
      assert((*pointer)->buffer->nvars > 0);
      assert((*pointer)->pos < (*pointer)->buffer->nvars);
      assert(((*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*pointer)->buffer->next != NULL));

      for( i = (*pointer)->pos + 1; i < (*pointer)->buffer->nvars; i++ )
      {
         (*vars)[curridx] = (*pointer)->buffer->vars[i];
         curridx += 1;
      }

      next = (*pointer)->buffer->next;
      if( next != NULL )
      {
         SCIP_CALL( GCGvarhistoryCaptureBuffer(next) );
         SCIP_CALL( GCGvarhistoryReleaseBuffer(scip, &(*pointer)->buffer) );

         (*pointer)->buffer = next;
         (*pointer)->pos = -1;
      }
      else
      {
         break;
      }
   }
   while( TRUE );

   assert(curridx == *nvars);

   (*pointer)->pos = (*pointer)->buffer->nvars - 1;
   assert((*pointer)->buffer->next == NULL);
   assert((*pointer)->pos >= 0);

   return SCIP_OKAY;
}

/** create a new history pointer to an empty existing buffer and captures it */
SCIP_RETCODE GCGvarhistoryCreate(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY**       pointer             /**< pointer to the history */
   )
{
   GCG_VARHISTORYBUFFER* buffer;

   assert(pointer != NULL);
   assert(*pointer == NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &buffer) );
   buffer->nvars = 0;
   buffer->next = NULL;
   buffer->nuses = 1;

   SCIP_CALL( SCIPallocBlockMemory(scip, pointer) );
   (*pointer)->buffer = buffer;
   (*pointer)->pos = -1;

   return SCIP_OKAY;
}

/** copy a pointer by creating a new one that points to the same buffer at the same position and capture it */
SCIP_RETCODE GCGvarhistoryCopyReference(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY**       pointer,            /**< pointer to the history */
   GCG_VARHISTORY*        source              /**< source pointer */
   )
{
   assert(pointer != NULL);
   assert(*pointer == NULL);
   assert(source != NULL);
   assert((source->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= (source->buffer->next != NULL));

   SCIP_CALL( GCGvarhistoryCaptureBuffer(source->buffer) );

   SCIP_CALL( SCIPallocBlockMemory(scip, pointer) );
   (*pointer)->buffer = source->buffer;
   (*pointer)->pos = source->pos;

   return SCIP_OKAY;
}

/** release the reference to the buffer and free the history pointer */
SCIP_RETCODE GCGvarhistoryFreeReference(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY**       pointer             /**< pointer to the history */
   )
{
   assert(scip != NULL);
   assert(pointer != NULL);
   assert(*pointer != NULL);
   assert(((*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*pointer)->buffer->next != NULL));

   if( (*pointer)->buffer != NULL )
   {
      SCIP_CALL( GCGvarhistoryReleaseBuffer(scip, &(*pointer)->buffer) );
   }

   SCIPfreeBlockMemory(scip, pointer);

   (*pointer) = NULL;

   return SCIP_OKAY;
}

/** add variable to history */
SCIP_RETCODE GCGvarhistoryAddVar(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORY*        pointer,            /**< pointer to the history */
   SCIP_VAR*              var                 /**< variable */
   )
{
   GCG_VARHISTORYBUFFER* buffer;

   assert(scip != NULL);
   assert(pointer != NULL);
   assert(var != NULL);
   // Check that the pointer is up to date
   assert(pointer->buffer->next == NULL);
   assert(pointer->pos == pointer->buffer->nvars - 1);

   SCIP_CALL( SCIPcaptureVar(scip, var) );

   if( pointer->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE )
   {
      SCIPdebugMessage("Creating new history buffer\n");
      SCIP_CALL( SCIPallocBlockMemory(scip, &buffer) );
      buffer->vars[0] = var;
      buffer->nvars = 1;
      buffer->next = NULL;
      buffer->nuses = 2; // one for the current buffer pointing to the new one and one for the current pointer

      pointer->buffer->next = buffer;

      SCIP_CALL( GCGvarhistoryReleaseBuffer(scip, &pointer->buffer) );
      pointer->buffer = buffer;
      pointer->pos = 0;

      return SCIP_OKAY;
   }
   else
   {
      SCIPdebugMessage("Adding to history buffer\n");
      assert(pointer->buffer->nvars < GCG_VARHISTORYBUFFER_SIZE);
      pointer->buffer->vars[pointer->buffer->nvars] = var;
      pointer->buffer->nvars += 1;
      pointer->pos += 1;
      return SCIP_OKAY;
   }
}

/** capture a reference to the history */
SCIP_RETCODE GCGvarhistoryCaptureBuffer(
   GCG_VARHISTORYBUFFER*  buffer              /**< buffer */
   )
{
   assert(buffer != NULL);
   assert((buffer->nvars == GCG_VARHISTORYBUFFER_SIZE) >= (buffer->next != NULL));
   buffer->nuses += 1;
   return SCIP_OKAY;
}

/** release a reference to the history */
SCIP_RETCODE GCGvarhistoryReleaseBuffer(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORYBUFFER** buffer              /**< buffer */
   )
{
   assert(scip != NULL);
   assert(buffer != NULL);
   assert(*buffer != NULL);
   assert((*buffer)->nuses > 0);
   assert(((*buffer)->nvars == GCG_VARHISTORYBUFFER_SIZE) >= ((*buffer)->next != NULL));

   (*buffer)->nuses -= 1;
   if( (*buffer)->nuses == 0 )
   {
      SCIP_CALL( historybufferFree(scip, buffer) );
   }
   (*buffer) = NULL;
   return SCIP_OKAY;
}
