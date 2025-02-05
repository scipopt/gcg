/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file    misc_varhistory.c
 * @brief   methods for managing variable history
 * @author  Til Mohr
 */

#include "misc_varhistory.h"
#include "struct_varhistory.h"

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
   assert(*var == NULL);
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
   GCG_VARHISTORY weak_copy;
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

   weak_copy.buffer = (*pointer)->buffer;
   weak_copy.pos = (*pointer)->pos;
   do
   {
      *nvars += weak_copy.buffer->nvars - weak_copy.pos - 1;

      if( weak_copy.buffer->next != NULL )
      {
         weak_copy.buffer = weak_copy.buffer->next;
         weak_copy.pos = -1;
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
