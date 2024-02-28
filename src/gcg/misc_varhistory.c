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
#include <scip/def.h>
#include <scip/pub_message.h>
#include <scip/scip.h>
#include <scip/type_retcode.h>
#include <scip/type_scip.h>

/** free a history buffer */
static
SCIP_RETCODE historybufferFree(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORYBUFFER** buffer              /**< pointer to the history buffer */
   )
{
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

   GCG_VARHISTORYBUFFER* next = (*buffer)->next;

   SCIPfreeBlockMemory(scip, buffer);

   if( next != NULL )
   {
      SCIP_CALL( GCGvarhistoryReleaseBuffer(scip, &next) );
   }

   return SCIP_OKAY;
}

/** check if there is a next history event */
SCIP_Bool GCGvarhistoryHasNext(
   GCG_VARHISTORYPOINTER* pointer            /**< pointer to the history */
   )
{
   if( pointer == NULL )
   {
      return FALSE;
   }

   assert(pointer->buffer != NULL);
   assert(0 <= pointer->pos);
   assert(pointer->pos < pointer->buffer->nvars);
   assert(pointer->buffer->nvars <= GCG_VARHISTORYBUFFER_SIZE);

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
   GCG_VARHISTORYPOINTER**pointer            /**< pointer to the history */
   )
{
   assert(scip != NULL);
   assert(pointer != NULL);

   if( (*pointer) == NULL )
   {
      return SCIP_OKAY;
   }

   assert((*pointer)->buffer != NULL);
   assert(0 <= (*pointer)->pos);
   assert((*pointer)->pos < (*pointer)->buffer->nvars);
   assert((*pointer)->buffer->nvars <= GCG_VARHISTORYBUFFER_SIZE);

   if( (*pointer)->pos < (*pointer)->buffer->nvars - 1 )
   {
      SCIPdebugMessage("Advancing history pointer\n");
      (*pointer)->pos += 1;
      return SCIP_OKAY;
   }

   assert((*pointer)->pos == (*pointer)->buffer->nvars - 1);

   if( (*pointer)->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE )
   {
      if( (*pointer)->buffer->next != NULL )
      {
         SCIPdebugMessage("Advancing history pointer to next buffer\n");
         SCIP_CALL( GCGvarhistoryCaptureBuffer((*pointer)->buffer->next) );
         SCIP_CALL( GCGvarhistoryReleaseBuffer(NULL, &(*pointer)->buffer) );

         (*pointer)->buffer = (*pointer)->buffer->next;
         (*pointer)->pos = 0;

         assert((*pointer)->buffer->nvars > 0);
         assert((*pointer)->buffer->vars[0] != NULL);

         return SCIP_OKAY;
      }
   }

   assert((*pointer)->buffer->next == NULL);

   (*pointer) = NULL;
   return SCIP_OKAY;
}

/** add variable to history */
SCIP_RETCODE GCGvarhistoryAddVar(
   SCIP*                  scip,               /**< SCIP data structure */
   GCG_VARHISTORYPOINTER* pointer,            /**< pointer to the history */
   SCIP_VAR*              var                 /**< variable */
   )
{
   assert(scip != NULL);
   assert(pointer != NULL);
   assert(var != NULL);

   SCIP_CALL( SCIPcaptureVar(scip, var) );

   if( pointer->buffer == NULL )
   {
      SCIPdebugMessage("Creating new history buffer\n");

      GCG_VARHISTORYBUFFER* buffer;
      SCIP_CALL( SCIPallocBlockMemory(scip, &buffer) );
      buffer->vars[0] = var;
      buffer->nvars = 1;
      buffer->next = NULL;
      buffer->nuses = 1;

      pointer->buffer = buffer;
      pointer->pos = 0;

      return SCIP_OKAY;
   }
   else
   {
      // Check that the pointer is up to date
      assert(pointer->buffer->next == NULL);
      assert(pointer->pos == pointer->buffer->nvars - 1);

      if( pointer->buffer->nvars == GCG_VARHISTORYBUFFER_SIZE )
      {
         SCIPdebugMessage("Creating new history buffer\n");
         GCG_VARHISTORYBUFFER* buffer;
         SCIP_CALL( SCIPallocBlockMemory(scip, &buffer) );
         buffer->vars[0] = var;
         buffer->nvars = 1;
         buffer->next = NULL;
         buffer->nuses = 1;

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
}

/** capture a reference to the history */
SCIP_RETCODE GCGvarhistoryCaptureBuffer(
   GCG_VARHISTORYBUFFER*  buffer              /**< buffer */
   )
{
   assert(buffer != NULL);
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
   (*buffer)->nuses -= 1;
   if( (*buffer)->nuses == 0 )
   {
      SCIP_CALL( historybufferFree(scip, buffer) );
   }
   (*buffer) = NULL;
   return SCIP_OKAY;
}
