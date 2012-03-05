/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: dec_stairheur.c,v 1.24 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   dec_stairheur.c
 * @ingroup DETECTORS
 * @brief  stairheur presolver
 * @author Mathias Luers
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>

#include "dec_stairheur.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/struct_var.h"

#define DEC_DETECTORNAME      "stairheur"       /**< name of the detector */
#define DEC_PRIORITY          -100              /**< priority of the detector */

/* Default parameter settings*/
#define DEFAULT_BLOCKS                    2     /**< number of blocks */
//#define DEFAULT_CONSWEIGHT                5     /**< weight for constraint hyperedges */
//#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
//#define DEFAULT_TIDY                      TRUE  /**< whether to clean up afterwards */
//#define DEFAULT_DUMMYNODES	              0.2   /**< percentage of dummy vertices*/

#define DEFAULT_MAXBLOCKS                 20    /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                 2     /**< value for the minimum number of blocks to be considered */

//#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
//#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
//#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_PRIORITY                  DEC_PRIORITY

#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)

#define GP_NAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_%d.gp", (name), (blocks), (cons), (dummy)


/*include the code for doubly linked lists */
#include <stdio.h>
#include <stdlib.h>

/** TODO:
 * currently, all vars from the first column where a linking var appears until the end of the block are considered as linking vars, although there might be empty columns. This could be changed so that these empty columns are considered as subscipvars and not linking vars.
 *
 * In some cases a block can consist of linking vars exclusively. This makes no real sense.
 *
 * For some instances the assertion regarding the consistency of the arrays ibegin and jbegin fails
 * */

struct node {
 void* data;
 struct node* prev;
 struct node* next;
};
typedef struct node NODE;


struct list {
   /* a doubly linked list */
   NODE* nil;
   int size;
};

/** a doubly linked list */
typedef struct list LIST;

struct iterator {
   NODE* node;
   LIST* list;
};

/** an iterator for traversing lists */
typedef struct iterator ITERATOR;

/*lists*/
LIST* SCIPlistCreate(SCIP* scip);
LIST* SCIPlistCreateInt(SCIP* scip, int from, int to);
LIST* SCIPlistCopyShallow(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistPushFront(SCIP* scip, LIST* list, void* data);
SCIP_Bool SCIPlistPopFront(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistPushBack(SCIP* scip, LIST* list, void* data);
SCIP_Bool SCIPlistPopBack(SCIP* scip, LIST* list);
NODE* SCIPlistInsert(SCIP* scip, ITERATOR* it, void* data);
SCIP_Bool SCIPlistAssign(ITERATOR* it, void* data);
SCIP_Bool SCIPlistAssignList(LIST* list1, LIST* list2);
SCIP_Bool SCIPlistIsEmpty(LIST* list);
SCIP_Bool SCIPlistErase(SCIP* scip, ITERATOR* it);
SCIP_Bool SCIPlistDelete(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistDeleteData(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistDeleteNested(SCIP* scip, LIST* list);
ITERATOR SCIPlistFind(ITERATOR first, ITERATOR last, void* value, SCIP_Bool (*comp_func)(void* a, void* b));
SCIP_Bool SCIPlistMove(ITERATOR position, ITERATOR it);
void SCIPlistMoveFront(ITERATOR position);
void SCIPlistMoveBack(ITERATOR position);
SCIP_RETCODE SCIPlistRearrange(SCIP* scip, LIST* list, LIST* order);
SCIP_Bool SCIPlistForeach(LIST* list, int(*func)(void*));
void SCIPlistPrint(LIST* list, int(*printfunc)(void*));

/*node*/
NODE* SCIPnodeCreate(SCIP* scip, void* data);

/*iteratos*/
ITERATOR SCIPiteratorBegin(LIST* list);
ITERATOR SCIPiteratorEnd(LIST* list);
SCIP_Bool SCIPiteratorNext(ITERATOR* it);
SCIP_Bool SCIPiteratorPrev(ITERATOR* it);
void SCIPiteratorEquals(ITERATOR* it1, ITERATOR* it2);
SCIP_Bool SCIPiteratorIsEqual(ITERATOR it1, ITERATOR it2);

/*print, compare callback*/
int printstring(void *s);
int printint(void *i);
static int compare (const void * a, const void * b);
static SCIP_Bool compare_int(void* a, void * b);

/** Creates an empty lists and returns a pointer to that list */
LIST* SCIPlistCreate(SCIP* scip)
{
   LIST* list;
   NODE* node;
   if( SCIPallocBlockMemory(scip, &list) == SCIP_NOMEMORY) return NULL;
   node = SCIPnodeCreate(scip, NULL);
   list->nil = node;
   //list->nil->next contains a pointer to the last node
   list->nil->next = list->nil;
   //list->nil->prev contains a pointer to the first node
   list->nil->prev = list->nil;
   list->nil->data = node->data;
   list->size = 0;
   return list;
}

/** Creates a list with integers running from 'from' to 'to'. */
LIST* SCIPlistCreateInt(SCIP* scip, int from, int to)
{
   LIST* list;
   int* data;
   int i;
   list = SCIPlistCreate(scip);
   for(i = from; i <= to; ++i)
   {
      SCIPallocMemory(scip, &data);
      *data=i;
      SCIPlistPushBack(scip, list, (void*) data);
   }
   return list;
}

/** Creates a copy of 'list', but does not copy the data stored in list->data. */
LIST* SCIPlistCopyShallow(SCIP* scip, LIST* list)
{
   LIST* copy;
   ITERATOR it;
   //case: list exists
   if(list)
   {
      copy = SCIPlistCreate(scip);
      for(it = SCIPiteratorBegin(list); ! SCIPiteratorIsEqual(it, SCIPiteratorEnd(list)); SCIPiteratorNext(&it))
      {
         SCIPlistPushBack(scip, copy, it.node->data);
      }
      return copy;
   }
   //case: list does not exist
   else
   {
      return NULL;
   }
}

/** Inserts a new node with a pointer to 'data' at the front of the list. */
SCIP_Bool SCIPlistPushFront(SCIP* scip, LIST* list, void* data)
{
   /** adds a node at the front of the list */
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return FALSE;
   //it points to the beginning of list
   it = SCIPiteratorBegin(list);
   SCIPlistInsert(scip, &it, data);
   return TRUE;
}

/** Removes the node at the front of the list. */
SCIP_Bool SCIPlistPopFront(SCIP* scip, LIST* list)
{
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (SCIPlistIsEmpty(list)) ) return FALSE;

   it = SCIPiteratorBegin(list);
   SCIPlistErase(scip, &it);
   return TRUE;
}

/** Inserts a new node with a pointer to 'data' at the back of the list. */
SCIP_Bool SCIPlistPushBack(SCIP* scip, LIST* list, void* data)
{
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return FALSE;
   //it points to the end of list
   it = SCIPiteratorEnd(list);
   SCIPlistInsert(scip, &it, data);
   return TRUE;
}

/** Removes the node at the back of the list. */
SCIP_Bool SCIPlistPopBack(SCIP* scip, LIST* list)
{
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (SCIPlistIsEmpty(list)) ) return FALSE;

   it = SCIPiteratorEnd(list);
   SCIPiteratorPrev(&it);
   SCIPlistErase(scip, &it);
   return TRUE;
}

/** Creates a new node before 'it->node' and returns the new node. */
NODE* SCIPlistInsert(SCIP* scip, ITERATOR* it, void* data)
{
   NODE* newnode;
   NODE* successor;
   NODE* predecessor;
   //case: 'it' points to no node or list
   if(it->node == NULL && it->list == NULL) return NULL;
   newnode=SCIPnodeCreate(scip, data);
   successor = it->node;
   predecessor = it->node->prev;
   //adjust all pointers
   //note: due to sentinel node no distinction of cases is necessary(begin, middle, end)
   predecessor->next = newnode;
   newnode->prev = predecessor;
   newnode->next = successor;
   successor->prev = newnode;

   ++it->list->size;
   return newnode;
}

/** Assigns new data to the node the iterator 'it' points to. */
SCIP_Bool SCIPlistAssign(ITERATOR* it, void* data)
{
   if(it && it->node)
   {
      it->node->data=data;
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** If target and origin have the same size, this function assigns the data pointers of origin to target. */
SCIP_Bool SCIPlistAssignList(LIST* target, LIST* origin)
{
   ITERATOR it1;
   ITERATOR it2;
   if(target && origin && target->size == origin->size)
   {
      for(it1 = SCIPiteratorBegin(target), it2 = SCIPiteratorBegin(origin); ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(target)); SCIPiteratorNext(&it1), SCIPiteratorNext(&it2))
      {
         SCIPlistAssign(&it1, it2.node->data);
      }
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** Creates a node with a void pointer to data and returns the node. */
NODE* SCIPnodeCreate(SCIP* scip, void* data)
{
   NODE* node;
   if( SCIPallocBlockMemory(scip, &node) == SCIP_NOMEMORY) return NULL;
   node->data = data;
   node->prev = NULL;
   node->next = NULL;
   return node;
}

/** Removes the node the iterator 'it' points to from the list.
 *
 * Returns TRUE if 'node' is in 'list'. Returns FALSE if 'node' is not in 'list'.
 *
 * The memory the data pointer points to is not deallocated.*/
SCIP_Bool SCIPlistErase(SCIP* scip, ITERATOR* it)
{
   NODE* successor;
   NODE* predecessor;
   //case: iterator, node or list do not exist or empty list
   if( (! it) || (! it->node) || (! it->list) || SCIPlistIsEmpty(it->list))
   {
      return FALSE;
   }
   successor = it->node->next;
   predecessor = it->node->prev;
   predecessor->next = successor;
   successor->prev = predecessor;
   SCIPfreeBlockMemory(scip, &it->node);
   it->node = successor;
   --it->list->size;
   return TRUE;
}

/** Deletes the entire list, but not the data stored in the list.
 *
 *  For deallocating memory of the list data, call 'list_deleta_data' before. */
SCIP_Bool SCIPlistDelete(SCIP* scip, LIST* list)
{
   if(! list)
   {
      return FALSE;
   }
   else
   {
      while(! SCIPlistIsEmpty(list))
      {
         SCIPlistPopBack(scip, list);
      }
      //deallocate sentinel node
//      SCIPfreeMemory(scip, &list->nil->data);
      SCIPfreeBlockMemory(scip, &list->nil);
      SCIPfreeBlockMemory(scip, &list);
      list = NULL;
      return TRUE;
   }
}

/** Deallocates the memory the data pointers of the list points to. */
SCIP_Bool SCIPlistDeleteData(SCIP* scip, LIST* list)
{
   ITERATOR it;
   if(! list)
   {
      return FALSE;
   }
   else
   {
      for( it = SCIPiteratorBegin(list); ! ( SCIPiteratorIsEqual(it, SCIPiteratorEnd(list)) ); SCIPiteratorNext(&it))
      {
         SCIPfreeMemory(scip, &it.node->data);
      }
      return TRUE;
   }
}

/** Deallocates all memory for a nested list including the data. */
SCIP_Bool SCIPlistDeleteNested(SCIP* scip, LIST* list)
{
   ITERATOR it1;
   if(! list)
   {
      return FALSE;
   }
   else
   {
      for(it1 = SCIPiteratorBegin(list);  ! ( SCIPiteratorIsEqual(it1, SCIPiteratorEnd(list)) ); SCIPiteratorNext(&it1))
         {
            SCIPlistDeleteData(scip, it1.node->data);
            SCIPlistDelete(scip, it1.node->data);
         }
      SCIPlistDelete(scip, list);
      return TRUE;
   }
}

/** Returns TRUE iff list exists and is empty. */
SCIP_Bool SCIPlistIsEmpty(LIST* list)
{
   //list is empty if the next pointer of the sentinel points to itself
   if( list && list->nil->next == list->nil)
   {
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** Runs the function 'func' for every node of the list. */
SCIP_Bool SCIPlistForeach(LIST* list, int(*func)(void*))
{
   ITERATOR it;
   for(it = SCIPiteratorBegin(list); ! ( SCIPiteratorIsEqual(it, SCIPiteratorEnd(list)) ); SCIPiteratorNext(&it))
   {
      if(func(it.node->data) != 0)
      {
         return FALSE;
      }
   }
   return TRUE;
}

/** Prints the contents of the list. */
void SCIPlistPrint(
      LIST* list,             /**< list to print */
      int(*printfunc)(void*)  /**< function to print contents of the list */
      )
{
   printf("( ");
   SCIPlistForeach(list, printfunc);
   printf(")\n");
}

/** Searches in the interval [first, last) for 'value'.
 * @param first position to start the search
 * @param last position before the search is stopped
 * @param value value which shall be found
 * @param comp_func a predicate function
 * @return Iterator pointing to the position of the first appearance of 'value'. Iterator points to the position after the end of the list, if 'value' is not found. */
ITERATOR SCIPlistFind(ITERATOR first, ITERATOR last, void* value, SCIP_Bool (*comp_func)(void* a, void* b))
{
   for (;
        //first != last AND end of list is not reached yet
        ! (SCIPiteratorIsEqual(first, last)) && ! (SCIPiteratorIsEqual(first, SCIPiteratorEnd(first.list)));
        SCIPiteratorNext(&first)
        )
   {
      if ( comp_func(value, first.node->data) )
      {
         break;
      }
   }
   return first;
}

/** Moves the element at 'position' to the location in front of element 'it'
 *
 * Both elements must be in the same list. */
SCIP_Bool SCIPlistMove(ITERATOR position, ITERATOR it)
{
   NODE* old_successor;
   NODE* old_predecessor;
   NODE* new_successor;
   NODE* new_predecessor;
   if(position.node && it.node && position.list == it.list)
   {
      //case 'it' and 'position' point to the same node or 'it' points to the succeeding node of 'position': no moving
      if(position.node == it.node || position.node->next == it.node)
      {
         return TRUE;
      }
      //all other cases
      //adjust pointers at the old position
      old_successor = position.node->next;
      old_predecessor = position.node->prev;
      old_predecessor->next = old_successor;
      old_successor->prev = old_predecessor;

      //adjust pointers at the new position
      new_successor = it.node;
      new_predecessor = it.node->prev;
      new_predecessor->next = position.node;
      position.node->prev = new_predecessor;
      position.node->next = new_successor;
      new_successor->prev = position.node;

      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** Moves the element at 'position' to the front of the list. */
void SCIPlistMoveFront(ITERATOR position)
{
   SCIPlistMove(position, SCIPiteratorBegin(position.list));
}

/** Moves the element at 'position' to the end of the list. */
void SCIPlistMoveBack(ITERATOR position)
{
   SCIPlistMove(position, SCIPiteratorEnd(position.list));
}

/** Rearranges elements of list according to the ordering of order.
 *
 * example: list = (a b c d); order = (3 2 4 1)
 * after calling SCIPlistRearrange(list, order): list = (c b d a)
 * both lists must have the same size
 * order must have elements from 1 to list->size */
SCIP_RETCODE SCIPlistRearrange(SCIP* scip, LIST* list, LIST* order)
{
   LIST* new_list;
   ITERATOR it1, it2;
   int i;
   if( list && order && list->size == order->size)
   {
      new_list = SCIPlistCreate(scip);
      for(it1 = SCIPiteratorBegin(order); ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(order)); SCIPiteratorNext(&it1))
      {
         for( it2 = SCIPiteratorBegin(list), i = 1; i < *(int*)it1.node->data; ++i)
         {
            SCIPiteratorNext(&it2);
         }
         SCIPlistPushBack(scip, new_list, it2.node->data);
      }
      SCIPlistAssignList(list, new_list);
      SCIPlistDelete(scip, new_list);
      return SCIP_OKAY;
   }
   else
   {
      return SCIP_ERROR;
   }
}

/** Print function for strings.
 *
 * Can be used as a parameter in SCIPlistPrint(LIST* list, int(*printfunc)(void*)). */
int printstring(void *s)
{
   printf("%s\n", (char *)s);
   return 0;
}

/** Print function for integers.
 *
 * Can be used as a parameter in SCIPlistPrint(LIST* list, int(*printfunc)(void*)). */
int printint(void *i)
{
   printf("%i ", *(int*)i);
   return 0;
}

/** Returns an iterator pointing to the first element of list. */
ITERATOR SCIPiteratorBegin(LIST* list)
{
   ITERATOR it = {NULL, NULL};
   if(list)
   {
      it.list = list;
      it.node = list->nil->next;
   }
   return it;
}

/** Returns an iterator pointing to the past-the-end element of list. */
ITERATOR SCIPiteratorEnd(LIST* list)
{
   //the past-the-end element of list is the sentinal node 'nil'
   ITERATOR it = {NULL, NULL};
   if(list)
   {
      it.list = list;
      it.node = list->nil;
   }
   return it;
}

/** Sets the iterator 'it' to the next element. */
SCIP_Bool SCIPiteratorNext(ITERATOR* it)
{
   if(it && it->list && it->node)
   {
      it->node = it->node->next;
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** Sets the iterator 'it' to the previous element. */
SCIP_Bool SCIPiteratorPrev(ITERATOR* it)
{
   if(it && it->list && it->node)
   {
      it->node = it->node->prev;
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** Assigns it1 = it2. */
void SCIPiteratorEquals(ITERATOR* it1, ITERATOR* it2)
{
   it1->list = it2->list;
   it1->node = it2->node;
}

/** Returns TRUE if it1 points to the same element as it2. */
SCIP_Bool SCIPiteratorIsEqual(ITERATOR it1, ITERATOR it2)
{
   if(it1.list == it2.list && it1.node == it2.node)
   {
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/*
 * Data structures
 */

/** A struct that contains 4 hashmaps, which maps variables and constraints to their position in the constraint matrix (Ax<=b) and vice versa */
struct IndexMap
{
   /** index in problem -> constraint */
   SCIP_HASHMAP* indexcons;
   /** constraint -> index in problem */
   SCIP_HASHMAP* consindex;
   /** index in problem -> variable */
   SCIP_HASHMAP* indexvar;
   /** variable -> index in problem */
   SCIP_HASHMAP* varindex;
};
typedef struct IndexMap INDEXMAP;

/** detector data */
struct DEC_DetectorData
{
   DECDECOMP* decdecomp;
   SCIP_VAR*** varsperblock;
   int* nvarsperblock;
   SCIP_CONS*** consperblock;
   int *nconsperblock;
   SCIP_VAR** linkingvars;
   int nlinkingvars;
   SCIP_CONS** linkingconss;
   int nlinkingconss;
   int blocks;
   int maxblocks;
   int minblocks;
   SCIP_CONS** relevantConss; //array with all non-empty constraints
   int nRelevantConss;        //number of relevants constraints
   INDEXMAP* indexmap;
   int* ibegin; //array, ibegin[i]: index of first nonzero entry in row i
   int* iend;   //array, iend[i]: index of last nonzero entry in row i
   int* jbegin; //array, jbegin[j]: index of first nonzero entry in column j
   int* jend;   //array, jend[j]: index of last nonzero entry in column j
   int* jmin;   //array, jmin[i]: index of first nonzero column of the i-th row
   int* jmax;   //array, jmax[i]: the last nonzero entry among all rows prior to and including the i-th row
   int* minV;   //array, minV[i]: number of linking variables corresponding to a partitioning after the i-th row
   int* width;  //array, width[i]: width of the band (of nonzero entries after ROC) at row i
   int* hashmapindices;  //array with integers running from 0 to maximum(nvars, ncons)+1 (for usage of hash maps)
   LIST* rowsWithConstrictions;
   LIST* blockedAfterrow;
   SCIP_CLOCK* clock;
   SCIP_Bool found;
   int priority;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** Allocates memory for an indexmap. */
static
SCIP_RETCODE indexmapCreate(
      SCIP* scip,          /**< SCIP data structure  */
      INDEXMAP** indexmap, /**< address of the pointer which shall store the index map*/
      int nconss,          /**< number of constraints */
      int nvars            /**< number of variables */
      )
{
   INDEXMAP* imap;
   assert(scip != NULL);
   assert(nconss > 0);
   assert(nvars > 0);
   SCIPallocMemory(scip, &imap);
   assert(imap != NULL);

   SCIP_CALL(SCIPhashmapCreate(&imap->indexvar, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&imap->varindex, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&imap->indexcons, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&imap->consindex, SCIPblkmem(scip), nconss));

   *indexmap = imap;
   return SCIP_OKAY;
}

/** Deallocates memory of indexmap. */
static
void indexmapFree(SCIP* scip, INDEXMAP* indexmap)
{
   SCIPhashmapFree(&indexmap->indexvar);
   SCIPhashmapFree(&indexmap->varindex);
   SCIPhashmapFree(&indexmap->indexcons);
   SCIPhashmapFree(&indexmap->consindex);
   SCIPfreeMemory(scip, &indexmap);
}

/** Predicate function for sorting arrays.
 *
 * Returns the value of a - b (after casting into integers). */
static
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

/** Predicate function to compare two integers.
 *
 * Returns true if the value of a is equal to the value of b (a==b). */
static
SCIP_Bool compare_int(void* a, void * b)
{
   return ( *(int*)a == *(int*)b ? TRUE : FALSE );
}

/** Returns the maximum of a and b: max(a,b). */
static
int maximum(int a, int b)
{
   return (a > b ? a : b);
}

/** Returns the value of the maximum in the array a.
 *
 *  Or 0 if a is empty or invalid.*/
static
int maxArray(int* a, int num_elements)
{
   int i;
   int max;
   if(num_elements > 0 && a != NULL)
   {
      max = a[0];
      for (i = 1; i<num_elements; i++)
      {
         if (a[i] > max)
         {
            max = a[i];
         }
      }
      return(max);
   }
   //case: empty array
   else
   {
      return 0;
   }
}

/** Returns the value of the minimum in the array a.
 *
 *  Or 0 if a is empty or invalid. */
static
int minArray(int* a, int num_elements)
{
   int i;
   int min;
   if(num_elements > 0 && a != NULL)
   {
      min = a[0];
      for (i = 1; i<num_elements; i++)
      {
         if (a[i] < min)
         {
            min = a[i];
         }
      }
      return(min);
   }
   //case: empty array
   else
   {
      return 0;
   }
}

/** Returns the value of the minimum in the list between the iterators it1 and it2
 *
 *  Or -1 if a is empty or invalid. */
static
int minList(ITERATOR first, ITERATOR last)
{
   int min;
   //first is valid
   if(first.list && first.node && ! SCIPlistIsEmpty(first.list))
   {
      min = *(int*)first.node->data;
      SCIPiteratorNext(&first);
      for (;
           //first != last AND end of list is not reached yet
           ! (SCIPiteratorIsEqual(first, last)) && ! (SCIPiteratorIsEqual(first, SCIPiteratorEnd(first.list)));
           SCIPiteratorNext(&first)
           )
      {
         if ( *(int*)first.node->data < min )
         {
            min = *(int*)first.node->data;
         }
      }
      return min;
   }
   //first invalid
   else
   {
      return -1;
   }
}

/** Switches the data the pointers p1 and p2 points to. */
static
void switchPointers(void** p1, void** p2)
{
   void* p3; /* local for swap */
    p3 = *p2;
    *p2 = *p1;
    *p1= p3;
}

//debug ?
/** Creates a data and a gnuplot file for the initial problem.
 * @param scip < SCIP data structure
 * @param detectordata < presolver data data structure
 * @param filename name of the output files (without any filename extension) */
static
void plotInitialProblem(SCIP* scip, DEC_DETECTORDATA* detectordata, char* filename)
{
   FILE* output;
   char datafile[256];
   char gpfile[256];
   char pdffile[256];
   int i;
   int j;
   int* varindex;
   int* consindex;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   //filenames
   sprintf(datafile, "%s.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotProblem!\n");
   }
   else
   {
      for(i = 0; i < detectordata->nRelevantConss; ++i)
      {
         cons = detectordata->relevantConss[i];
         consindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->consindex, (void*) cons);
         assert(consindex != NULL);
         for(j = 0; j < SCIPgetNVarsXXX(scip, cons); ++j)
         {
            var = SCIPgetVarsXXX(scip, cons)[j];
            varindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->varindex, (void*) var);
            assert(varindex != NULL);
            fprintf(output, "%i %i\n", *varindex, *consindex);
         }
      }
   }
   fclose(output);

   //write Gnuplot file
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nunset xtics\nunset ytics\nunset border\nset pointsize 0.05\nset xrange [0:%i]\nset yrange[%i:0]\nplot '%s' lt 0 pt 5 notitle", pdffile, SCIPgetNVars(scip), detectordata->nRelevantConss, datafile);
   fclose(output);
}

//debug ?
static
/** Creates a data and a gnuplot file for the blocked problem.
 * @param scip < SCIP data structure
 * @param detectordata < presolver data data structure
 * @param filename name of the output files (without any filename extension) */
void plotBlocking(SCIP* scip, DEC_DETECTORDATA* detectordata, char* filename)
{
   FILE* output;
   char datafile[256];
   char gpfile[256];
   char pdffile[256];
   int i;
   int j;
   int k;
   int nvars;
   int* varindex;
   int* consindex;
   SCIP_VAR** vars;
   SCIP_CONS* cons;
   //filenames
   sprintf(datafile, "%s.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotBlocking!\n");
   }
   else
   {
      //loop over all blocks
      for(i = 0; i < detectordata->blocks; ++i)
      {
         //loop over all constraints in block
         for(j = 0; j < detectordata->nconsperblock[i]; ++j)
         {
            cons = detectordata->consperblock[i][j];
            consindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->consindex, (void*) cons);
            assert(consindex != NULL);
            nvars = SCIPgetNVarsXXX(scip, cons);
            vars = SCIPgetVarsXXX(scip, cons);
            //loop over all vars in constraint
            for(k = 0; k < nvars; ++k)
            {
               varindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->varindex, (void*) vars[k]);
               assert(varindex != NULL);
               fprintf(output, "%i %i\n", *varindex, *consindex);
            }//loop over all vars in constraint
         }//loop over all constraints in block
         fprintf(output, "\n");
      }//loop over all blocks
   }
   fclose(output);

   //write Gnuplot file
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nunset xtics\nunset ytics\nunset border\nset style line 1 lt 0 lw 1 pt 5\nset style line 2 lt 9 lw 1 pt 5\nset pointsize 0.05\nset xrange [0:%i]\nset yrange[%i:0]\nplot for [i=0:%i:1] '%s' every :::i::(i+1) linestyle (i%%2+1) notitle", pdffile, SCIPgetNVars(scip), detectordata->nRelevantConss, detectordata->blocks-1, datafile);
   fclose(output);
}

//debug ?
/** Creates a data and a gnuplot file for the graph representing the array minV (number of linking variables).
 * @param detectordata < presolver data data structure
 * @param filename name of the output files (without any filename extension) */
static
void plotMinV(SCIP* scip, DEC_DETECTORDATA* detectordata, char* filename)
{
   FILE* output;
   char datafile[256];
   char blockingfile[256];
   char gpfile[256];
   char pdffile[256];
   int i;
   ITERATOR it1;
   //filenames
   sprintf(datafile, "%s.dat", filename);
   sprintf(blockingfile, "%s_blocked_at.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);

   //datafile
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotMinV!\n");
   }
   else
   {
      //write data to datafile
      for(i = 0; i < detectordata->nRelevantConss -1; ++i)
      {
         fprintf(output, "%i\n", detectordata->minV[i]);
      }
   }
   fclose(output);

   //blocking points
   output = fopen(blockingfile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for blocking output in plotMinV!\n");
   }
   else
   {
      //write data to blockingfile
      for(it1 = SCIPiteratorBegin(detectordata->blockedAfterrow); ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(detectordata->blockedAfterrow)); SCIPiteratorNext(&it1))
      {
         fprintf(output, "%i %i\n", *(int*)it1.node->data - 1, detectordata->minV[*(int*)it1.node->data-1]);
      }
   }
   fclose(output);
   //write Gnuplot file
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nset style line 1 lt 1 lc rgb \"black\"\nplot '%s' title '# verb. Variablen' ls 1 with lines, \\\n '%s' lt 0 pt 4 with points title \"Blockgrenze\"", pdffile, datafile, blockingfile);
   fclose(output);
}


static
void writeParams(SCIP* scip, DEC_DETECTORDATA* detectordata, char* paramfile, int ROC_iterations, int tau, double time)
{
   FILE* output;
   int i;
   int nvars;
   int ncons;
   int nonzeros;
   int zeros;
   float sparsity;
   int minimum_linking_vars;
   output=fopen(paramfile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotMinV!\n");
   }
   else
   {
      nvars = SCIPgetNOrigVars(scip);
//      ncons = SCIPgetNConss(scip);
      ncons = detectordata->nRelevantConss;
      nonzeros = 0;
      for(i = 0; i < ncons; ++i)
      {
         nonzeros += SCIPgetNVarsXXX(scip,  SCIPgetConss(scip)[i]);
      }
      zeros = nvars*ncons - nonzeros;
      sparsity = (float) nonzeros / (nvars*ncons);
      SCIPdebugMessage("minList.\n");
      minimum_linking_vars = minList(SCIPiteratorBegin(detectordata->rowsWithConstrictions), SCIPiteratorEnd(detectordata->rowsWithConstrictions));
      SCIPdebugMessage("minList.\n");
      fprintf(output, "# of rows\n%i\n", ncons);
      fprintf(output, "# of columns\n%i\n", nvars);
      fprintf(output, "# of nonzeros\n%i\n", nonzeros);
      fprintf(output, "# of zeros\n%i\n", zeros);
      fprintf(output, "# sparsity\n%f\n", sparsity);
      fprintf(output, "# detection time in seconds\n%f\n", time);
      fprintf(output, "# tau\n%i\n", tau);
      fprintf(output, "# of blocks\n%i\n", detectordata->blocks);
      fprintf(output, "# of iterations\n%i\n", ROC_iterations);
      fprintf(output, "# of minimum linking vars\n%i\n", minimum_linking_vars);
      fprintf(output, "# of linking vars\n%i\n", detectordata->nlinkingvars);
      for(i = 0; i < detectordata->blocks; ++i)
      {
         fprintf(output, "block # %i\n", i+1);
         fprintf(output, "# nonlinking vars\n%i\n", detectordata->nvarsperblock[i]);
         fprintf(output, "# cons per block\n%i\n", detectordata->nconsperblock[i]);
      }
      fclose(output);
   }
}

/** Scans all constraints of the constraint array of the scip object,
 * and stores pointers to all constraints that have at least one variable in detectordata->relevantConss.
 * Thus it removes all empty constraints.
 */
static
SCIP_RETCODE findRelevantConss(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   SCIP_CONS** cons_array;
   LIST* relevantConssIndices;
   int i;
   int* data;
   ITERATOR it1;
   cons_array = SCIPgetConss(scip);
   relevantConssIndices = SCIPlistCreate(scip);
   for(i = 0; i < SCIPgetNConss(scip); ++i)
   {
      if(SCIPgetNVarsXXX(scip, cons_array[i]) > 0)
      {
         SCIPallocMemory(scip, &data);
         *data = i;
         SCIPlistPushBack(scip, relevantConssIndices, data);
      }
   }
   //debug
//   SCIPlistPrint(relevantConssIndices, printint);
   //allocate memory for detectordata->relevantConss and store pointers of relevant conss
   detectordata->nRelevantConss = relevantConssIndices->size;
   SCIPdebugMessage("nRelevantConss: %i \n", detectordata->nRelevantConss);
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->relevantConss, detectordata->nRelevantConss));
   for(i = 0, it1 = SCIPiteratorBegin(relevantConssIndices); i < detectordata->nRelevantConss; ++i, SCIPiteratorNext(&it1))
   {
      detectordata->relevantConss[i] = cons_array[*(int*)it1.node->data];
   }
   SCIPlistDeleteData(scip, relevantConssIndices);
   SCIPlistDelete(scip, relevantConssIndices);

   return SCIP_OKAY;
}

/** Creates a nested list with the indices of the nonzero entries of each row.
 *
 * example:
 * constraint matrix:
 *
 *  1 1 0 1 0
 *
 *  0 1 1 0 0
 *
 *  0 0 0 0 1
 *
 *  resulting list:
 *  ( (1 2 4)
 *    (2 3)
 *    (5)    )
 */
static
LIST* rowindices_list(
      SCIP* scip,                      /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      SCIP_HASHMAP* indexcons,         /**< hashmap index -> constraint */
      SCIP_HASHMAP* varindex           /**< hashmap variable -> index*/
      )
{
   //create the rowindices list
   int i;
   int j;
   int* data;
   LIST* rowindices;
   LIST* rowindices_row;
   int ncons; //number of constraints of the problem
   int nconsvars; //number of variables in a constraint
   int* probindices;
   int* hashmapindex;
   SCIP_CONS* cons; //one constraint of the problem
   SCIP_VAR** consvars; //array of variables that occur in a constraint (unequal zero)

   rowindices = SCIPlistCreate(scip);
   ncons = detectordata->nRelevantConss;
   for(i = 0; i < ncons; ++i)
   {
      hashmapindex = &detectordata->hashmapindices[i+1];
      cons = (SCIP_CONS*) SCIPhashmapGetImage(indexcons, (void*) hashmapindex);
      nconsvars = SCIPgetNVarsXXX(scip, cons);
      consvars = SCIPgetVarsXXX(scip, cons);
      //allocate memory for the array of probindices
      SCIPallocMemoryArray(scip, &probindices, nconsvars);
      //fill the array with the indices of the variables of the current constraint
      for(j = 0; j < nconsvars; ++j)
      {
//         probindices[j] = SCIPvarGetProbindex(consvars[j])+1;
         probindices[j] = *(int*) SCIPhashmapGetImage(varindex, consvars[j]);
      }
      //sort the elements of probindices ('<')
      qsort(probindices, nconsvars, sizeof(int), compare);
      //store a copy of the elements of probindices in the list rowindices_row
      rowindices_row = SCIPlistCreate(scip);
      for(j = 0; j < nconsvars; ++j)
      {
         SCIPallocMemory(scip, &data);
         *data = probindices[j];
         SCIPlistPushBack(scip, rowindices_row, data);
      }
      //deallocate memory
      SCIPfreeMemoryArray(scip, &probindices);
      //add rowindices_row to the list rowindices
      SCIPlistPushBack(scip, rowindices, rowindices_row);
   }
   return rowindices;
}

/** Creates a nested list with the indices of the nonzero entries of each column.
 *
 * example:
 *
 * constraint matrix:
 *
 *  1 1 0 1 0
 *
 *  0 1 1 0 0
 *
 *  0 0 0 0 1
 *
 *  resulting list:
 *  ( (1)
 *    (1 2)
 *    (2)
 *    (1)
 *    (3)    )
 */
static
LIST* columnindices_list(
      SCIP* scip,                      /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata,  /**< detector data data structure */
      LIST* rowindices                 /**< A list with the row indices (achieved from calling rowindices_list() ) */
      )
{
   LIST** columnindices_array;
   LIST* columnindices;
   LIST* rowindices_row;
   int* data;
   int position;
   int nvars;
   int i;
   ITERATOR it1;
   ITERATOR it2;
   nvars = SCIPgetNVars(scip);
   //create the columnindices_array with pointers to empty lists
   SCIPallocMemoryArray(scip, &columnindices_array, nvars);
   for(i = 0; i < nvars; ++i)
   {
      columnindices_array[i] = SCIPlistCreate(scip);
   }

   for(it1 = SCIPiteratorBegin(rowindices), i = 0; ! ( SCIPiteratorIsEqual(it1, SCIPiteratorEnd(rowindices)) ); SCIPiteratorNext(&it1), ++i)
   {
      rowindices_row = it1.node->data;
      for(it2 = SCIPiteratorBegin(rowindices_row); ! ( SCIPiteratorIsEqual(it2, SCIPiteratorEnd(rowindices_row)) ); SCIPiteratorNext(&it2))
      {
         SCIPallocMemory(scip, &data);
         *data = i+1;
         position = *(int*)(it2.node->data)-1;
         SCIPlistPushBack(scip, columnindices_array[position], data);
      }
   }
   //create a columnindices list instead of an array
   columnindices = SCIPlistCreate(scip);
   for(i = 0; i < nvars; ++i)
   {
      SCIPlistPushBack(scip, columnindices, columnindices_array[i]);
   }
   //deallocate memory
   SCIPfreeMemoryArray(scip, &columnindices_array);
   return columnindices;
}

/** Does the row ordering of the ROC2 algorithm.
 *
 * It also works for the column ordering. In this case the terms row<->column have to be exchanged.
 *
 * @param columnindices A list of the nonzero entries in each column.
 * @param nrows The number of rows of the constraint matrix (=number of relevant constraints)
 * @return A list with the new row order. E.g. (2 3 1) means the second row comes first now, and so on. */
static
LIST* rowOrdering(SCIP* scip, LIST* columnindices, int nrows)
{
   LIST* roworder;
   LIST* new_roworder;
   ITERATOR it1;
   ITERATOR it2;
   ITERATOR it3;
   ITERATOR it4;

   //create a list for the order of the rows ( 1 2 3 ... nrows )
   roworder = SCIPlistCreateInt(scip, 1, nrows);
   new_roworder = SCIPlistCopyShallow(scip, roworder);

   for( it1 = SCIPiteratorEnd(columnindices), SCIPiteratorPrev(&it1); it1.node != it1.list->nil; SCIPiteratorPrev(&it1))
   {
      for( it2 = SCIPiteratorEnd(roworder), SCIPiteratorPrev(&it2); it2.node != it2.list->nil; SCIPiteratorPrev(&it2) )
      {
         it3 = SCIPlistFind(SCIPiteratorBegin(it1.node->data), SCIPiteratorEnd(it1.node->data), it2.node->data, compare_int);
         if(it3.node->data != NULL)
         {
            it4 = SCIPlistFind(SCIPiteratorBegin(new_roworder), SCIPiteratorEnd(new_roworder), it2.node->data, compare_int);
            SCIPlistMoveFront(it4);
         }
         else
         {
            ;
         }
      }
      SCIPlistAssignList(roworder, new_roworder);
   }
   //deallocate memory
   SCIPlistDelete(scip, new_roworder);
   return roworder;
}

/** Stores the first and last entry of the i-th column(row) in begin[i] and end[i] respectively.
 *
 * @param begin Array to store the first nonzero entry of the i-th column (row)
 * @param end Array to store the last nonzero entry of the i-th column (row)
 * @param indices columnindices list (rowindices list) */
static
SCIP_RETCODE formIndexArray(int* begin, int* end, LIST* indices)
{
   ITERATOR it1;
   ITERATOR it2;
   int i;
   assert(begin != NULL && end != NULL && indices != NULL);
   for(it1 = SCIPiteratorBegin(indices), i = 0; ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(indices)); SCIPiteratorNext(&it1), ++i)
   {
      //case: list not empty
      if(! SCIPlistIsEmpty(it1.node->data))
      {
         it2 = SCIPiteratorBegin(it1.node->data);
         begin[i] = *(int*)it2.node->data;
         it2 = SCIPiteratorEnd(it1.node->data);
         SCIPiteratorPrev(&it2);
         end[i] = *(int*)it2.node->data;
      }
      //case: list empty
      else
      {
         begin[i] = 0;
         end[i] = 0;
      }
   }
   return SCIP_OKAY;
}

//static
//SCIP_Bool IndexArrayChanges(SCIP* scip, DEC_DETECTORDATA* detectordata)
//{
//   /*returns TRUE if at least one entry of a new index array differs from the corresponding entry of the old index array*/
//   int i;
//   //any change in ibegin or iend?
//   for(i = 0; i < detectordata->nRelevantConss; ++i)
//   {
//      if(   detectordata->ibegin_old[i] != detectordata->ibegin_new[i]
//         || detectordata->iend_old[i] != detectordata->iend_new[i] )
//      {
//         return TRUE;
//      }
//   }
//   //any change in jbegin or jend?
//   for(i = 0; i < SCIPgetNVars(scip); ++i)
//   {
//      if(   detectordata->jbegin_old[i] != detectordata->jbegin_new[i]
//         || detectordata->jend_old[i] != detectordata->jend_new[i] )
//      {
//         return TRUE;
//      }
//   }
//   //case: all entries of old and new are equal
//   return FALSE;
//}

/**Returns FALSE if at least one entry of new_array and old_array are different.*/
static
SCIP_Bool arraysAreEqual(int* new_array, int* old_array, int num_elements)
{
   int i;
   for(i = 0; i < num_elements; ++i)
   {
      if( new_array[i] != old_array[i] )
      {
         return FALSE;
      }
   }
   //case: all entries of old and new are equal
   return TRUE;
}

/**Permutes the order of rows and columns in inputmap and stores the result in outputmap.
 *
 *  One call of this function is equivalent to one iteration of the ROC2-algortihm. */
static
SCIP_RETCODE rankOrderClustering(
      SCIP*             scip,          /**< SCIP data structure */
      DEC_DETECTORDATA*  detectordata, /**< presolver data data structure */
      INDEXMAP* inputmap,              /**< indexmap for input */
      INDEXMAP* outputmap              /**< indexmap for output */
      )
{
   LIST* roworder;
   LIST* columnorder;
   LIST* rowindices;
   LIST* columnindices;
   ITERATOR it1;
   int nvars;
   int ncons;
   int i;
   int position;
   int* hashmapindex;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(detectordata != NULL);
   nvars = SCIPgetNVars(scip);
   ncons = detectordata->nRelevantConss;
   //create the lists containing the positions of nonzero entries; row and column ordering
   rowindices = rowindices_list(scip, detectordata, inputmap->indexcons, inputmap->varindex);
   columnindices = columnindices_list(scip, detectordata, rowindices);
   roworder = rowOrdering(scip, columnindices, ncons);
   SCIPlistRearrange(scip, rowindices, roworder);
   columnorder = rowOrdering(scip, rowindices, nvars);

   //consindex and indexcons
   for(it1 = SCIPiteratorBegin(roworder), i = 0; ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(roworder)) && i < ncons; ++i, SCIPiteratorNext(&it1))
   {
      position = (*(int*)it1.node->data);
      hashmapindex = &detectordata->hashmapindices[position];
      cons = SCIPhashmapGetImage(inputmap->indexcons, (void*) hashmapindex);
      assert ( cons != NULL);
      //consindex
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( SCIPhashmapExists(outputmap->consindex, (void*) cons));
      SCIPhashmapSetImage(outputmap->consindex, (void*) cons, (void*) hashmapindex);
      //indexcons
      assert( SCIPhashmapExists(outputmap->indexcons, (void*) hashmapindex ));
      SCIPhashmapSetImage(outputmap->indexcons, (void*) hashmapindex, cons);
   }
   //varindex and indexvar
   for(it1 = SCIPiteratorBegin(columnorder), i = 0; ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(columnorder)) &&i < nvars; ++i, SCIPiteratorNext(&it1))
   {
      position = (*(int*)it1.node->data);
      hashmapindex = &detectordata->hashmapindices[position];
      var = (SCIP_VAR*) SCIPhashmapGetImage(inputmap->indexvar, (void*) hashmapindex);
      assert ( var != NULL);

      //varindex
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( SCIPhashmapExists(outputmap->varindex, (void*) var) );
      SCIPhashmapSetImage(outputmap->varindex, (void*) var, (void*) hashmapindex);
      //indexvar
      assert( SCIPhashmapExists(outputmap->indexvar, (void*) hashmapindex ));
      SCIPhashmapSetImage(outputmap->indexvar, (void*) hashmapindex, var);
   }
   //deallocate memory
   SCIPlistDeleteData(scip, roworder);
   SCIPlistDelete(scip, roworder);
   SCIPlistDeleteData(scip, columnorder);
   SCIPlistDelete(scip, columnorder);
   SCIPlistDeleteNested(scip, rowindices);
   SCIPlistDeleteNested(scip, columnindices);
   return SCIP_OKAY;
}

/** finds rows with local minima regarding the number of linking variables and stores them in detectordata->rowsWithConstrictions */
static
void rowsWithConstriction(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   //if blocking is performed after row i+1; local minima
   int i;
   int* data;
   for(i = 1; i < detectordata->nRelevantConss - 2; ++i)
   {
      //is minV[i] a local minimum?    < or <=   ? What does make more sense?
      if(detectordata->minV[i] < detectordata->minV[i-1] && detectordata->minV[i] < detectordata->minV[i+1])
      {
         SCIPallocMemory(scip, &data);
         *data = i+1;
         SCIPlistPushBack(scip, detectordata->rowsWithConstrictions, data);
      }
   }
}

/** assigns variables to a block, divided into linking variables and nonlinking variables */
static
void assignVarsToBlock(
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      int block,                       /**< number of current block */
      int first_var,                   /**< index of first variable in the block*/
      int last_var,                    /**< index of last variable in the block*/
      int first_linkingvar             /**< index of first linking variable*/
      )
{
   int i;
   int j;
   int* hashmapindex;
   SCIP_VAR* var;
   //assign the subscipvars (=nonlinking vars)
   assert(first_linkingvar >= first_var);
   detectordata->nvarsperblock[block-1] = first_linkingvar - first_var;
   for(i = first_var, j = 0; i < first_linkingvar; ++i, ++j)
   {
      hashmapindex = &detectordata->hashmapindices[i];
      var = (SCIP_VAR*) SCIPhashmapGetImage(detectordata->indexmap->indexvar, (void*) hashmapindex);
      assert(var != NULL);
      detectordata->varsperblock[block-1][j] = var;
   }
   //assign linking vars
   for(i = first_linkingvar; i <= last_var; ++i)
   {
      hashmapindex = &detectordata->hashmapindices[i];
      var = (SCIP_VAR*) SCIPhashmapGetImage(detectordata->indexmap->indexvar, (void*) hashmapindex);
      assert(var != NULL);
      detectordata->linkingvars[detectordata->nlinkingvars] = var;
      ++detectordata->nlinkingvars;
   }
}

/** assigns constraints in the interval [first_cons, last_cons] to 'block' */
static
void assignConsToBlock(SCIP* scip, DEC_DETECTORDATA* detectordata, int block, int first_cons, int last_cons)
{
   int i;
   int j;
   int* hashmapindex;
   SCIP_CONS* cons;
   //assign the constraints to the current block
   detectordata->nconsperblock[block-1] = last_cons - first_cons + 1;
   for(i = first_cons, j = 0; i <= last_cons; ++i, ++j)
   {
      hashmapindex = &detectordata->hashmapindices[i];
      cons = (SCIP_CONS*) SCIPhashmapGetImage(detectordata->indexmap->indexcons, (void*) hashmapindex);
      assert(cons != NULL);
      detectordata->consperblock[block-1][j] = cons;
   }
   SCIPlistPushBack(scip, detectordata->blockedAfterrow, &detectordata->hashmapindices[last_cons]);
}

/** returns the largest column index of a nonzero entry between rows [from_row, to_row] */
static
int getMaxColIndex(DEC_DETECTORDATA* detectordata, int from_row, int to_row)
{
   //some pointer arithmetic
   return maxArray(detectordata->iend + (from_row -1), to_row - from_row + 1);
}

/** returns the column index of the first nonzero entry in 'row' */
static
int getMinColIndex(DEC_DETECTORDATA* detectordata, int row)
{
   return detectordata->ibegin[row];
}

/** determines if a blocking at 'block_at_row' is a valid blocking
 *
 * @param detectordata detectordata data structure
 * @param prev_block_first_row first row of the previous block
 * @param prev_block_last_row last row of the previous block
 * @param block_at_row the row for which you want to determine if the blocking is valid
 * @return TRUE if blocking is valid, else FALSE
 */
static
SCIP_Bool isValidBlocking(DEC_DETECTORDATA* detectordata, int prev_block_first_row, int prev_block_last_row, int block_at_row)
{
   int last_column_prev_block;
   int first_column_current_block;

   last_column_prev_block = getMaxColIndex(detectordata, prev_block_first_row, prev_block_last_row);
   first_column_current_block = getMinColIndex(detectordata, block_at_row);
   return ( first_column_current_block > last_column_prev_block ? TRUE : FALSE);
}

/** this functions looks for rows to block at, which creates block of size min_block_size or bigger
 *
 * @param it_constrictions Iterator pointing to a list of constraints (detectordata->rowsWithConstrictions)
 * @param min_block_size minimum number of rows to be in a block
 * @param prev_block_last_row the last row of the preceding block
 * @return Iterator pointing to a node which contains a suitable row for blocking; If the iterator points after the last element, no candidate was found
 */
static
ITERATOR findBlockingCandidate(ITERATOR it_constrictions, int min_block_size, int prev_block_last_row)
{
   for(;;)
   {
      //end of the list?
      if( SCIPiteratorIsEqual(it_constrictions, SCIPiteratorEnd(it_constrictions.list)))
      {
         return it_constrictions;
      }
      //does a blocking to the next row forming a constriction comprise more rows than min_block_size
      if( (* (int*) (it_constrictions.node->data) - prev_block_last_row) >= min_block_size )
      {
         return it_constrictions;
      }
      //advance iterator to next element
      SCIPiteratorNext(&it_constrictions);
   }
}

/** this functions determines the next row to block at
 *
 * @param detectordata detectordata data structure
 * @param it_constrictions Iterator pointing to a list of constraints (detectordata->rowsWithConstrictions)
 * @param min_block_size minimum number of rows to be in a block
 * @param prev_block_first_row the first row of the preceding block
 * @param prev_block_last_row the last row of the preceding block
 * @return Iterator pointing to a node which contains a suitable row for blocking; If the iterator points after the last element, no row was found
 */
static
ITERATOR nextRowToBlockAt(DEC_DETECTORDATA* detectordata, ITERATOR it_constrictions, int min_block_size, int prev_block_first_row, int prev_block_last_row)
{
   assert(it_constrictions.list != NULL);
   assert(it_constrictions.node != NULL);

   //end of the constriction list?
   if( SCIPiteratorIsEqual(it_constrictions, SCIPiteratorEnd(it_constrictions.list)))
   {
      return it_constrictions;
   }

   for(;;)
   {
      //find a blocking candidate
      it_constrictions = findBlockingCandidate(it_constrictions, min_block_size, prev_block_last_row);
      //case: no candidate found
      if( SCIPiteratorIsEqual(it_constrictions, SCIPiteratorEnd(it_constrictions.list)))
      {
         break;
      }
      //case: candidate found
      else
      {
         //valid blocking
         if(isValidBlocking(detectordata, prev_block_first_row, prev_block_last_row, *(int*) it_constrictions.node->data))
         {
            break;
         }
         //invalid blocking
         else
         {
            SCIPiteratorNext(&it_constrictions);
         }
      }
   }
   return it_constrictions;
}

/** Tries to divide the problem into subproblems (blocks)*/
static
void blocking(
      SCIP* scip,                      /**< scip object */
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      int tau,                         /**< desired number of blocks */
      int nvars                        /**< number of variables in the problem*/
      )
{
   int block;
   int prev_block_first_row;
   int prev_block_last_row;
   int current_row;
   int min_block_size;
   //notation: i=current block; im1=i-1=previous block; ip1=i+1=next block
   int max_col_index_im1;
   int min_col_index_ip1;
   int max_col_index_i;
   ITERATOR it1;
   //debug
   SCIPdebugMessage("Starting Blocking...\n");
   SCIPdebugMessage("Max blocks: %i\n", detectordata->maxblocks);
   block = 1;
   prev_block_first_row = 0;
   prev_block_last_row = 0;
   max_col_index_im1 = 0;
   min_col_index_ip1 = 1;
   min_block_size = round( (float) detectordata->nRelevantConss / (2 * tau ));
   it1 = SCIPiteratorBegin(detectordata->rowsWithConstrictions);

   for( it1 = nextRowToBlockAt(detectordata, it1, min_block_size, prev_block_first_row, prev_block_last_row);
         ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(it1.list)) && block < detectordata->maxblocks;
         it1 = nextRowToBlockAt(detectordata, it1, min_block_size, prev_block_first_row, prev_block_last_row) )
   {
      current_row = * (int*) (it1.node->data);
      max_col_index_i = getMaxColIndex(detectordata, prev_block_last_row + 1, current_row);
      min_col_index_ip1 = getMinColIndex(detectordata, current_row + 1);
      SCIPdebugMessage("assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, current_row);
      SCIPdebugMessage("vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, max_col_index_i, min_col_index_ip1, max_col_index_i);
      //assign the variables and constraints to block
      assignVarsToBlock(detectordata, block, max_col_index_im1+1, max_col_index_i, min_col_index_ip1);
      assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, current_row);
      //update variables in the while loop
      max_col_index_im1 = max_col_index_i;
      prev_block_first_row = prev_block_last_row + 1;
      prev_block_last_row = current_row;
      ++block;
   }
   //assign the remaining (< M/2tau) cons and vars to the last block; no new linking vars are added
   //debug
   SCIPdebugMessage("last time: assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, detectordata->nRelevantConss);
   SCIPdebugMessage("last time: vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, nvars, nvars+1, nvars);
   assignVarsToBlock(detectordata, block, max_col_index_im1+1, nvars, nvars+1);
   assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, detectordata->nRelevantConss);
   SCIPlistPopBack(scip, detectordata->blockedAfterrow);
   detectordata->blocks = block;
}


/** copies the variable and block information to the decomp structure */
static
SCIP_RETCODE copyDetectorDataToDecomp(
      SCIP*             scip,         /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata, /**< presolver data data structure */
      DECDECOMP*        decomp        /**< DECOMP data structure */
      )
{
   int i;
   assert(scip != 0);
   assert(detectordata != 0);
   assert(decomp != 0);

   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipvars, detectordata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipconss, detectordata->blocks));

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->linkingvars, detectordata->linkingvars, detectordata->nlinkingvars));
   decomp->nlinkingvars = detectordata->nlinkingvars;
   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->linkingconss, detectordata->linkingconss, detectordata->nlinkingconss));
   decomp->nlinkingconss = detectordata->nlinkingconss;

   for( i = 0; i < detectordata->blocks; ++i )
   {
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipconss[i], detectordata->consperblock[i], detectordata->nconsperblock[i]));
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipvars[i], detectordata->varsperblock[i], detectordata->nvarsperblock[i]));
   }

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipconss, detectordata->nconsperblock, detectordata->blocks));
   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipvars, detectordata->nvarsperblock, detectordata->blocks));

   decomp->varindex = detectordata->indexmap->varindex;
   decomp->consindex = detectordata->indexmap->consindex;
   decomp->nblocks = detectordata->blocks;
   decomp->type = DEC_DECTYPE_STAIRCASE;
   return SCIP_OKAY;
}

static
DEC_DECL_INITDETECTOR(initStairheur)
{

   int i;
   int nvars;
   int nconss;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   SCIP_CALL(SCIPcreateWallClock(scip, &detectordata->clock));
   SCIP_CALL(SCIPstartClock(scip, detectordata->clock));
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nconsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nvarsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->linkingvars, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->linkingconss, nconss));
   for( i = 0; i < detectordata->maxblocks; ++i )
   {
      detectordata->nvarsperblock[i] = 0;
      detectordata->nconsperblock[i] = 0;
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock[i], nconss));
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock[i], nvars));
   }

   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->ibegin, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->iend, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jbegin, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jend, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jmin, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jmax, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->minV, nconss-1));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->width, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->hashmapindices, maximum(nvars, nconss) + 1));
   for(i = 0; i < maximum(nvars, nconss)+1; ++i)
   {
      detectordata->hashmapindices[i] = i;
   }
   detectordata->rowsWithConstrictions = SCIPlistCreate(scip);
   detectordata->blockedAfterrow = SCIPlistCreate(scip);

   detectordata->nlinkingvars = 0;
   detectordata->nlinkingconss = 0;
   /* create hash tables */
   indexmapCreate(scip, &detectordata->indexmap, nconss, nvars);
//   SCIP_CALL(SCIPhashmapCreate(&detectordata->indexvar, SCIPblkmem(scip), nvars));
//   SCIP_CALL(SCIPhashmapCreate(&detectordata->varindex, SCIPblkmem(scip), nvars));
//   SCIP_CALL(SCIPhashmapCreate(&detectordata->indexcons, SCIPblkmem(scip), nconss));
//   SCIP_CALL(SCIPhashmapCreate(&detectordata->consindex, SCIPblkmem(scip), nconss));
   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DEC_DECL_EXITDETECTOR(exitStairheur)
{
   int i;
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   if( !detectordata->found)
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   /* free presolver data */
   for( i = 0; i < detectordata->maxblocks; ++i)
   {
      SCIPfreeMemoryArray(scip, &detectordata->consperblock[i]);
      SCIPfreeMemoryArray(scip, &detectordata->varsperblock[i]);
   }
   SCIPfreeMemoryArray(scip, &detectordata->varsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nvarsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->consperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nconsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->linkingvars);
   SCIPfreeMemoryArray(scip, &detectordata->linkingconss);
   SCIPfreeMemoryArray(scip, &detectordata->relevantConss);

   SCIPfreeMemoryArray(scip, &detectordata->ibegin);
   SCIPfreeMemoryArray(scip, &detectordata->iend);
   SCIPfreeMemoryArray(scip, &detectordata->jbegin);
   SCIPfreeMemoryArray(scip, &detectordata->jend);
   SCIPfreeMemoryArray(scip, &detectordata->jmin);
   SCIPfreeMemoryArray(scip, &detectordata->jmax);
   SCIPfreeMemoryArray(scip, &detectordata->minV);
   SCIPfreeMemoryArray(scip, &detectordata->width);
   SCIPfreeMemoryArray(scip, &detectordata->hashmapindices);
   SCIP_CALL(SCIPfreeClock(scip, &detectordata->clock));
   //delete lists
   //data had to be deleted before because of SCIP memory management
   SCIPlistDelete(scip, detectordata->rowsWithConstrictions);
   SCIPlistDelete(scip, detectordata->blockedAfterrow);
   //free deep copied hash maps
   //DO NOT FREE varindex and consindex because they are only shallow copied and contain the final permutation
   //debug
   SCIPhashmapFree(&detectordata->indexmap->indexvar);
   SCIPhashmapFree(&detectordata->indexmap->indexcons);
   SCIPfreeMemory(scip, &detectordata->indexmap);
   SCIPfreeMemory(scip, &detectordata);
   return SCIP_OKAY;
}

static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildStair)
{
   int i;
   int* hashmapindex;
   int ncons; //number of constraints in the problem
   int nvars; //number of variables in the problem
   SCIP_VAR** vars_array;
   SCIP_VAR* var;
   SCIP_CONS** cons_array;
   SCIP_CONS* cons;
//   DEC_DETECTOR* stairheur;
//   DEC_DETECTORDATA* detectordata;
   LIST* rowindices;
   LIST* columnindices;
   INDEXMAP* indexmap_permuted;
   int* ibegin_permuted;
   int* iend_permuted;
   int* jbegin_permuted;
   int* jend_permuted;
   int n;
   int v;
   int tau;
   int ROC_iterations;

   assert(scip != NULL);
//   stairheur = DECfindDetector(scip, DEC_DETECTORNAME);
//   detectordata = DECdetectorGetData(stairheur);
   assert(detectordata != NULL);
//   assert(strcmp(DECdetectorGetName(stairheur), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("%s\n", SCIPgetProbName(scip));
   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   //remove empty constraints
   SCIP_CALL( findRelevantConss(scip, detectordata) );
   SCIP_CALL(DECdecdecompCreate(scip, &detectordata->decdecomp));
   nvars = SCIPgetNVars(scip);
   vars_array = SCIPgetVars(scip);
   ncons = detectordata->nRelevantConss;
   cons_array = detectordata->relevantConss;
#ifdef SCIP_DEBUG
   {
      //remove '/' from problem name
//      const char* pname;
//      const char* ext;
//      char fname[256];
//      pname = strrchr(SCIPgetProbName(scip), '/');
//      if( pname == NULL )
//      {
//         pname = SCIPgetProbName(scip);
//      }
//      else
//      {
//         pname = pname+1;
//      }
//      ext = strrchr(SCIPgetProbName(scip), '.');
//      if( ext == NULL )
//      {
//         strcpy(fname, pname);
//      }
//      else
//      {
//         strncpy(fname, pname, ext-pname);
//         fname[ext-pname]='\0';
//      }
//      strcat(fname, "_permuted.lp");
//      sprintf(fname, "permuted_%s", pname);
//
//      SCIPdebugMessage("Writing Orig Problem to %s\n", fname);
//      SCIPwriteOrigProblem  (scip, fname, "lp", TRUE);
//      SCIPdebugMessage("Writing Trans Problem\n");
//      SCIPwriteTransProblem  (scip, fname, "lp", TRUE);


      SCIPdebugMessage("ncons: %i \n", ncons);
      SCIPdebugMessage("nvars: %i \n", nvars);
      SCIPdebugMessage("initializing hash maps\n");
   }
#endif
   //initialize hash maps for variables: indexvar_old, indexvar_new, varindex_old, varindex_new
   for(i = 0; i < nvars; ++i)
   {
      var = vars_array[i];
      //careful: hashmapindex+1, because '0' is treated as an empty hashmap entry, which causes an error
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( ! SCIPhashmapExists(detectordata->indexmap->indexvar, (void*) hashmapindex));
      SCIPhashmapInsert(detectordata->indexmap->indexvar, (void*) hashmapindex, (void*) var);
      assert( ! SCIPhashmapExists(detectordata->indexmap->varindex, (void*) var));
      SCIPhashmapInsert(detectordata->indexmap->varindex, (void*) var, (void*) hashmapindex);
   }
   //initialize hash maps for constraints: indexcons_old, indexcons_new, consindex_old, consindex_new
   for(i = 0; i < ncons; ++i)
   {
      cons = cons_array[i];
      //careful: i+1, because '0' is treated as an empty hashmap entry, which causes an error
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( ! SCIPhashmapExists(detectordata->indexmap->indexcons, (void*) hashmapindex));
      SCIPhashmapInsert(detectordata->indexmap->indexcons, (void*) hashmapindex, (void*) cons);
      assert( ! SCIPhashmapExists(detectordata->indexmap->consindex, (void*) cons));
      SCIPhashmapInsert(detectordata->indexmap->consindex, (void*) cons, (void*) hashmapindex);
   }
   SCIPdebugMessage("initializing hash maps DONE.\n");
   SCIPdebugMessage("initializing index arrays...\n");
#ifndef NDEBUG
   {
      //remove '/' from problem name
      const char* pname;
      char filename[256];
      pname = strrchr(SCIPgetProbName(scip), '/');
      if( pname == NULL )
      {
         pname = SCIPgetProbName(scip);
      }
      else
      {
         pname = pname+1;
      }
      sprintf(filename, "%s_initial_problem", pname);
      plotInitialProblem(scip, detectordata, filename);
   }
#endif
   //initialize index arrays ibegin, iend, jbegin, jend
   rowindices = rowindices_list(scip, detectordata, detectordata->indexmap->indexcons, detectordata->indexmap->varindex);
   columnindices = columnindices_list(scip, detectordata, rowindices);
   formIndexArray(detectordata->ibegin, detectordata->iend, rowindices);
   formIndexArray(detectordata->jbegin, detectordata->jend, columnindices);

   indexmapCreate(scip, &indexmap_permuted, ncons, nvars);
   for(i = 0; i < nvars; ++i)
   {
      var = vars_array[i];
      //careful: hashmapindex+1, because '0' is treated as an empty hashmap entry, which causes an error
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( ! SCIPhashmapExists(indexmap_permuted->indexvar, (void*) hashmapindex));
      SCIPhashmapInsert(indexmap_permuted->indexvar, (void*) hashmapindex, (void*) var);
      assert( ! SCIPhashmapExists(indexmap_permuted->varindex, (void*) var));
      SCIPhashmapInsert(indexmap_permuted->varindex, (void*) var, (void*) hashmapindex);
   }
   //debug: temporary hash maps: how to initialize???
   for(i = 0; i < ncons; ++i)
   {
      cons = cons_array[i];
      //careful: i+1, because '0' is treated as an empty hashmap entry, which causes an error
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( ! SCIPhashmapExists(indexmap_permuted->indexcons, (void*) hashmapindex));
      SCIPhashmapInsert(indexmap_permuted->indexcons, (void*) hashmapindex, (void*) cons);
      assert( ! SCIPhashmapExists(indexmap_permuted->consindex, (void*) cons));
      SCIPhashmapInsert(indexmap_permuted->consindex, (void*) cons, (void*) hashmapindex);
   }
   SCIP_CALL(SCIPallocMemoryArray(scip, &ibegin_permuted, ncons));
   SCIP_CALL(SCIPallocMemoryArray(scip, &iend_permuted, ncons));
   SCIP_CALL(SCIPallocMemoryArray(scip, &jbegin_permuted, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &jend_permuted, nvars));

   //ROC2 algorithm
   SCIPdebugMessage("starting ROC2 algortihm\n");
   i = 0;
   do
   {
      ++i;
      SCIPdebugMessage("Iteration # %i of ROC2\n", i);
      rankOrderClustering(scip, detectordata, detectordata->indexmap, indexmap_permuted);
      //form the new index arrays after the permutation
      SCIPlistDeleteNested(scip, rowindices);
      SCIPlistDeleteNested(scip, columnindices);
      rowindices = rowindices_list(scip, detectordata, indexmap_permuted->indexcons, indexmap_permuted->varindex);
      columnindices = columnindices_list(scip, detectordata, rowindices);
      formIndexArray(ibegin_permuted, iend_permuted, rowindices);
      formIndexArray(jbegin_permuted, jend_permuted, columnindices);
      //switch between index arrays containing new and old indices
      switchPointers( (void*) &detectordata->ibegin, (void*) &ibegin_permuted);
      switchPointers( (void*) &detectordata->iend, (void*) &iend_permuted);
      switchPointers( (void*) &detectordata->jbegin, (void*) &jbegin_permuted);
      switchPointers( (void*) &detectordata->jend, (void*) &jend_permuted);
      //switch between hash maps containing new and old indices
      switchPointers( (void*) &detectordata->indexmap, (void*) &indexmap_permuted);
   }
   while( ! (arraysAreEqual(detectordata->ibegin, ibegin_permuted, ncons)
          && arraysAreEqual(detectordata->ibegin, ibegin_permuted, ncons)
          && arraysAreEqual(detectordata->ibegin, ibegin_permuted, ncons)
          && arraysAreEqual(detectordata->ibegin, ibegin_permuted, ncons)));
   ROC_iterations = i;
#ifndef NDEBUG
   {
      //remove '/' from problem name
      const char* pname;
      char filename[256];
      pname = strrchr(SCIPgetProbName(scip), '/');
      if( pname == NULL )
      {
         pname = SCIPgetProbName(scip);
      }
      else
      {
         pname = pname+1;
      }
      sprintf(filename, "%s_ROC", pname);
      plotInitialProblem(scip, detectordata, filename);
   }
#endif
   //check conditions for arrays ibegin and jbegin: ibegin[i]<=ibegin[i+k] for all positive k
#ifdef SCIP_DEBUG
   for(i = 0; i < detectordata->nRelevantConss - 1; ++i)
   {
      assert(detectordata->ibegin[i] <= detectordata->ibegin[i+1]);
   }
   for(i = 0; i < nvars - 1; ++i)
   {
      assert(detectordata->jbegin[i] <= detectordata->jbegin[i+1]);
   }
#endif
   //arrays jmin, jmax and minV
   SCIPdebugMessage("calculating index arrays\n");
   detectordata->jmin[0] = detectordata->ibegin[0];
   detectordata->jmax[0] = detectordata->iend[0];
   detectordata->width[0] = detectordata->iend[0] - detectordata->ibegin[0];
   for(i = 1; i < ncons; ++i)
   {
      detectordata->width[i] = detectordata->iend[i] - detectordata->ibegin[i];
      detectordata->jmin[i] = detectordata->ibegin[i];
      detectordata->jmax[i] = maximum(detectordata->iend[i], detectordata->jmax[i-1]);
      detectordata->minV[i-1]=1 + (detectordata->jmax[i-1] - detectordata->jmin[i]);
   }
   n = maxArray(detectordata->width, ncons);
   v = minArray(detectordata->width, ncons);
   tau = round((nvars - v)/(n - v));
   if(tau > detectordata->maxblocks)
   {
      tau = detectordata->maxblocks;
   }

#ifndef NDEBUG
   SCIPdebugMessage("<N> <n> <v> <tau>: <%i> <%i> <%i> <%i>\n", nvars, n, v, tau);
//   printf("minV = [ ");
//   for(i = 0; i < ncons-1; ++i)
//   {
//      printf("%i ", detectordata->minV[i]);
//   }
//   printf("]\n");
//   printf("jmin = [ ");
//   for(i = 0; i < ncons; ++i)
//   {
//      printf("%i ", detectordata->jmin[i]);
//   }
//   printf("]\n");
//   printf("jmax = [ ");
//   for(i = 0; i < ncons; ++i)
//   {
//      printf("%i ", detectordata->jmax[i]);
//   }
//   printf("]\n");
#endif
   //continue only if tau >= 2
   if(tau < 2)
   {
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   rowsWithConstriction(scip, detectordata);
//   SCIPdebugMessage("rows with constriction:");
//   SCIPlistPrint(detectordata->rowsWithConstrictions, printint); //debug

   //BLOCKING
   blocking(scip, detectordata, tau, nvars);
   SCIPdebugMessage("BLOCKING DONE.\n");

   //debug plot the blocking  plot for [i=1:2:1] 'test.dat' every :::i::i lt i pt 5
#ifndef NDEBUG
   {
   //remove '/' from problem name
   const char* pname;
   char filename1[256];
   char filename2[256];
   char paramfile[256];
   pname = strrchr(SCIPgetProbName(scip), '/');
   if( pname == NULL )
   {
      pname = SCIPgetProbName(scip);
   }
   else
   {
      pname = pname+1;
   }
   sprintf(filename1, "%s_blocking", pname);
   sprintf(filename2, "%s_minV", pname);
   sprintf(paramfile, "%s.params", pname);
   SCIPdebugMessage("plotBlocking.\n");
   plotBlocking(scip, detectordata, filename1);
   SCIPdebugMessage("plotMinV.\n");
   plotMinV(scip, detectordata, filename2);
   SCIP_CALL(SCIPstopClock(scip, detectordata->clock));
   SCIPdebugMessage("writing Params.\n");
   writeParams(scip, detectordata, paramfile, ROC_iterations, tau, SCIPgetClockTime(scip, detectordata->clock));
   }
#endif
   //fill detectordata for a single block for testing
//   detectordata->blocks = 1;
//   detectordata->varsperblock[0] = vars_array;
//   detectordata->nvarsperblock[0] = nvars;
//   detectordata->consperblock[0] = SCIPgetConss(scip);
//   detectordata->nconsperblock[0] = ncons;
//   detectordata->linkingconss = NULL;
//   detectordata->nlinkingconss = 0;
   detectordata->found = TRUE;
   SCIP_CALL(copyDetectorDataToDecomp(scip, detectordata, detectordata->decdecomp));

   //deallocate memory
   SCIPlistDeleteData(scip, detectordata->rowsWithConstrictions);
//   SCIPlistDeleteData(scip, detectordata->blockedAfterrow);
   SCIPlistDeleteNested(scip, rowindices);
   SCIPlistDeleteNested(scip, columnindices);
   indexmapFree(scip, indexmap_permuted);

   SCIPfreeMemoryArray(scip, &ibegin_permuted);
   SCIPfreeMemoryArray(scip, &iend_permuted);
   SCIPfreeMemoryArray(scip, &jbegin_permuted);
   SCIPfreeMemoryArray(scip, &jend_permuted);

   //debug Ein Hack von Martin
   SCIP_CALL(SCIPallocMemoryArray(scip, decdecomps, 1));
   (*decdecomps)[0] = detectordata->decdecomp;
   *ndecdecomps = 1;
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

static
DEC_DECL_GETPRIORITY(getPriority)
{
   DEC_DETECTOR* arrowheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   arrowheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(arrowheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(arrowheur), DEC_DETECTORNAME) == 0);
   return detectordata->priority;
}

/** creates the stairheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionStairheur(
   SCIP*                 scip              /**< SCIP data structure */

   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL(SCIPallocMemory(scip, &detectordata));

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->decdecomp = NULL;
   SCIP_CALL(DECincludeDetector(scip, DEC_DETECTORNAME, detectordata, detectAndBuildStair, initStairheur, exitStairheur, getPriority));

   /* add stairheur presolver parameters */
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/priority", "Priority of detector", &detectordata->priority, FALSE, DEFAULT_PRIORITY, INT_MIN, INT_MAX, NULL, NULL));
   return SCIP_OKAY;
}
