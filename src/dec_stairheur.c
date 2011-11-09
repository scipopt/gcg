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
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "dec_stairheur.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "scip_misc.h"
#include "scip/struct_var.h"

#define DEC_DETECTORNAME      "stairheur"   /**< name of the detector */
#define DEC_PRIORITY          -100              /**< priority of the detector */

/* Default parameter settings*/
#define DEFAULT_BLOCKS                    2     /**< number of blocks */
#define DEFAULT_CONSWEIGHT                5     /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
#define DEFAULT_TIDY                      TRUE  /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES	              0.2   /**< percentage of dummy vertices*/

#define DEFAULT_MAXBLOCKS                 10    /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                 2     /**< value for the minimum number of blocks to be considered */

#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_PRIORITY                  DEC_PRIORITY

#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)

#define GP_NAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_%d.gp", (name), (blocks), (cons), (dummy)


/*include the code for doubly linked lists */
#include <stdio.h>
#include <stdlib.h>

/** TODO:
 * deallocate all memory (lists, iterators, data stored in lists...)
 *
 * */

typedef enum { false, true } bool;

struct node {
 void* data;
 struct node* prev;
 struct node* next;
};
typedef struct node NODE;

struct list {
   /** a doubly linked list
    * the first and last node are referred to by nil->next and nil->prev, respectively.
    */
   NODE* nil;
   int size;
};
typedef struct list LIST;

struct iterator {
   NODE* node;
   LIST* list;
};
typedef struct iterator ITERATOR;

struct reverse_iterator {
   NODE* node;
   LIST* list;
};
typedef struct reverse_iterator REVERSE_ITERATOR;

LIST* list_create_empty(void);
LIST* list_create(void* data);
LIST* list_int_create(int from, int to);
LIST* list_flat_copy(LIST* list);
bool list_push_front(LIST* list, void* data);
bool list_pop_front(LIST* list);
bool list_push_back(LIST* list, void* data);
bool list_pop_back(LIST* list);
NODE* list_insert(ITERATOR* it, void* data);
bool list_assign(ITERATOR* it, void* data);
bool list_assign_list(LIST* list1, LIST* list2);
bool list_empty(LIST* list);
bool list_erase(ITERATOR* it);
bool list_delete(LIST* list);
bool list_delete_data(LIST* list);
ITERATOR list_find(ITERATOR first, ITERATOR last, void* value, bool (*comp_func)(void* a, void* b));
bool list_move(ITERATOR position, ITERATOR it);
void list_move_front(ITERATOR position);
void list_move_back(ITERATOR position);
SCIP_RETCODE list_rearrange(LIST* list, LIST* order);
bool list_foreach(LIST* list, int(*func)(void*));

NODE* node_create(void* data);
int printstring(void *s);
int printint(void *i);
void list_prin(LIST* list, int(*func)(void*));
static int compare (const void * a, const void * b);
static bool compare_int(void* a, void * b);

ITERATOR* iterator_create(void);
ITERATOR iterator_begin(LIST* list);
ITERATOR iterator_end(LIST* list);
bool iterator_next(ITERATOR* it);
bool iterator_prev(ITERATOR* it);
void iterator_equals(ITERATOR it1, ITERATOR it2);
bool iterator_is_equal(ITERATOR it1, ITERATOR it2);

//REVERSE_ITERATOR* reverse_iterator_create(void);
//REVERSE_ITERATOR reverse_iterator_begin(LIST* list);
//REVERSE_ITERATOR reverse_iterator_end(LIST* list);
//bool reverse_iterator_next(REVERSE_ITERATOR* it);
//bool reverse_iterator_prev(REVERSE_ITERATOR* it);
//void reverse_iterator_equals(REVERSE_ITERATOR it1, REVERSE_ITERATOR it2);
//bool reverse_iterator_is_equal(REVERSE_ITERATOR it1, REVERSE_ITERATOR it2);






LIST* list_create_empty()
{
   LIST* list;
   NODE* node;
   if(! (list = malloc(sizeof(LIST)))) return NULL;
   //debugging
//   int* i;
//   i = malloc(sizeof(int));
//   *i = 0;
   //node = node_create(i);
   node = node_create(NULL);
   list->nil = node;
   list->nil->next = list->nil;
   list->nil->prev = list->nil;
   list->nil->data = node->data;
   list->size = 0;
   return list;
}

LIST* list_create(void* data)
{
   LIST* list;
   list = list_create_empty();
   list_push_front(list, data);
   return list;
}

LIST* list_int_create(int from, int to)
{
   LIST* list;
   int* data;
   int i;
   list = list_create_empty();
   for(i = from; i <= to; ++i)
   {
      data = (int*) malloc(sizeof(int));
      *data=i;
      list_push_back(list, (void*) data);
   }
   return list;
}

LIST* list_flat_copy(LIST* list)
{
   /**creates a copy of list, but does not copy the data stored in list->data */
   LIST* copy;
   ITERATOR it;
   //case: list exists
   if(list)
   {
      copy = list_create_empty();
      for(it = iterator_begin(list); ! iterator_is_equal(it, iterator_end(list)); iterator_next(&it))
      {
         list_push_back(copy, it.node->data);
      }
      return copy;
   }
   //case: list does not exist
   else
   {
      return NULL;
   }
}

bool list_push_front(LIST* list, void* data)
{
   /** adds a node at the front of the list */
   NODE* newnode;
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return false;
   newnode = node_create(data);
   //it points to the beginning of list
   it = iterator_begin(list);
   list_insert(&it, data);
   return true;
}

bool list_pop_front(LIST* list)
{
   /** removes the node at the front of the list */
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (list_empty(list)) ) return false;

   it = iterator_begin(list);
   list_erase(&it);
   return true;
}

bool list_push_back(LIST* list, void* data)
{
   /** adds a node at the back of the list */
   NODE* newnode;
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return false;
   newnode = node_create(data);
   //it points to the beginning of list
   it = iterator_end(list);
   list_insert(&it, data);
   return true;
}

bool list_pop_back(LIST* list)
{
   /** removes the node at the back of the list */
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (list_empty(list)) ) return false;

   it = iterator_end(list);
   iterator_prev(&it);
   list_erase(&it);
   return true;
}

NODE* list_insert(ITERATOR* it, void* data)
{
   /** creates a new node before 'it->node' and returns the new node*/
   NODE* newnode;
   NODE* successor;
   NODE* predecessor;
   //case: 'it' points to no node or list
   if(it->node == NULL && it->list == NULL) return NULL;
   newnode=node_create(data);
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

bool list_assign(ITERATOR* it, void* data)
{
   /** assign new data to 'it->node' */
   if(it && it->node)
   {
      it->node->data=data;
      return true;
   }
   else
   {
      return false;
   }
}


bool list_assign_list(LIST* target, LIST* origin)
{
   /**if target and origin have the same size, this function assigns the data pointers of origin to target */
   ITERATOR it1;
   ITERATOR it2;
   if(target && origin && target->size == origin->size)
   {
      for(it1 = iterator_begin(target), it2 = iterator_begin(origin); ! iterator_is_equal(it1, iterator_end(target)); iterator_next(&it1), iterator_next(&it2))
      {
         list_assign(&it1, it2.node->data);
      }
      return true;
   }
   else
   {
      return false;
   }
}

NODE* node_create(void* data)
{
   /** creates a node with a void pointer to data and returns the node */
   NODE* node;
   if(! (node = malloc(sizeof(NODE)))) return NULL;
   node->data = data;
   node->prev = NULL;
   node->next = NULL;
   return node;
}

bool list_erase(ITERATOR* it)
{
   /** removes the node 'it->node' from the list 'it->list'
    * returns true if 'node' is in 'list'
    * returns false if 'node' is not in 'list'
    * the memory the data pointer points to is not deallocated*/
   NODE* successor;
   NODE* predecessor;
   //case: iterator, node or list do not exist or empty list
   if( (! it) || (! it->node) || (! it->list) || list_empty(it->list))
   {
      return false;
   }
   successor = it->node->next;
   predecessor = it->node->prev;
   predecessor->next = successor;
   successor->prev = predecessor;
   free(it->node);
   it->node = successor;
   --it->list->size;
   return true;
}

bool list_delete(LIST* list)
{
   /** deletes the entire list
    *  for deallocating memory of the list data, call 'list_deleta_data' before */
   if(! list)
   {
      return false;
   }
   else
   {
      while(! list_empty(list))
      {
         list_pop_back(list);
      }
      //deallocate sentinel node
      free(list->nil);
      free(list);
      list = NULL;
      return true;
   }
}

bool list_delete_data(LIST* list)
{
   /** deallocates the memory the data pointers of the list points to */
   ITERATOR it;
   if(! list)
   {
      return false;
   }
   else
   {
      for( it = iterator_begin(list); ! ( iterator_is_equal(it, iterator_end(list)) ); iterator_next(&it))
      {
         free(it.node->data);
      }
      return true;
   }
}

bool list_empty(LIST* list)
{
   /** returns true iff list exists and has 0 nodes (only sentinel node)*/
//   if( list && list->size == 0)
   //list is empty if the next pointer of the sentinel points to itself
   if( list && list->nil->next == list->nil)
   {
      return true;
   }
   else
   {
      return false;
   }
}

bool list_foreach(LIST* list, int(*func)(void*))
{
   /** runs the function 'func' for every node of the list 'node' */
//   NODE* node = list->first;
//   while(node) {
//      if(func(node->data)!=0) return false;
//      node=node->next;
//   }
//   return true;

   ITERATOR it;
   for(it = iterator_begin(list); ! ( iterator_is_equal(it, iterator_end(list)) ); iterator_next(&it))
   {
      if(func(it.node->data) != 0)
      {
         return false;
      }
   }
   return true;
}

void list_prin(LIST* list, int(*func)(void*))
{
   printf("( ");
   list_foreach(list, func);
   printf(")\n");
}

ITERATOR list_find(ITERATOR first, ITERATOR last, void* value, bool (*comp_func)(void* a, void* b))
{
   for (; ! (iterator_is_equal(first, last)); iterator_next(&first))
   {
      if ( comp_func(value, first.node->data) )
      {
         break;
      }
   }
   return first;
}

bool list_move(ITERATOR position, ITERATOR it)
{
   /** moves the element at 'position' to the location in front of element 'it'
    * both elements must be in the same list
    */
   NODE* old_successor;
   NODE* old_predecessor;
   NODE* new_successor;
   NODE* new_predecessor;
   if(position.node && it.node && position.list == it.list)
   {
      //case 'it' and 'position' point to the same node or 'it' points to the succeeding node of 'position': no moving
      if(position.node == it.node || position.node->next == it.node)
      {
         return true;
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

      return true;
   }
   else
   {
      return false;
   }
}

void list_move_front(ITERATOR position)
{
   /** moves the element at 'position' to the front of the list */
   list_move(position, iterator_begin(position.list));
}

void list_move_back(ITERATOR position)
{
   /** moves the element at 'position' to the end of the list */
   list_move(position, iterator_end(position.list));
}

SCIP_RETCODE list_rearrange(LIST* list, LIST* order)
{
   /** rearranges elements of list according to the ordering of order
    * example: list = (a b c d); order = (3 2 4 1)
    * after calling list_rearrange(list, order): list = (c b d a)
    * both lists must have same size!
    * order must have elements from 1 to list->size
    */
   LIST* new_list;
   ITERATOR it1, it2;
   int i;
   if( list && order && list->size == order->size)
   {
      new_list = list_create_empty();
      for(it1 = iterator_begin(order); ! iterator_is_equal(it1, iterator_end(order)); iterator_next(&it1))
      {
         for( it2 = iterator_begin(list), i = 1; i < *(int*)it1.node->data; ++i)
         {
            iterator_next(&it2);
         }
         list_push_back(new_list, it2.node->data);
      }
      list_assign_list(list, new_list);
      list_delete(new_list);
      return SCIP_OKAY;
   }
   else
   {
      return SCIP_ERROR;
   }
}

int printstring(void *s)
{
   printf("%s\n", (char *)s);
   return 0;
}

int printint(void *i)
{
//   int* n;
//   n=(int*) i;
//   printf("%i\n", *n);
//   return 0;
   printf("%i ", *(int*)i);
   return 0;
}

//forward iterators
ITERATOR* iterator_create()
{
   ITERATOR* it;
   if(! (it = malloc(sizeof(ITERATOR)))) return NULL;
   it->list = NULL;
   it->node = NULL;
   return it;
}

ITERATOR iterator_begin(LIST* list)
{
   ITERATOR it = {NULL, NULL};
   if(list)
   {
      it.list = list;
      it.node = list->nil->next;
   }
   return it;
}

ITERATOR iterator_end(LIST* list)
{
   ITERATOR it = {NULL, NULL};
   if(list)
   {
      it.list = list;
      it.node = list->nil;
   }
   return it;
}


bool iterator_next(ITERATOR* it)
{
   /** sets it to the next element of the list
    */
   if(it && it->list && it->node)
   {
      it->node = it->node->next;
      return true;
   }
   else
   {
      return false;
   }
}

bool iterator_prev(ITERATOR* it)
{
   /** sets it to the previous element of the list
    */
   if(it && it->list && it->node)
   {
      it->node = it->node->prev;
      return true;
   }
   else
   {
      return false;
   }
}

void iterator_equals(ITERATOR it1, ITERATOR it2)
{
   /** assigns it1 = it2 */
   it1.list = it2.list;
   it1.node = it2.node;
}

bool iterator_is_equal(ITERATOR it1, ITERATOR it2)
{
   if(it1.list == it2.list && it1.node == it2.node)
   {
      return true;
   }
   else
   {
      return false;
   }
}

////reverse iterators
//REVERSE_ITERATOR* reverse_iterator_create()
//{
//   REVERSE_ITERATOR* rit;
//   if(! (rit = malloc(sizeof(REVERSE_ITERATOR)))) return NULL;
//   rit->list = NULL;
//   rit->node = NULL;
//   return rit;
//}
//
//REVERSE_ITERATOR reverse_iterator_begin(LIST* list)
//{
//   REVERSE_ITERATOR rit = {NULL, NULL};
//   if(list)
//   {
//      rit.list = list;
//      rit.node = list->nil->prev;
//   }
//   return rit;
//}
//
//REVERSE_ITERATOR reverse_iterator_end(LIST* list)
//{
//   REVERSE_ITERATOR rit = {NULL, NULL};
//   if(list)
//   {
//      rit.list = list;
//      rit.node = list->nil;
//   }
//   return rit;
//}
//
//
//bool reverse_iterator_next(REVERSE_ITERATOR* rit)
//{
//   /** sets rit to the next element of the list
//    */
//   if(rit && rit->list && rit->node)
//   {
//      rit->node = rit->node->prev;
//      return true;
//   }
//   else
//   {
//      return false;
//   }
//}
//
//bool reverse_iterator_prev(REVERSE_ITERATOR* rit)
//{
//   /** sets rit to the previous element of the list
//    */
//   if(rit && rit->list && rit->node)
//   {
//      rit->node = rit->node->next;
//      return true;
//   }
//   else
//   {
//      return false;
//   }
//}
//
//void reverse_iterator_equals(REVERSE_ITERATOR rit1, REVERSE_ITERATOR rit2)
//{
//   /** assigns rit1 = rit2 */
//   rit1.list = rit2.list;
//   rit1.node = rit2.node;
//}
//
//bool reverse_iterator_is_equal(REVERSE_ITERATOR rit1, REVERSE_ITERATOR rit2)
//{
//   if(rit1.list == rit2.list && rit1.node == rit2.node)
//   {
//      return true;
//   }
//   else
//   {
//      return false;
//   }
//}

/*
 * Data structures
 */
///** score data structure **/
//struct Dec_StairheurScores
//{
//   SCIP_Real borderscore;
//   SCIP_Real minkequicutscore;
//   SCIP_Real equicutscorenormalized;
//   SCIP_Real densityscore;
//   SCIP_Real linkingscore;
//};
typedef struct Dec_StairheurScores DEC_STAIRHEURSCORES;

/** detector data */
struct DEC_DetectorData
{
   DECDECOMP* decdecomp;
   SCIP_VAR*** varsperblock;
   int* nvarsperblock;
   SCIP_CONS*** consperblock;
   int *nconsperblock;
   SCIP_CONS** linkingconss;
   int nlinkingconss;
   int blocks;
   int maxblocks;
   int minblocks;
   SCIP_HASHMAP*  probindex_to_var;
   SCIP_HASHMAP*  varindex;
   SCIP_HASHMAP*  consindex;
   SCIP_Bool found;
   int priority;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

static
bool compare_int(void* a, void * b)
{
   return ( *(int*)a == *(int*)b ? true : false );
}


static
LIST* rowindices_list(SCIP* scip)
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
   SCIP_CONS** cons_array; //array of constraints of the problem
   SCIP_CONS* cons; //one constraint of the problem
   SCIP_VAR** consvars; //array of variables that occur in a constraint (unequal zero)

   rowindices = list_create_empty();
   cons_array = SCIPgetConss(scip);
   ncons = SCIPgetNConss(scip);
   for(i = 0; i < ncons; ++i)
   {
      cons = cons_array[i];
      nconsvars = SCIPgetNVarsXXX(scip, cons);
      consvars = SCIPgetVarsXXX(scip, cons);
      //allocate memory for the array of probindices
      probindices=(int*)malloc(nconsvars*sizeof(int));
      //fill the array with the indices of the variables of the current constraint
      for(j = 0; j < nconsvars; ++j)
      {
         probindices[j] = SCIPvarGetProbindex(consvars[j])+1;
      }
      //sort the elements of probindices ('<')
      qsort(probindices, nconsvars, sizeof(int), compare);
      //store a copy of the elements of probindices in the list rowindices_row
      rowindices_row = list_create_empty();
      for(j = 0; j < nconsvars; ++j)
      {
         data = (int*) malloc(sizeof(int));
         *data = probindices[j];
         list_push_back(rowindices_row, data);
      }
      //deallocate memory
      free(probindices);
      //add rowindices_row to the list rowindices
      list_push_back(rowindices, rowindices_row);
   }
   return rowindices;
}

static
LIST* columnindices_list(SCIP* scip, LIST* rowindices)
{
   /** creates a nested list which contains the position of non-zero entries of the boundary constraints */
   LIST** columnindices_array;
   LIST* columnindices;
   LIST* rowindices_row;
   int* data;
   int position;
   int nvars;
   int ncons;
   int i;
   ITERATOR it1;
   ITERATOR it2;
   nvars = SCIPgetNVars(scip);
   ncons = SCIPgetNConss(scip);
   //create the columnindices_array with pointers to empty lists
   columnindices_array = malloc(nvars * sizeof(LIST*));
   for(i = 0; i < nvars; ++i)
   {
      columnindices_array[i] = list_create_empty();
   }

   for(it1 = iterator_begin(rowindices), i = 0; ! ( iterator_is_equal(it1, iterator_end(rowindices)) ); iterator_next(&it1), ++i)
   {
      rowindices_row = it1.node->data;
      for(it2 = iterator_begin(rowindices_row); ! ( iterator_is_equal(it2, iterator_end(rowindices_row)) ); iterator_next(&it2))
      {
         data = (int*) malloc(sizeof(int));
         *data = i+1;
         position = *(int*)(it2.node->data)-1;
         list_push_back(columnindices_array[position], data);
      }
   }
   //create a columnindices list instead of an array
   columnindices = list_create_empty();
   for(i = 0; i < nvars; ++i)
   {
      list_push_back(columnindices, columnindices_array[i]);
   }
   //deallocate memory
   free(columnindices_array);
   return columnindices;
}

static
LIST* row_ordering(LIST* roworder, LIST* columnindices)
{
   /** does the row ordering of the ROC2 algortihm */
   LIST* new_roworder;
   ITERATOR it1;
   ITERATOR it2;
   ITERATOR it3;
   ITERATOR it4;
   new_roworder = list_flat_copy(roworder);

   for( it1 = iterator_end(columnindices), iterator_prev(&it1); it1.node != it1.list->nil; iterator_prev(&it1))
   {
      for( it2 = iterator_end(roworder), iterator_prev(&it2); it2.node != it2.list->nil; iterator_prev(&it2) )
      {
         it3 = list_find(iterator_begin(it1.node->data), iterator_end(it1.node->data), it2.node->data, compare_int);
         if(it3.node->data != NULL)
         {
            it4 = list_find(iterator_begin(new_roworder), iterator_end(new_roworder), it2.node->data, compare_int);
            list_move_front(it4);
         }
         else
         {
            ;
         }
      }
      list_assign_list(roworder, new_roworder);
   }
   //deallocate memory
   list_delete(new_roworder);
   return roworder;
}

static
SCIP_RETCODE RankOrderClustering(
      SCIP*             scip,       /**< SCIP data structure */
      DEC_DETECTORDATA*  detectordata /**< presolver data data structure */
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
   SCIP_CONS** cons_array;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(detectordata != NULL);

   nvars = SCIPgetNVars(scip);
   ncons = SCIPgetNConss(scip);

   printf("ncons: %i \n", ncons); //debug
   printf("nvars: %i \n", nvars); //debug
   //create a list for the order of the rows ( 1 2 3 ... ncons ) and columns ( 1 2 3 ... nvars )
   roworder = list_int_create(1, ncons);
   columnorder = list_int_create(1, nvars);

   //create the rowindices list from the scip object;
   rowindices = rowindices_list(scip);
//   printf("rowindices:\n");
//   for(it1 = iterator_begin(rowindices); ! ( iterator_is_equal(it1, iterator_end(rowindices)) ); iterator_next(&it1))
//   {
//      list_prin(it1.node->data, printint);
//   }
   //create the columnindices list
   columnindices = columnindices_list(scip, rowindices);

   //print for debugging
//   printf("columnorder: ");
//   list_prin(columnorder, printint);
//   printf("roworder: ");
//   list_prin(roworder, printint);
//   printf("rowindices:\n");
//   for(it1 = iterator_begin(rowindices); ! ( iterator_is_equal(it1, iterator_end(rowindices)) ); iterator_next(&it1))
//   {
//      list_prin(it1.node->data, printint);
//   }
//   printf("columnindices:\n");
//   for(it1 = iterator_begin(columnindices); ! ( iterator_is_equal(it1, iterator_end(columnindices)) ); iterator_next(&it1))
//   {
//      list_prin(it1.node->data, printint);
//   }

   //row ordering
   roworder = row_ordering(roworder, columnindices);
   //rearrange the rowindices list according to the new rowordering
   list_rearrange(rowindices, roworder);
   //column ordering
   columnorder = row_ordering(columnorder, rowindices);

   printf("row ordering:");
   list_prin(roworder, printint);
   printf("column ordering:");
   list_prin(columnorder, printint);

//   printf("ROC #1\n"); //debug
   //store the row and column ordering in consindex and varindex of data structure detectordata
//   assert( SCIPhashmapRemoveAll(detectordata->consindex) );
//   assert( SCIPhashmapRemoveAll(detectordata->varindex) );
   cons_array = SCIPgetConss(scip);

   //consindex
   for(it1 = iterator_begin(roworder), i = 0; i < ncons; ++i, iterator_next(&it1))
   {
      position = (*(int*)it1.node->data)-1;
//      assert( ! SCIPhashmapExists(detectordata->consindex, cons_array[position]) );
//      SCIPhashmapInsert(detectordata->consindex, cons_array[position], (void*) (i+1));
      assert( SCIPhashmapExists(detectordata->consindex, cons_array[position]) );
      SCIPhashmapSetImage(detectordata->consindex, cons_array[position], (void*) (i+1));
   }
   //varindex
   for(it1 = iterator_begin(columnorder), i = 0; i < nvars; ++i, iterator_next(&it1))
   {
      position = (*(int*)it1.node->data);
      var = (SCIP_VAR*) SCIPhashmapGetImage(detectordata->probindex_to_var, (void*) position);
      assert ( var != NULL);
//      assert( ! SCIPhashmapExists(detectordata->varindex, (void*) var) );
//      SCIPhashmapInsert( detectordata->varindex, var, (void*) (i+1) );
      assert( SCIPhashmapExists(detectordata->varindex, (void*) var) );
      SCIPhashmapSetImage(detectordata->varindex,  var, (void*) (i+1));
   }


   //deallocate memory
   list_delete_data(roworder);
   list_delete(roworder);
   list_delete_data(columnorder);
   list_delete(columnorder);
   for(it1 = iterator_begin(rowindices); ! ( iterator_is_equal(it1, iterator_end(rowindices)) ); iterator_next(&it1))
   {
      list_delete_data(it1.node->data);
      list_delete(it1.node->data);
   }
   for(it1 = iterator_begin(columnindices);  ! ( iterator_is_equal(it1, iterator_end(columnindices)) ); iterator_next(&it1))
   {
      list_delete_data(it1.node->data);
      list_delete(it1.node->data);
   }
   list_delete(rowindices);
   list_delete(columnindices);

   return SCIP_OKAY;
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

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->linkingconss, detectordata->linkingconss, detectordata->nlinkingconss));
   decomp->nlinkingconss = detectordata->nlinkingconss;

   for( i = 0; i < detectordata->blocks; ++i )
   {
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipconss[i], detectordata->consperblock[i], detectordata->nconsperblock[i]));
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipvars[i], detectordata->varsperblock[i], detectordata->nvarsperblock[i]));
   }

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipconss, detectordata->nconsperblock, detectordata->blocks));
   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipvars, detectordata->nvarsperblock, detectordata->blocks));

   decomp->varindex = detectordata->varindex;
   decomp->consindex = detectordata->consindex;
   decomp->nblocks = detectordata->blocks;
   decomp->type = DEC_STAIRCASE;
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

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nconsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nvarsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->linkingconss, nconss));
   for( i = 0; i < detectordata->maxblocks; ++i )
   {
      detectordata->nvarsperblock[i] = 0;
      detectordata->nconsperblock[i] = 0;
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock[i], nconss));
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock[i], nvars));
   }
   detectordata->nlinkingconss = 0;
   /* create hash tables */
   SCIP_CALL(SCIPhashmapCreate(&detectordata->probindex_to_var, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->varindex, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->consindex, SCIPblkmem(scip), nconss));
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
   SCIPfreeMemoryArray(scip, &detectordata->linkingconss);
   SCIPfreeMemory(scip, &detectordata);
   return SCIP_OKAY;
}

static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildStair)
{
   int i;
   int probindex;
   int ncons; //number of constraints in the problem
   int nvars; //number of variables in the problem
   SCIP_VAR** vars_array;
   SCIP_VAR* var;
   SCIP_CONS** cons_array;
   DEC_DETECTOR* stairheur;
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);
   stairheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(stairheur);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(stairheur), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);

   nvars = SCIPgetNVars(scip);
   ncons = SCIPgetNConss(scip);
   vars_array = SCIPgetVars(scip);
   cons_array = SCIPgetConss(scip);

   //create a hashmap probindex_to_var: key: probindex, value: variable with this probindex
   for(i = 0; i < nvars; ++i)
   {
      var = vars_array[i];
      //careful: probindex+1, because '0' is treated as an empty hashmap entry, which causes an error
      probindex = (SCIPvarGetProbindex(var)+1);
      assert( ! SCIPhashmapExists(detectordata->probindex_to_var, (void*) probindex));
      SCIPhashmapInsert(detectordata->probindex_to_var, (void*) probindex, (void*) var);
   }
   //varindex, consindex
   //store the row and column ordering in consindex and varindex of data structure detectordata
   //consindex
   for(i = 0; i < ncons; ++i)
   {
      assert( ! SCIPhashmapExists(detectordata->consindex, cons_array[i]) );
      SCIPhashmapInsert(detectordata->consindex, cons_array[i], (void*) (i+1));
   }
   //varindex
   for(i = 0; i < nvars; ++i)
   {
      var = vars_array[i];
      probindex = SCIPvarGetProbindex(var);
      assert ( var != NULL);
      assert( ! SCIPhashmapExists(detectordata->varindex, (void*) var) );
      SCIPhashmapInsert( detectordata->varindex, var, (void*) (probindex + 1) );
   }

   //ROC2 algorithm
   for(i = 0; i < 1; ++i)
   {
      RankOrderClustering(scip, detectordata);
   }

   //fill detectordata for a single block for testing
   detectordata->blocks = 1;
   detectordata->varsperblock[0] = vars_array;
   detectordata->nvarsperblock[0] = nvars;
   detectordata->consperblock[0] = SCIPgetConss(scip);
   detectordata->nconsperblock[0] = ncons;
   detectordata->linkingconss = NULL;
   detectordata->nlinkingconss = 0;
   detectordata->found = TRUE;
   SCIP_CALL(copyDetectorDataToDecomp(scip, detectordata, detectordata->decdecomp));
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** set the decomp structure */
static
DEC_DECL_SETSTRUCTDECOMP(StairheurSetDecomp)
{
   DEC_DETECTOR* stairheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   stairheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(stairheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(stairheur), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Setting decdecomp\n");
   detectordata->decdecomp = decdecomp;
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
   SCIP_CALL(DECincludeDetector(scip, DEC_DETECTORNAME, detectordata, detectAndBuildStair, StairheurSetDecomp, initStairheur, exitStairheur, getPriority));

   /* add stairheur presolver parameters */
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/priority", "Priority of detector", &detectordata->priority, FALSE, DEFAULT_PRIORITY, INT_MIN, INT_MAX, NULL, NULL));
   return SCIP_OKAY;
}
