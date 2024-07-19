#include "symmetry/type_nauty.h"

void struct_graph::add_vertex(
   int id
   )
{

}

void struct_graph::add_edge(
   int v1,
   int v2
   )
{

}

unsigned int struct_graph::get_nof_vertices(
   )
{
   return 0;
}

void struct_graph::find_automorphisms(
   void* ptrhook,
   void (*fhook)(void*, unsigned int, const unsigned int*),
   unsigned int searchnodelimit,
   unsigned int generatorlimit
   )
{
   sparsenauty(&graph, lab, ptn, orbits, &options, &stats, NULL);
}
