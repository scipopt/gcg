#include "type_bliss.h"

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
   )
{
#ifdef BLISS_PATCH_PRESENT
   graph.set_search_limits(searchnodelimit, generatorlimit);
#endif

#if BLISS_VERSION_MAJOR >= 1 || BLISS_VERSION_MINOR >= 76
   auto report = [&](unsigned int n, const unsigned int* aut) {
      fhook((void*)ptrhook, n, aut);
   };

   auto term = [&]() {
      return (generatorlimit > 0 && bstats.get_nof_generators() >= (long unsigned int) generatorlimit);
   };

   graph.find_automorphisms(bstats, report, term);
#else
   graph.find_automorphisms(bstats, fhook, ptrhook);
#endif
}
