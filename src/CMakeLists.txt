include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include(GNUInstallDirs)
#
# interface function for setting common library properties
#
function(setLibProperties targetname outputname)
    set_target_properties(${targetname} PROPERTIES
        OUTPUT_NAME ${outputname}
        MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endfunction(setLibProperties)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED on)

set(gcgsources
    gcg/benders_gcg.c
    gcg/bendersplugins.c
    gcg/branchgcg.c
    gcg/branch_bpstrong.c
    gcg/branch_compbnd.c
    gcg/branch_empty.c
    gcg/branch_generic.c
    gcg/branch_orig.c
    gcg/branch_relpsprob.c
    gcg/branch_ryanfoster.c
    gcg/class_conspartition.cpp
    gcg/class_detprobdata.cpp
    gcg/class_indexpartition.cpp
    gcg/class_partialdecomp.cpp
    gcg/class_pricingcontroller.cpp
    gcg/class_pricingtype.cpp
    gcg/class_stabilization.cpp
    gcg/class_varpartition.cpp
    gcg/clscons_consnamelevenshtein.cpp
    gcg/clscons_consnamenonumbers.cpp
    gcg/clscons_gamsdomain.cpp
    gcg/clscons_gamssymbol.cpp
    gcg/clscons_miplibconstypes.cpp
    gcg/clscons_nnonzeros.cpp
    gcg/clscons_scipconstypes.cpp
    gcg/clscons.c
    gcg/clsvar_gamsdomain.cpp
    gcg/clsvar_gamssymbol.cpp
    gcg/clsvar_objvalues.cpp
    gcg/clsvar_objvaluesigns.cpp
    gcg/clsvar_scipvartypes.cpp
    gcg/clsvar.c
    gcg/colpool.c
    gcg/cons_decomp.cpp
    gcg/cons_integralorig.c
    gcg/cons_masterbranch.c
    gcg/cons_origbranch.c
    gcg/dec_compgreedily.cpp
    gcg/dec_connected_noNewLinkingVars.cpp
    gcg/dec_connectedbase.cpp
    gcg/dec_consclass.cpp
    gcg/dec_constype.cpp
    gcg/dec_dbscan.cpp
    gcg/dec_densemasterconss.cpp
    gcg/dec_generalmastersetcover.cpp
    gcg/dec_generalmastersetpack.cpp
    gcg/dec_generalmastersetpart.cpp
    gcg/dec_mastersetcover.cpp
    gcg/dec_mastersetpack.cpp
    gcg/dec_mastersetpart.cpp
    gcg/dec_mst.cpp
    gcg/dec_neighborhoodmaster.cpp
    gcg/dec_postprocess.cpp
    gcg/dec_staircase_lsp.cpp
    gcg/dec_stairheur.cpp
    gcg/dec_varclass.cpp
    gcg/decomp.c
    gcg/dialog_explore.cpp
    gcg/dialog_gcg.c
    gcg/dialog_graph.cpp
    gcg/dialog_master.c
    gcg/disp_gcg.c
    gcg/disp_master.c
    gcg/event_bestsol.c
    gcg/event_display.c
    gcg/event_mastersol.c
    gcg/event_relaxsol.c
    gcg/event_mastersepacut.c
    gcg/event_solvingstats.c
    gcg/gcg.cpp
    gcg/gcg_general.c
    gcg/gcgcgcut.c
    gcg/gcgcol.c
    gcg/gcggithash.c
    gcg/gcgheur.c
    gcg/gcgplugins.c
    gcg/gcgpqueue.c
    gcg/gcgsepa.c
    gcg/gcgsort.c
    gcg/gcgvar.c
    gcg/gcgvarhistory.c
    gcg/heur_gcgcoefdiving.c
    gcg/heur_gcgdins.c
    gcg/heur_gcgfeaspump.c
    gcg/heur_gcgfracdiving.c
    gcg/heur_gcgguideddiving.c
    gcg/heur_gcglinesdiving.c
    gcg/heur_gcgpscostdiving.c
    gcg/heur_gcgrens.c
    gcg/heur_gcgrins.c
    gcg/heur_gcgrounding.c
    gcg/heur_gcgshifting.c
    gcg/heur_gcgsimplerounding.c
    gcg/heur_gcgveclendiving.c
    gcg/heur_gcgzirounding.c
    gcg/heur_greedycolsel.c
    gcg/heur_ipcolgen.c
    gcg/heur_mastercoefdiving.c
    gcg/heur_masterdiving.c
    gcg/heur_masterfracdiving.c
    gcg/heur_masterlinesdiving.c
    gcg/heur_mastervecldiving.c
    gcg/heur_origdiving.c
    gcg/heur_relaxcolsel.c
    gcg/heur_restmaster.c
    gcg/heur_setcover.c
    gcg/heur_xpcrossover.c
    gcg/heur_xprins.c
    gcg/mastersepacut.c
    gcg/extendedmasterconsdata.c
    gcg/masterplugins.c
    gcg/misc.c
    gcg/miscvisualization.cpp
    gcg/nodesel_master.c
    gcg/objdialog.cpp
    gcg/params_visu.c
    gcg/presol_roundbound.c
    gcg/pricer_gcg.cpp
    gcg/pricestore_gcg.c
    gcg/pricingcb.c
    gcg/pricingjob.c
    gcg/pricingprob.c
    gcg/reader_blk.cpp
    gcg/reader_cls.cpp
    gcg/reader_dec.cpp
    gcg/reader_gp.cpp
    gcg/reader_ref.c
    gcg/reader_tex.cpp
    gcg/relax_gcg.c
    gcg/scip_misc.c
    gcg/score_bender.cpp
    gcg/score_border.cpp
    gcg/score_classic.cpp
    gcg/score_fawh.cpp
    gcg/score_forswh.cpp
    gcg/score_maxwhite.cpp
    gcg/score_spfawh.cpp
    gcg/score_spfwh.cpp
    gcg/score_strong.cpp
    gcg/score.c
    gcg/sepa_basis.c
    gcg/sepa_original.c
    gcg/sepa_subsetrow.c
    gcg/solver_gcg.cpp
    gcg/solver_knapsack.c
    gcg/solver_mip.c
    gcg/solver.c
    gcg/stat.c
    gcg/zerohalf_selector.c
)

set(gcggraphsources
    graph/graph_gcg.cpp
    graph/graph_tclique.cpp
    graph/inst.cpp
    graph/weights.cpp
)

set(gcgheaders
    gcg/branchgcg.h
    gcg/branch_empty.h
    gcg/branch_bpstrong.h
    gcg/branch_compbnd.h
    gcg/branch_generic.h
    gcg/branch_orig.h
    gcg/branch_relpsprob.h
    gcg/branch_ryanfoster.h
    gcg/class_conspartition.h
    gcg/class_indexpartition.h
    gcg/gcgvarhistory.h
    gcg/miscvisualization.h
    gcg/class_pricingcontroller.h
    gcg/class_pricingtype.h
    gcg/class_partialdecomp.h
    gcg/class_detprobdata.h
    gcg/class_stabilization.h
    gcg/class_varpartition.h
    gcg/clscons.h
    gcg/clsvar.h
    gcg/colpool.h
    gcg/cons_decomp.h
    gcg/cons_decomp.hpp
    gcg/cons_integralorig.h
    gcg/cons_masterbranch.h
    gcg/cons_origbranch.h
    gcg/dec_compgreedily.h
    gcg/dec_connected_noNewLinkingVars.h
    gcg/dec_connectedbase.h
    gcg/dec_consclass.h
    gcg/dec_constype.h
    gcg/dec_dbscan.h
    gcg/dec_densemasterconss.h
    gcg/dec_generalmastersetcover.h
    gcg/dec_generalmastersetpack.h
    gcg/dec_generalmastersetpart.h
    gcg/dec_mastersetcover.h
    gcg/dec_mastersetpack.h
    gcg/dec_mastersetpart.h
    gcg/dec_mst.h
    gcg/dec_neighborhoodmaster.h
    gcg/dec_postprocess.h
    gcg/dec_staircase_lsp.h
    gcg/dec_stairheur.h
    gcg/dec_varclass.h
    gcg/decomp.h
    gcg/def.h
    gcg/dialog_gcg.h
    gcg/dialog_graph.h
    gcg/dialog_master.h
    gcg/disp_gcg.h
    gcg/disp_master.h
    gcg/event_bestsol.h
    gcg/event_display.h
    gcg/event_mastersol.h
    gcg/event_relaxsol.h
    gcg/event_mastersepacut.h
    gcg/event_solvingstats.h
    gcg/gcg.h
    gcg/gcgcol.h
    gcg/gcg_general.h
    gcg/gcggithash.h
    gcg/gcgplugins.h
    gcg/gcgpqueue.h
    gcg/gcgsort.h
    gcg/heur_gcgcoefdiving.h
    gcg/heur_gcgdins.h
    gcg/heur_gcgfeaspump.h
    gcg/heur_gcgfracdiving.h
    gcg/heur_gcgguideddiving.h
    gcg/heur_gcglinesdiving.h
    gcg/heur_gcgpscostdiving.h
    gcg/heur_gcgrens.h
    gcg/heur_gcgrins.h
    gcg/heur_gcgrounding.h
    gcg/heur_gcgshifting.h
    gcg/heur_gcgsimplerounding.h
    gcg/heur_gcgveclendiving.h
    gcg/heur_gcgzirounding.h
    gcg/heur_greedycolsel.h
    gcg/heur_ipcolgen.h
    gcg/heur_mastercoefdiving.h
    gcg/heur_masterdiving.h
    gcg/heur_masterfracdiving.h
    gcg/heur_masterlinesdiving.h
    gcg/heur_mastervecldiving.h
    gcg/heur_origdiving.h
    gcg/heur_relaxcolsel.h
    gcg/heur_restmaster.h
    gcg/heur_setcover.h
    gcg/heur_xpcrossover.h
    gcg/heur_xprins.h
    gcg/mastersepacut.h
    gcg/extendedmasterconsdata.h
    gcg/masterplugins.h
    gcg/nodesel_master.h
    gcg/objdialog.h
    gcg/objpricer_gcg.h
    gcg/params_visu.h
    gcg/presol_roundbound.h
    gcg/pricer_gcg.h
    gcg/pricestore_gcg.h
    gcg/pub_pricingcb.h
    gcg/pricingjob.h
    gcg/pricingprob.h
    gcg/pub_branchgcg.h
    gcg/pub_clscons.h
    gcg/pub_clsvar.h
    gcg/pub_colpool.h
    gcg/pub_decomp.h
    gcg/pub_gcg.h
    gcg/pub_gcg.hpp
    gcg/pub_gcgcol.h
    gcg/pub_gcgheur.h
    gcg/pub_gcgsepa.h
    gcg/pub_gcgpqueue.h
    gcg/pub_gcgvar.h
    gcg/pub_extendedmasterconsdata.h
    gcg/pub_mastersepacut.h
    gcg/pub_pricingcb.h
    gcg/pub_pricingjob.h
    gcg/pub_pricingprob.h
    gcg/pub_score.h
    gcg/pub_solver.h
    gcg/reader_blk.h
    gcg/reader_cls.h
    gcg/reader_dec.h
    gcg/reader_gp.h
    gcg/reader_ref.h
    gcg/reader_tex.h
    gcg/relax_gcg.h
    gcg/scip_misc.h
    gcg/score_bender.h
    gcg/score_border.h
    gcg/score_classic.h
    gcg/score_fawh.h
    gcg/score_forswh.h
    gcg/score_maxwhite.h
    gcg/score_spfawh.h
    gcg/score_spfwh.h
    gcg/score_strong.h
    gcg/score.h
    gcg/sepa_basis.h
    gcg/sepa_original.h
    gcg/sepa_subsetrow.h
    gcg/solver.h
    gcg/solver_cliquer.h
    gcg/solver_gcg.h
    gcg/solver_knapsack.h
    gcg/solver_mip.h
    gcg/solver_xyz.h
    gcg/stat.h
    gcg/struct_branchgcg.h
    gcg/struct_colpool.h
    gcg/struct_decomp.h
    gcg/struct_detector.h
    gcg/struct_gcg.h
    gcg/struct_gcgcol.h
    gcg/struct_gcgpqueue.h
    gcg/struct_mastersepacut.h
    gcg/struct_extendedmasterconsdata.h
    gcg/struct_pricestore_gcg.h
    gcg/struct_pricingcb.h
    gcg/struct_pricingjob.h
    gcg/struct_pricingprob.h
    gcg/struct_score.h
    gcg/struct_sepagcg.h
    gcg/struct_solver.h
    gcg/struct_vardata.h
    gcg/struct_gcgvarhistory.h
    gcg/type_branchgcg.h
    gcg/type_classifier.h
    gcg/type_colpool.h
    gcg/type_consclassifier.h
    gcg/type_decomp.h
    gcg/type_detector.h
    gcg/type_gcg.h
    gcg/type_gcgcol.h
    gcg/type_gcgpqueue.h
    gcg/type_mastersepacut.h
    gcg/type_gcgpricer.h
    gcg/type_locks.h
    gcg/type_extendedmasterconsdata.h
    gcg/type_masterdiving.h
    gcg/type_origdiving.h
    gcg/type_parameter.h
    gcg/type_pricestore_gcg.h
    gcg/type_pricetype.h
    gcg/type_pricingcb.h
    gcg/type_pricingjob.h
    gcg/type_pricingprob.h
    gcg/type_pricingstatus.h
    gcg/type_sepagcg.h
    gcg/type_score.h
    gcg/type_solver.h
    gcg/type_varclassifier.h
    gcg/zerohalf_selector.h
    gcg/type_gcgvarhistory.h
)

set(gcggraphheaders
    graph/bipartitegraph.h
    graph/bipartitegraph_def.h
    graph/bridge.h
    graph/columngraph.h
    graph/columngraph_def.h
    graph/graph.h
    graph/graph_def.h
    graph/graph_gcg.h
    graph/graph_interface.h
    graph/graph_tclique.h
    graph/graphalgorithms.h
    graph/graphalgorithms_def.h
    graph/hypercolgraph.h
    graph/hypercolgraph_def.h
    graph/hypergraph.h
    graph/hypergraph_def.h
    graph/hyperrowcolgraph.h
    graph/hyperrowcolgraph_def.h
    graph/hyperrowgraph.h
    graph/hyperrowgraph_def.h
    graph/matrixgraph.h
    graph/matrixgraph_def.h
    graph/rowgraph.h
    graph/rowgraph_def.h
    graph/weights.h
)

if(NOT SYM STREQUAL "none")
   set(gcgsources ${gcgsources} ${symsources}
      symmetry/automorphism.cpp
      gcg/dec_isomorph.cpp
   )
   set(gcgheaders ${gcgheaders}
      symmetry/automorphism.h
      symmetry/automorphism.hpp
      symmetry/pub_automorphism.h
      gcg/dec_isomorph.h
   )
endif()

if(SYM STREQUAL "bliss" OR SYM STREQUAL "sbliss")
   set(gcgsources ${gcgsources}
      symmetry/automorphism_bliss.cpp
   )
endif()

if(SYM STREQUAL "nauty" OR SYM STREQUAL "snauty")
   set(gcgsources ${gcgsources}
      symmetry/automorphism_nauty.cpp
   )
endif()

if(CLIQUER_FOUND)
   set(gcgsources ${gcgsources}
      gcg/solver_cliquer.c
   )
endif()

if(CPLEX_FOUND)
   set(gcgsources ${gcgsources}
      gcg/solver_cplex.c
   )
   set(gcgheaders ${gcgheaders}
      gcg/solver_cplex.h
   )
endif()

if(HIGHS_FOUND)
   set(gcgsources ${gcgsources}
      gcg/solver_highs.c
   )
   set(gcgheaders ${gcgheaders}
      gcg/solver_highs.h
   )
endif()

if(HMETIS_FOUND AND NOT WIN32)
    add_definitions(-DHMETIS_HEADER)
    configure_file(hmetis.h.in gcg/hmetis.h)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    set(gcgsources ${gcgsources}
        gcg/dec_hcgpartition.cpp
        gcg/dec_hrcgpartition.cpp
        gcg/dec_hrgpartition.cpp
    )
    set(gcgheaders ${gcgheaders}
        gcg/dec_hcgpartition.h
        gcg/dec_hrcgpartition.h
        gcg/dec_hrgpartition.h
        ${CMAKE_CURRENT_BINARY_DIR}/gcg/hmetis.h
    )
endif()

if(JANSSON_FOUND)
    set(gcgsources ${gcgsources}
        gcg/reader_jdec.cpp
    )
    set(gcgheaders ${gcgheaders}
        gcg/reader_jdec.h
    )
endif()

if(OPENMP_FOUND)
   set(gcgheaders ${gcgheaders}
      gcg/struct_locks.h
   )
endif()

# all source files should be compiled with a c++ compiler
if(CXXONLY)
    set_source_files_properties(main.c ${gcgsources} ${gcggraphsources} PROPERTIES LANGUAGE CXX)

    # for the clang compiler this suppresses the warnings about treating 'c' input as 'c++' when CXXONLY is enabled
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
       add_compile_options(-x c++)
    endif()
endif()

source_group("gcg" FILES ${gcgsources} ${gcgheaders})
source_group("graph" FILES ${gcggraphsources} ${gcggraphheaders})
source_group("main" FILES main.c)

add_library(libgcg ${gcgsources} ${gcggraphsources} ${gcgheaders} ${gcggraphheaders})

if(MSVC)
   # msvc otherwise is not smart enough and tries to link the gcg.exe binary as a library
   setLibProperties(libgcg "libgcg")
else()
   setLibProperties(libgcg "gcg")
endif()

target_link_libraries(libgcg
    PRIVATE
    ${SCIP_LIBRARIES}
    ${SYM_LIBRARIES}
    ${CLIQUER_LIBRARIES}
    ${CPLEX_LIBRARIES}
    $<$<BOOL:${OPENMP_FOUND}>:OpenMP::OpenMP_C>
    $<$<BOOL:${OPENMP_FOUND}>:OpenMP::OpenMP_CXX>
    $<$<BOOL:${GSL_FOUND}>:GSL::gsl>
    $<$<BOOL:${JANSSON_FOUND}>:jansson::jansson>
    $<$<BOOL:${HIGHS_FOUND}>:highs::highs>)

if(SHARED)
    set(SCIP_NEEDED 0)
    target_link_libraries(libgcg PRIVATE ${SCIP_PIC_LIBRARIES} ${SYM_PIC_LIBRARIES})
    add_executable(gcg main.c ${gcgsources} ${gcggraphsources} ${gcgheaders} ${gcggraphheaders})
else()
    set(SCIP_NEEDED ${SCIP_FOUND})
    target_link_libraries(libgcg PRIVATE ${SCIP_LIBRARIES} ${SYM_LIBRARIES})
    add_executable(gcg main.c)
    target_link_libraries(gcg libgcg)
endif()

add_executable(GCG::GCG ALIAS gcg)
add_library(GCG::libgcg ALIAS libgcg)

target_compile_options(gcg PRIVATE ${GCG_COMP_OPTIONS})
target_compile_options(libgcg PRIVATE ${GCG_COMP_OPTIONS})

include(GenerateExportHeader)
generate_export_header(libgcg BASE_NAME gcg EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/gcg/gcg_export.h)
target_compile_definitions(gcg PRIVATE GCG_STATIC_DEFINE)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if (SANITIZE_ADDRESS OR SANITIZE_MEMORY OR SANITIZE_THREAD OR SANITIZE_UNDEFINED)
       find_package(Sanitizers)
       add_sanitizers(gcg)
    endif()
endif()

target_link_libraries(gcg
    ${SCIP_LIBRARIES}
    ${SYM_LIBRARIES}
    ${CLIQUER_LIBRARIES}
    ${CPLEX_LIBRARIES}
    $<$<BOOL:${OPENMP_FOUND}>:OpenMP::OpenMP_C>
    $<$<BOOL:${OPENMP_FOUND}>:OpenMP::OpenMP_CXX>
    $<$<BOOL:${GSL_FOUND}>:GSL::gsl>
    $<$<BOOL:${JANSSON_FOUND}>:jansson::jansson>
    $<$<BOOL:${HIGHS_FOUND}>:highs::highs>)

add_dependencies(libgcg gcg_update_githash)
add_dependencies(gcg gcg_update_githash)

set_target_properties(libgcg PROPERTIES
    VERSION ${GCG_VERSION_MAJOR}.${GCG_VERSION_MINOR}.${GCG_VERSION_PATCH}.${GCG_VERSION_SUB}
    SOVERSION ${GCG_VERSION_MAJOR}.${GCG_VERSION_MINOR}
    INSTALL_RPATH_USE_LINK_PATH TRUE)

# set the install rpath to the installed destination
set_target_properties(gcg PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
    INSTALL_RPATH_USE_LINK_PATH TRUE)

# install the header files of gcg
install(FILES ${gcgheaders} DESTINATION include/gcg)
# separate installation of graph headers to preserve the folder structure
# see e.g. here: https://stackoverflow.com/a/11097012
install(FILES ${gcggraphheaders} DESTINATION include/graph)

# install the binary and the library to appropriate locations and add them to an export group
install(TARGETS gcg libgcg EXPORT gcg-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Add all targets to the build-tree export set
export(TARGETS gcg libgcg
  FILE "${CMAKE_BINARY_DIR}/gcg-targets.cmake")

#make scip dir absolute for the config file
if(SCIP_FOUND)
    get_filename_component(CONF_SCIP_DIR ${SCIP_DIR} REALPATH BASE_DIR ${CMAKE_SOURCE_DIR})
endif()

#configure the config file for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src")
configure_file(${PROJECT_SOURCE_DIR}/gcg-config.cmake.in
  "${CMAKE_BINARY_DIR}/gcg-config.cmake" @ONLY)

#configure the config file for the install
set(CONF_INCLUDE_DIRS "\${CMAKE_CURRENT_LIST_DIR}/../../../include")
configure_file(${PROJECT_SOURCE_DIR}/gcg-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gcg-config.cmake" @ONLY)

# install the targets of the gcg export group and the config file so that other projects
# can link easily against gcg
install(EXPORT gcg-targets FILE gcg-targets.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gcg)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gcg-config.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gcg)
