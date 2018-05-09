include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#
# interface function for setting common library properties
#
function(setLibProperties targetname outputname)
    set_target_properties(${targetname} PROPERTIES
        OUTPUT_NAME ${outputname}
        MACOSX_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endfunction(setLibProperties)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED on)

set(gcgsources
    graph/graph_tclique.cpp
    graph/inst.cpp
    graph/weights.cpp
    branch_empty.c
    branch_generic.c
    branch_orig.c
    branch_relpsprob.c
    branch_ryanfoster.c
    branch_xyz.c
    colpool.c
    cons_decomp.c
    cons_integralorig.c
    cons_masterbranch.c
    cons_origbranch.c
    dec_connected.c
    dec_cutpacking.c
    dec_random.c
    dec_staircase.c
    dec_xyz.c
    decomp.c
    dialog_gcg.c
    dialog_master.c
    disp_gcg.c
    disp_master.c
    event_bestsol.c
    event_display.c
    event_mastersol.c
    event_relaxsol.c
    event_solvingstats.c
    gcgcol.c
    gcggithash.c
    gcgheur.c
    gcgplugins.c
    gcgpqueue.c
    gcgvar.c
    githash.c
    heur_gcgcoefdiving.c
    heur_gcgdins.c
    heur_gcgfeaspump.c
    heur_gcgfracdiving.c
    heur_gcgguideddiving.c
    heur_gcglinesdiving.c
    heur_gcgpscostdiving.c
    heur_gcgrens.c
    heur_gcgrins.c
    heur_gcgrounding.c
    heur_gcgshifting.c
    heur_gcgsimplerounding.c
    heur_gcgveclendiving.c
    heur_gcgzirounding.c
    heur_greedycolsel.c
    heur_mastercoefdiving.c
    heur_masterdiving.c
    heur_masterfracdiving.c
    heur_masterlinesdiving.c
    heur_mastervecldiving.c
    heur_origdiving.c
    heur_relaxcolsel.c
    heur_restmaster.c
    heur_setcover.c
    heur_xpcrossover.c
    heur_xprins.c
    masterplugins.c
    misc.c
    nodesel_master.c
    pricestore_gcg.c
    reader_blk.c
    reader_dec.c
    reader_gp.c
    reader_ref.c
    relax_gcg.c
    scip_misc.c
    sepa_basis.c
    sepa_master.c
    solver_cplex.c
    solver_knapsack.c
    solver_mip.c
    solver_xyz.c
    stat.c
    class_pricingtype.cpp
    class_stabilization.cpp
    dec_arrowheur.cpp
    dec_colors.cpp
    dec_consname.cpp
    dec_stairheur.cpp
    dialog_graph.cpp
    objdialog.cpp
    pricer_gcg.cpp
)

set(gcgheaders
    graph/bipartitegraph.h
    graph/bipartitegraph_def.h
    graph/bridge.h
    graph/columngraph.h
    graph/columngraph_def.h
    graph/graph.h
    graph/graph_def.h
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
    branch_empty.h
    branch_generic.h
    branch_orig.h
    branch_relpsprob.h
    branch_ryanfoster.h
    branch_xyz.h
    class_pricingtype.h
    class_stabilization.h
    colpool.h
    cons_decomp.h
    cons_integralorig.h
    cons_masterbranch.h
    cons_origbranch.h
    dec_arrowheur.h
    dec_colors.h
    dec_connected.h
    dec_consname.h
    dec_cutpacking.h
    dec_isomorph.h
    dec_random.h
    dec_staircase.h
    dec_stairheur.h
    dec_xyz.h
    decomp.h
    dialog_gcg.h
    dialog_graph.h
    dialog_master.h
    disp_gcg.h
    disp_master.h
    event_bestsol.h
    event_display.h
    event_mastersol.h
    event_relaxsol.h
    event_solvingstats.h
    gcg.h
    gcgcol.h
    gcggithash.h
    gcgplugins.h
    gcgpqueue.h
    heur_gcgcoefdiving.h
    heur_gcgdins.h
    heur_gcgfeaspump.h
    heur_gcgfracdiving.h
    heur_gcgguideddiving.h
    heur_gcglinesdiving.h
    heur_gcgpscostdiving.h
    heur_gcgrens.h
    heur_gcgrins.h
    heur_gcgrounding.h
    heur_gcgshifting.h
    heur_gcgsimplerounding.h
    heur_gcgveclendiving.h
    heur_gcgzirounding.h
    heur_greedycolsel.h
    heur_mastercoefdiving.h
    heur_masterdiving.h
    heur_masterfracdiving.h
    heur_masterlinesdiving.h
    heur_mastervecldiving.h
    heur_origdiving.h
    heur_relaxcolsel.h
    heur_restmaster.h
    heur_setcover.h
    heur_xpcrossover.h
    heur_xprins.h
    hmetis.h
    masterplugins.h
    nodesel_master.h
    objdialog.h
    objpricer_gcg.h
    pricestore_gcg.h
    pricer_gcg.h
    pub_bliss.h
    pub_colpool.h
    pub_decomp.h
    pub_gcgcol.h
    pub_gcgheur.h
    pub_gcgpqueue.h
    pub_gcgvar.h
    reader_blk.h
    reader_dec.h
    reader_gp.h
    reader_ref.h
    relax_gcg.h
    scip_misc.h
    sepa_basis.h
    sepa_master.h
    solver_cplex.h
    solver_knapsack.h
    solver_mip.h
    solver_xyz.h
    stat.h
    struct_branchgcg.h
    struct_decomp.h
    struct_detector.h
    struct_gcgcol.h
    struct_gcgpqueue.h
    struct_solver.h
    struct_vardata.h
    type_branchgcg.h
    type_decomp.h
    type_detector.h
    type_gcgcol.h
    type_gcgpqueue.h
    type_masterdiving.h
    type_origdiving.h
    type_solver.h
)

set(symsources
    ${sym}
)

set(symheaders
    bliss_automorph.h
    dec_isomorph.h
)

set(cliquersources
    ${cliquer}
)

if(HMETIS_FOUND)
   include_directories(${CMAKE_CURRENT_BINARY_DIR})
   add_definitions(-DHMETIS_HEADER)
   configure_file(hmetis.h.in hmetis.h)
endif()

# all source files should be compiled with a c++ compiler
if(CXXONLY)
    set_source_files_properties(main.c ${gcgsources} ${symsources} ${cliquersources} PROPERTIES LANGUAGE CXX)

    # for the clang compiler this suppresses the warnings about treating 'c' input as 'c++' when CXXONLY is enabled
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
       add_compile_options(-x c++)
    endif()
endif()

add_library(libgcg ${gcgsources} ${symsources} ${cliquersources})
setLibProperties(libgcg "gcg")

target_link_libraries(libgcg PRIVATE ${SCIP_LIBRARIES} ${SYM_LIBRARIES})

if(SHARED)
    target_link_libraries(libgcg PRIVATE ${SCIP_PIC_LIBRARIES} ${SYM_PIC_LIBRARIES} ${CLIQUER_PIC_LIBRARIES})
    add_executable(gcg main.c ${gcgsources} ${symsources} ${cliquersources})
else()
    target_link_libraries(libgcg PRIVATE ${SCIP_LIBRARIES} ${SYM_LIBRARIES} ${CLIQUER_LIBRARIES})
    add_executable(gcg main.c)
    target_link_libraries(gcg libgcg)
endif()

if(CMAKE_BUILD_TYPE EQUAL "Debug")
    find_package(Sanitizers)
    add_sanitizers(gcg)
endif()

target_compile_definitions(gcg PRIVATE EXTERN=extern)

target_link_libraries(gcg ${SCIP_LIBRARIES} ${SYM_LIBRARIES} ${CLIQUER_LIBRARIES})

add_dependencies(libgcg gcg_update_githash)
add_dependencies(gcg gcg_update_githash)

set_target_properties(libgcg PROPERTIES
    VERSION ${GCG_VERSION_MAJOR}.${GCG_VERSION_MINOR}.${GCG_VERSION_PATCH}.${GCG_VERSION_SUB}
    SOVERSION ${GCG_VERSION_MAJOR}.${GCG_VERSION_MINOR}
    INSTALL_RPATH_USE_LINK_PATH TRUE)

# set the install rpath to the installed destination
set_target_properties(gcg PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib"
    INSTALL_RPATH_USE_LINK_PATH TRUE)

# install the header files of gcg
install(FILES ${gcgheaders} DESTINATION include/gcg)
install(FILES ${symheaders} DESTINATION include/gcg)

# install the binary and the library to appropriate locations and add them to an export group
install(TARGETS gcg libgcg EXPORT gcg-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include)

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
install(EXPORT gcg-targets FILE gcg-targets.cmake DESTINATION lib/cmake/gcg)
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/gcg-config.cmake" DESTINATION lib/cmake/gcg)

