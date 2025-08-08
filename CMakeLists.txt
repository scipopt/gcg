cmake_minimum_required(VERSION 3.11)

# INTERPROCEDURAL_OPTIMIZATION is enforced when enabled.
if(POLICY CMP0069)
    cmake_policy(SET CMP0069 NEW)
endif()

# option() honors normal variables.
if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

# MSVC runtime library flags are selected by an abstraction.
if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()

# Avoid warnings about removed FindBoost.cmake until our minimum required version is >= 3.30.
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_DEBUG} ${CMAKE_C_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS_RELEASE}")

set(GCG_VERSION_MAJOR 4)
set(GCG_VERSION_MINOR 0)
set(GCG_VERSION_PATCH 0)
set(GCG_VERSION_SUB 0)
set(GCG_VERSION_API 21)

project(GCG
    VERSION ${GCG_VERSION_MAJOR}.${GCG_VERSION_MINOR}.${GCG_VERSION_PATCH}
    LANGUAGES CXX C
    DESCRIPTION "GCG is a generic decomposition solver for mixed-integer linear programs (MILPs)."
    )


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

if(SCIPOptSuite_BINARY_DIR)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${SCIPOptSuite_BINARY_DIR}/bin)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${SCIPOptSuite_BINARY_DIR}/lib)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${SCIPOptSuite_BINARY_DIR}/lib)
endif()

# path to e.g. FindGMP module
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

# make 'Release' the default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(MSVC)
    # GCG requires OpenMP features which are not supported by the native MSVC
    set(DEFAULT_OPENMP OFF)
else()
    set(DEFAULT_OPENMP ON)
endif()

option(OPENMP "should gcg be compiled with openmp" ${DEFAULT_OPENMP})
option(STATIC_GMP "prefer static gmp library" OFF)
option(GSL "should gcg be compiled with GSL" ON)
option(GMP "should gmp be linked" ON)
option(CLIQUER "should cliquer be linked" ON)
option(CPLEX "should Cplex be used for solving pricing MIPs" OFF)
option(HIGHS "should HiGHS be used for solving pricing MIPs" OFF)
option(MT "use static runtime libraries for Visual Studio compiler" OFF)
option(GCG_DEV_BUILD "compile SCIP and SoPlex located in lib/scip-git and lib/soplex-git" OFF)
option(HMETIS "should hmetis be called" OFF)
option(JSON "should JSON library be linked" ON)
option(LTO "Enable Link Time Optimization, aka Interprocedural Optimization" ON)

if(GCG_DEV_BUILD)
    # @todo: fix this
    set(LTO OFF)
endif()

if(LTO AND (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo"))
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result OUTPUT output LANGUAGES C CXX)
    if(result)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
        message(STATUS "LTO enabled")
    else()
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
        message(STATUS "${output}")
    endif()
endif()

# OPENMP
if(OPENMP)
   find_package(OpenMP)
   if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(THREADSAFE ON FORCE)
   endif()
endif()

# GSL
if(GSL)
   find_package(GSL 2.0)
   if(GSL_FOUND)
      include_directories(${GSL_INCLUDE_DIRS})
      add_definitions(-DWITH_GSL)
      set(GSL_PIC_LIBRARIES ${GSL_LIBRARIES})
   else()
      set(GSL_LIBRARIES)
      set(GSL_PIC_LIBRARIES)
   endif()
endif()

# GMP
if(GMP)
   find_package(GMP)
endif()
if(GMP_FOUND)
   include_directories(${GMP_INCLUDE_DIRS})
   add_definitions(-DWITH_GMP)
else()
   set(GMP_LIBRARIES)
endif()

# CLIQUER
if(CLIQUER)
   find_package(CLIQUER)
endif()
if(CLIQUER_FOUND)
   include_directories(${CLIQUER_INCLUDE_DIRS})
   add_definitions(-DWITH_CLIQUER)
   set(CLIQUER_PIC_LIBRARIES ${CLIQUER_LIBRARIES})
else()
   set(CLIQUER_LIBRARIES)
   set(CLIQUER_PIC_LIBRARIES)
endif()

# CPLEX
if(CPLEX)
   find_package(CPLEX REQUIRED)
endif()
if(CPLEX_FOUND)
   add_definitions(-DWITH_CPLEXSOLVER)
   include_directories(${CPLEX_INCLUDE_DIRS})
   set(CPLEX_PIC_LIBRARIES ${CPLEX_LIBRARIES})
else()
   set(CPLEX_LIBRARIES)
   set(CPLEX_PIC_LIBRARIES)
endif()

# HiGHS
if(HIGHS)
   find_package(HIGHS REQUIRED)
endif()
if(HIGHS_FOUND)
   find_package(Threads)
   add_definitions(-DWITH_HIGHS)
   include_directories(${HIGHS_INCLUDE_DIRS})
   set(HIGHS_LIBRARIES highs::highshighs) # I would expect that this should be set by highs-config.cmake
   set(HIGHS_PIC_LIBRARIES ${HIGHS_LIBRARIES})
else()
   set(HIGHS_LIBRARIES)
   set(HIGHS_PIC_LIBRARIES)
endif()

# HMETIS
if(HMETIS)
   find_package(HMETIS)
endif()
if(HMETIS_FOUND)
   add_definitions(-DWITH_HMETIS)
endif()

# JSON
if(JSON)
   find_package(JANSSON)
endif()
if(JANSSON_FOUND)
   include_directories(${JANSSON_INCLUDE_DIRS})
   add_definitions(-DWITH_JSON)
   set(JANSSON_PIC_LIBRARIES ${JANSSON_LIBRARIES})
else()
   set(JANSSON_LIBRARIES)
   set(JANSSON_PIC_LIBRARIES)
endif()

set(SYM snauty CACHE STRING "options for symmetry computation")  #create the variable
set_property(CACHE SYM PROPERTY STRINGS sbliss bliss nauty snauty none )  #define list of values GUI will offer for the variable

#set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

#set defines for Windows
if(WIN32)
    add_definitions(-DNO_SIGACTION)
    add_definitions(-DNO_STRTOK_R)
endif()
if(MSVC)
    add_definitions(/wd4100)
    add_definitions(/wd4244)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    add_compile_options(/bigobj)
endif()

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

# create a target for updating the current git hash
file(WRITE ${PROJECT_BINARY_DIR}/gcg_update_githash.cmake "
find_program(GIT git)
if(EXISTS \${DST})
   file(STRINGS \${DST} GITHASH_OLD)
   string(REGEX REPLACE \"#define GCG_GITHASH \\\"(.*)\\\"\" \"\\\\1\" GITHASH_OLD \"\${GITHASH_OLD}\")
endif()
if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REGEX REPLACE \"^.*-g\" \"\" GITHASH \"\${GITHASH}\")
   if(NOT \${GITHASH} STREQUAL \"\${GITHASH_OLD}\")
      file(WRITE \${DST} \"#define GCG_GITHASH \\\"\${GITHASH}\\\"\n\")
   endif()
else()
   set(GITHASH \${GITHASH_OLD})
   if(NOT GITHASH)
      message(STATUS \"Compiling without git information\")
      set(GITHASH \"NoGitInfo\")
   endif()
   file(WRITE \${DST} \"#define GCG_GITHASH \\\"\${GITHASH}\\\"\n\")
endif()
message(STATUS \"Git hash: \" \${GITHASH})
")
add_custom_target(gcg_update_githash
   COMMAND ${CMAKE_COMMAND}
      -DDST=${PROJECT_SOURCE_DIR}/src/gcg/githash.c
      -P ${PROJECT_BINARY_DIR}/gcg_update_githash.cmake)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)

# set function visibility default to hidden
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# SCIP
if(NOT GCG_DEV_BUILD)
   find_package(SCIP CONFIG REQUIRED)
   if(SCIP_FOUND)
      # SCIP headers can be directly included
      include_directories(${SCIP_INCLUDE_DIRS})

      if(ZLIB AND NOT TARGET ZLIB::ZLIB)
         find_package(ZLIB REQUIRED)
      endif()
      if(READLINE)
         find_package(Readline REQUIRED)
      endif()
   endif()
endif()

if(NOT DEFINED GCG_DEV_BUILD_SOPLEX_PATH)
   set(GCG_DEV_BUILD_SOPLEX_PATH ${PROJECT_SOURCE_DIR}/lib/soplex-git)
endif()

if(NOT DEFINED GCG_DEV_BUILD_SCIP_PATH)
   set(GCG_DEV_BUILD_SCIP_PATH ${PROJECT_SOURCE_DIR}/lib/scip-git)
endif()

if((NOT SCIP_FOUND) AND (EXISTS ${GCG_DEV_BUILD_SCIP_PATH}) AND (EXISTS ${GCG_DEV_BUILD_SOPLEX_PATH}))
   message(STATUS "Using '${GCG_DEV_BUILD_SOPLEX_PATH}'")
   message(STATUS "Using '${GCG_DEV_BUILD_SCIP_PATH}'")
   set(SOPLEX_DIR ${CMAKE_BINARY_DIR})
   add_subdirectory(${GCG_DEV_BUILD_SOPLEX_PATH})
   add_subdirectory(${GCG_DEV_BUILD_SCIP_PATH})
   set(SCIP_LIBRARIES libscip)
elseif(GCG_DEV_BUILD)
   message(FATAL_ERROR "GCG_DEV_BUILD option is set but could not find SCIP and SoPlex in the lib/ subdirectory. Please specify GCG_DEV_BUILD_SOPLEX_PATH and GCG_DEV_BUILD_SCIP_PATH")
endif()

#search the selected symmetry computation program
if(SYM STREQUAL "bliss" OR SYM STREQUAL "sbliss")
   if(NOT TARGET Bliss::libbliss)
      find_package(Bliss CONFIG HINTS ${BLISS_DIR})
      if(NOT TARGET Bliss::libbliss)
         find_package(Bliss REQUIRED)
         include_directories(${BLISS_INCLUDE_DIRS})
      endif()
   endif()
   set(BLISS_LIBRARIES Bliss::libbliss)
   set(SYM_LIBRARIES ${BLISS_LIBRARIES} ${GMP_LIBRARIES})
   set(SYM_PIC_LIBRARIES ${BLISS_LIBRARIES} ${GMP_LIBRARIES})
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BLISS_DEFINITIONS} -DWITH_BLISS")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BLISS_DEFINITIONS} -DWITH_BLISS")
   # if sbliss: sassy is a header-only library and is already included by SCIP
elseif(SYM STREQUAL "nauty" OR SYM STREQUAL "snauty")
   if(NOT GCG_DEV_BUILD AND NOT SCIPOptSuite_SOURCE_DIR)
      find_package(NAUTY REQUIRED)
      include_directories(${NAUTY_INCLUDE_DIRS})
      set(SYM_LIBRARIES ${NAUTY_LIBRARIES})
      set(SYM_PIC_LIBRARIES ${NAUTY_LIBRARIES})
   else()
      message(STATUS "Using nauty provided by SCIP")
      if(GCG_DEV_BUILD)
         set(NAUTY_SRC_DIR ${GCG_DEV_BUILD_SCIP_PATH}/src/nauty/)
      elseif(SCIPOptSuite_SOURCE_DIR)
         set(NAUTY_SRC_DIR ${SCIPOptSuite_SOURCE_DIR}/scip/src/nauty/)
      endif()
      set(symsources
         ${NAUTY_SRC_DIR}/nauty.c
         ${NAUTY_SRC_DIR}/nautil.c
         ${NAUTY_SRC_DIR}/nausparse.c
         ${NAUTY_SRC_DIR}/schreier.c
         ${NAUTY_SRC_DIR}/naurng.c
      )
      set(SYM_LIBRARIES)
      set(SYM_PIC_LIBRARIES)
   endif()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NAUTY_DEFINITIONS} -DWITH_NAUTY")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${NAUTY_DEFINITIONS} -DWITH_NAUTY")
   # if snauty: sassy is a header-only library and is already included by SCIP
else()
   message(STATUS "GCG will be compiled without symmetry computation")
   set(SYM_LIBRARIES)
   set(SYM_PIC_LIBRARIES)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BLISS_DEFINITIONS} -DNO_AUT_LIB")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BLISS_DEFINITIONS} -DNO_AUT_LIB")
endif()

# if changing these flags, also update GCCWARN/GXXWARN in make/make.project
set(ADD_C_AND_CXX_FLAGS -pedantic -Wall -W -Wpointer-arith -Wcast-align -Wwrite-strings -Wshadow  -Wredundant-decls -Wdisabled-optimization -Wno-long-long -Wno-unknown-pragmas -Wno-unused-parameter)
set(ADD_C_FLAGS -Wsign-compare -Wstrict-prototypes -Wmissing-declarations -Wmissing-prototypes -Wno-override-init -Wno-uninitialized -Wno-maybe-uninitialized)
set(ADD_CXX_FLAGS -Wnon-virtual-dtor -Wreorder -Woverloaded-virtual -Wsign-promo -Wsynth -Wcast-qual -Wno-strict-overflow -Wno-psabi)

# disable fused floating point contraction to enhance reproducibility across compilers and architectures
if(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/fp:precise)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    add_compile_options(-fp-model=precise)
    # to prevent symbols from svml exported twice we use the -shared-intel flag, see https://git.zib.de/integer/scip/issues/2872
    # this happens with the cplex static lib
    if(LPS MATCHES "cpx")
        add_compile_options(-shared-intel) # for c and cxx
    endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # require at least gcc 4.8
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
        message(WARNING "GCC version not supported, should be at least 4.8!")
    endif()
    add_compile_options(-ffp-contract=off ${ADD_C_AND_CXX_FLAGS})
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${ADD_C_FLAGS}>$<$<COMPILE_LANGUAGE:CXX>:${ADD_CXX_FLAGS}>")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # clang specific suppressions
    set(ADD_C_AND_CXX_FLAGS ${ADD_C_AND_CXX_FLAGS} -Wno-unknown-warning-option)
    set(ADD_C_FLAGS ${ADD_C_FLAGS} -Wno-initializer-overrides)
    add_compile_options(-ffp-contract=off ${ADD_C_AND_CXX_FLAGS})
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${ADD_C_FLAGS}>$<$<COMPILE_LANGUAGE:CXX>:${ADD_CXX_FLAGS}>")
endif()

# TODO: In SCIP, license is in COPYING file and not in LICENSE.
# TODO: SCIP has src/scipbuildflag.c which adds a function that returns the used build flags. Do we also want this?
# TODO: Symmetry handling is separate from the remaining code: there is a symmetry/* code that is compiled differently depending on whether BLISS was found or not. the other tools (presol_sym*) is always included, but can query the symmetry code whether it can actually detect symmetry. This is slightly cleaner from a build system point of view.

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${SCIP_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${SCIP_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${SCIP_VERSION_PATCH}")
set(CPACK_PACKAGE_VENDOR "Operations Research, RWTH Aachen University and Zuse Institute Berlin (ZIB)")
set(CPACK_PACKAGE_CONTACT "http://www.or.rwth-aachen.de/gcg")
include(CPack)

# go to src/ and compile the code
add_subdirectory(src)

#
# we set the GCG_DIR variable explicitly for the following examples/applications and unittests that depend on GCG.
#
set(GCG_DIR ${CMAKE_BINARY_DIR})

#
# add GCG tests
#
include(CTest)
if(BUILD_TESTING)
   if(NOT TARGET check)
      add_custom_target(check)
   endif()
   add_subdirectory(check)

   #
   # add GCG documentation
   #
   if(NOT TARGET doc)
      add_custom_target(doc)
   endif()
   add_subdirectory(doc)

   #
   # add unit tests as a single target. Including tests will add the unit tests as single executables
   #
   add_custom_target(gcg_unittests)
   #add_subdirectory(tests EXCLUDE_FROM_ALL)
   add_custom_target(gcg_all_executables DEPENDS gcg gcg_unittests)

   enable_testing()
endif()

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
   include(FeatureSummary)
   feature_summary(WHAT ALL)
endif()
