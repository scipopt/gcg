include(CTest)

add_dependencies(check
   gcg_check)

if(MSVC)
    add_custom_target(gcg_check
        COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIG> -R "-default" --output-on-failure --schedule-random
        DEPENDS gcg)
else()
    add_custom_target(gcg_check
        COMMAND ${CMAKE_CTEST_COMMAND} -R "-default" --output-on-failure --schedule-random
        DEPENDS gcg)
endif()



file(RELATIVE_PATH GCG_REL_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gcg)
add_custom_target(gcg_cluster
   COMMAND ${CMAKE_COMMAND} -E env "GCG_BINARY=${GCG_REL_PATH}" ${PROJECT_SOURCE_DIR}/check/cmake_wrapper.sh cluster
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/check/
   DEPENDS gcg)

add_custom_target(gcg_visu
   COMMAND ${CMAKE_COMMAND} -E env "GCG_BINARY=${GCG_REL_PATH}" ${PROJECT_SOURCE_DIR}/check/cmake_wrapper.sh visu
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/check/
   DEPENDS gcg)

add_custom_target(gcg_test
   COMMAND ${CMAKE_COMMAND} -E env "GCG_BINARY=${GCG_REL_PATH}" ${PROJECT_SOURCE_DIR}/check/cmake_wrapper.sh test
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/check/
   DEPENDS gcg)

#
# define the instance sets and settings.
#
# semicolon '\;' is used to split an instance and its optimal objective value
# For infeasible instances, '+infinity' is used (or '-infinity' in case of maximization)

# In order to run a sub-suite, use -R. `make check` just runs anything matching default.*.
#

#
# MIPLIB instances
#

set(settings_MIPLIB
    "default"
)

set(instances_MIPLIB
   "instances/miplib/noswot.mps\;instances/miplib/noswot.dec\;-41"
   )

#
# BPP instances
#

set(settings_BPP
    "default"
)

set(instances_BPP
   "instances/bpp/N1C1W4_M.BPP.lp\;instances/bpp/N1C1W4_M.BPP.dec\;41"
   "instances/bpp/N1C2W2_O.BPP.lp\;instances/bpp/N1C2W2_O.BPP.dec\;29"
   "instances/bpp/N1C3W1_A.lp\;instances/bpp/N1C3W1_A.blk\;16"
   )

#
# CPMP instances
#

set(settings_CPMP
    "default"
)

set(instances_CPMP
   "instances/cpmp/p1250-2.lp\;instances/cpmp/p1250-2.blk\;412"
   "instances/cpmp/p1650-2.txt.lp\;\;336"
   "instances/cpmp/p2050-1.txt.lp\;instances/cpmp/p2050-1.txt.dec\;266"
   )

#
# CS instances
#

set(settings_CS
    "default"
)

set(instances_CS
   "instances/cs/TEST0055.lp\;instances/cs/TEST0055.dec\;11"
   "instances/cs/TEST0059.lp\;instances/cs/TEST0059.dec\;11"
   )

#
# GAP instances
#

set(settings_GAP
    "default"
)

set(instances_GAP
   "instances/gap/gap4_2.txt.lp\;\;644"
   "instances/gap/gap8_4.txt.lp\;instances/gap/gap8_4.txt.dec\;1117"
   )

#
# MKP instances
#

set(settings_MKP
    "default"
)

set(instances_MKP
   "instances/mkp/strong_d_45_15_12.lp\;\;-8703"
   "instances/mkp/strong_s_75_15_14.lp\;\;-7158"
   "instances/mkp/strong_s_75_15_18.lp\;\;-8360"
   )

#
# add a test to build the GCG binary that all further tests depend on
#
add_test(NAME gcg-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target gcg
        )

#
# avoid that several build jobs try to concurrently build the GCG library
# note that this ressource lock name is not the actual libgcg target
#
set_tests_properties(gcg-build
                    PROPERTIES
                        RESOURCE_LOCK libgcg
                    )


#
# macro to split an instance into its relevant information
# - path
# - decomp
# - optval
# - basename
#
macro(split_instance instance)
    list(GET instance 0 path)
    list(GET instance 1 decpath)
    list(GET instance 2 optval)
    get_filename_component(basename ${path} NAME)
endmacro(split_instance)

macro(add_instancetests instances settings prefix)
#
# loop over the instances
#
    foreach(instance ${${instances}})
        split_instance(instance)
        #
        # loop over all settings
        #
        foreach(setting ${${settings}})
            #
            # treat the instance as a list of three values
            #
            set(GCG_TEST_COMMAND
                $<TARGET_FILE:gcg>
                -f ${PROJECT_SOURCE_DIR}/check/${path}
                -s ${PROJECT_SOURCE_DIR}/settings/${setting}.set
                -o ${optval} ${optval}
            )
            if(NOT "${decpath}" STREQUAL "")
                list(APPEND GCG_TEST_COMMAND -d ${PROJECT_SOURCE_DIR}/check/${decpath})
            endif()

            add_test(NAME GCG-${prefix}-${setting}-${basename}
                    COMMAND ${GCG_TEST_COMMAND}
                    )
            set_tests_properties(
                GCG-${prefix}-${setting}-${basename}
                PROPERTIES
                    PASS_REGULAR_EXPRESSION "Validation         : Success"
                    DEPENDS gcg-build
                    RESOURCE_LOCK libscip
            )
            if(MSVC)
                set_tests_properties(
                    GCG-${prefix}-${setting}-${basename}
                    PROPERTIES
                        ENVIRONMENT_MODIFICATION "PATH=path_list_append:${SCIP_BINARY_DIR}/bin/$<CONFIG>/"
                )
            endif()
        endforeach(setting)
    endforeach(instance)
endmacro(add_instancetests)

add_instancetests(instances_MIPLIB settings_MIPLIB "MIPLIB")
add_instancetests(instances_BPP settings_BPP "BPP")
add_instancetests(instances_CPMP settings_CPMP "CPMP")
add_instancetests(instances_CS settings_CS "CS")
add_instancetests(instances_GAP settings_GAP "GAP")
add_instancetests(instances_MKP settings_MKP "MKP")

#this test cannot work, since gcg's main.c does not handle SIGTERMS
#
## add a test for handling of the SIGTERM signal. The test uses the timeout command that
## is only available on Linux, that is available on MAC OS as "gtimeout" after installing
## the 'coreutils' package
##
#if (UNIX)
#    if (APPLE)
#        set(timeoutcommand "gtimeout")
#    else ()
#        set(timeoutcommand "timeout")
#    endif (APPLE)
#    add_test(NAME signal-handling-sigterm-blend2
#            COMMAND ${timeoutcommand} -sSIGTERM 1 $<TARGET_FILE:gcg> -f ${CMAKE_CURRENT_SOURCE_DIR}/instances/miplib/noswot.mps
#            )
#    set_tests_properties(signal-handling-sigterm-blend2
#                        PROPERTIES
#                            DEPENDS gcg-build
#                            PASS_REGULAR_EXPRESSION "termination signal received"
#                            #
#                            # I assume that this test takes longer than 2 seconds
#                            #
#                            )
#endif (UNIX)

