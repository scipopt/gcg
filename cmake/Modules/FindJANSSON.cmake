# ${JANSSON_INCLUDE_DIRS} contains the paths to jansson.h if Jansson is found.
# ${JANSSON_LIBRARIES} contains libjansson if jansson is found.
# creates target jansson::jansson

# Check whether environment variable JANSSON_DIR was set.
if(NOT JANSSON_DIR)
    set(ENV_JANSSON_DIR $ENV{JANSSON_DIR})
    if(ENV_JANSSON_DIR)
        set(JANSSON_DIR $ENV{JANSSON_DIR} CACHE PATH "Path to Jansson directory")
    endif()
endif()

find_path(JANSSON_INCLUDE_DIRS
    NAMES jansson.h
    HINTS ${JANSSON_DIR}
    PATH_SUFFIXES include)

if(STATIC_JANSSON)
    message(STATUS "Attempting to find libjansson.a")
    find_library(JANSSON_LIBRARY
        NAMES libjansson.a jansson
        HINTS ${JANSSON_DIR}
        PATH_SUFFIXES lib)
else()
    find_library(JANSSON_LIBRARY
        NAMES jansson
        HINTS ${JANSSON_DIR}
        PATH_SUFFIXES lib)
    find_library(JANSSON_LIBRARY_DEBUG
        NAMES jansson_d
        HINTS ${JANSSON_DIR}
        PATH_SUFFIXES debug/lib lib)
endif()

set(JANSSON_LIBRARIES ${JANSSON_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JANSSON DEFAULT_MSG JANSSON_INCLUDE_DIRS JANSSON_LIBRARIES)

if(WIN32)
    find_file(JANSSON_LIBRARY_DLL
        NAMES jansson.dll
        HINTS ${JANSSON_DIR}
        PATH_SUFFIXES bin lib)
    find_file(JANSSON_LIBRARY_DLL_DEBUG
        NAMES jansson_d.dll janssond.dll
        HINTS ${JANSSON_DIR}
        PATH_SUFFIXES debug/bin debug debug/lib)
endif()

if(JANSSON_FOUND AND NOT TARGET jansson::jansson)
    if(EXISTS "${JANSSON_LIBRARY_DLL}")
        add_library(jansson::jansson SHARED IMPORTED)
        set_target_properties(jansson::jansson PROPERTIES
            IMPORTED_LOCATION_RELEASE "${JANSSON_LIBRARY_DLL}"
            IMPORTED_IMPLIB "${JANSSON_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${JANSSON_INCLUDE_DIRS}"
            IMPORTED_CONFIGURATIONS Release
            IMPORTED_LINK_INTERFACE_LANGUAGES "C")
        if(EXISTS "${JANSSON_LIBRARY_DLL_DEBUG}")
            set_property(TARGET jansson::jansson APPEND PROPERTY IMPORTED_CONFIGURATIONS Debug )
            set_target_properties(jansson::jansson PROPERTIES
                IMPORTED_LOCATION_DEBUG "${JANSSON_LIBRARY_DLL_DEBUG}"
                IMPORTED_IMPLIB_DEBUG "${JANSSON_LIBRARY_DEBUG}" )
        endif()
    else()
        add_library(jansson::jansson UNKNOWN IMPORTED)
        set_target_properties(jansson::jansson PROPERTIES
            IMPORTED_LOCATION "${JANSSON_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${JANSSON_INCLUDE_DIR}")
    endif()
endif()
