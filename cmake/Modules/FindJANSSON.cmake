# ${JANSSON_INCLUDE_DIRS} contains the paths to jansson.h if Jansson is found.
# ${JANSSON_LIBRARIES} contains libjansson if jansson is found.

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
endif()

SET(JANSSON_LIBRARIES ${JANSSON_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JANSSON DEFAULT_MSG JANSSON_INCLUDE_DIRS JANSSON_LIBRARIES)
