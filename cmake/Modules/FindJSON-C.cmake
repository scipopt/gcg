# ${JSON-C_INCLUDE_DIRS} contains the paths to json.h if json-c is found.
# ${JSON-C_LIBRARIES} contains libjson-c if json-c is found.

# Check whether environment variable JSON-C_DIR was set.
if(NOT JSON-C_DIR)
  set(ENV_JSON-C_DIR $ENV{JSON-C_DIR})
  if(ENV_JSON-C_DIR)
    set(JSON-C_DIR $ENV{JSON-C_DIR} CACHE PATH "Path to json-c directory")
  endif()
endif()

find_path(JSON-C_INCLUDE_DIRS
    NAMES json-c/json.h
    HINTS ${JSON-C_DIR}
    PATH_SUFFIXES include)

if(STATIC_JSON-C)
    message(STATUS "Attempting to find libjson-c.a")
    find_library(JSON-C_LIBRARY
        NAMES libjson-c.a json-c
        HINTS ${JSON-C_DIR}
        PATH_SUFFIXES lib)
else()
    find_library(JSON-C_LIBRARY
        NAMES json-c
        HINTS ${JSON-C_DIR}
        PATH_SUFFIXES lib)
endif()

SET(JSON-C_LIBRARIES ${JSON-C_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JSON-C DEFAULT_MSG JSON-C_INCLUDE_DIRS JSON-C_LIBRARIES)
