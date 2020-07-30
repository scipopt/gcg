# ${YAML_INCLUDE_DIRS} contains the paths to yaml.h if YAML is found.
# ${YAML_LIBRARIES} contains libyaml if YAML is found.

# Check whether environment variable YAML_DIR was set.
if(NOT YAML_DIR)
  set(ENV_YAML_DIR $ENV{YAML_DIR})
  if(ENV_YAML_DIR)
    set(YAML_DIR $ENV{YAML_DIR} CACHE PATH "Path to yaml directory")
  endif()
endif()

find_path(YAML_INCLUDE_DIRS
    NAMES yaml.h
    HINTS ${YAML_DIR}
    PATH_SUFFIXES include)

if(STATIC_YAML)
    message(STATUS "Attempting to find libyaml.a")
    find_library(YAML_LIBRARY
        NAMES libyaml.a yaml
        HINTS ${YAML_DIR}
        PATH_SUFFIXES lib)
else()
    find_library(YAML_LIBRARY
        NAMES yaml
        HINTS ${YAML_DIR}
        PATH_SUFFIXES lib)
endif()

SET(YAML_LIBRARIES ${YAML_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAML DEFAULT_MSG YAML_INCLUDE_DIRS YAML_LIBRARIES)
