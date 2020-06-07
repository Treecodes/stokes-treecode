cmake_minimum_required (VERSION 2.8.0)

if (USER_PARAM_n)
  message (STATUS "USER_PARAM_n = " ${USER_PARAM_n})
else ()
  set(USER_PARAM_n 3)
  message (STATUS "USER_PARAM_n = " ${USER_PARAM_n} " (default)")
endif ()

if (USER_PARAM_N)
  message (STATUS "USER_PARAM_N = " ${USER_PARAM_N})
else ()
  set(USER_PARAM_N 10000)
  message (STATUS "USER_PARAM_N = " ${USER_PARAM_N} " (default)")
endif ()

if (USER_PARAM_N0)
  message (STATUS "USER_PARAM_N0 = " ${USER_PARAM_N0})
else ()
  set(USER_PARAM_N0 2000)
  message (STATUS "USER_PARAM_N0 = " ${USER_PARAM_N0} " (default)")
endif ()

if (USER_PARAM_theta)
  message (STATUS "USER_PARAM_theta = " ${USER_PARAM_theta})
else ()
  set(USER_PARAM_theta 0.7)
  message (STATUS "USER_PARAM_theta = " ${USER_PARAM_theta} " (default)")
endif ()

if (USER_PARAM_eps)
  message (STATUS "USER_PARAM_eps = " ${USER_PARAM_eps})
else ()
  set(USER_PARAM_eps 0.02)
  message (STATUS "USER_PARAM_eps = " ${USER_PARAM_eps} " (default)")
endif ()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"  )

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")

# add the executable
add_executable(ctreecode KITC_MRS.cpp)
