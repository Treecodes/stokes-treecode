cmake_minimum_required (VERSION 2.8.0)

if (USER_PARAM_P)
  message (STATUS "USER_PARAM_P = " ${USER_PARAM_P})
else ()
  set(USER_PARAM_P 3)
  message (STATUS "USER_PARAM_P = " ${USER_PARAM_P} " (default)")
endif ()

if (USER_PARAM_NUM_PARS)
  message (STATUS "USER_PARAM_NUM_PARS = " ${USER_PARAM_NUM_PARS})
else ()
  set(USER_PARAM_NUM_PARS 10000)
  message (STATUS "USER_PARAM_NUM_PARS = " ${USER_PARAM_NUM_PARS} " (default)")
endif ()

if (USER_PARAM_N0)
  message (STATUS "USER_PARAM_N0 = " ${USER_PARAM_N0})
else ()
  set(USER_PARAM_N0 2000)
  message (STATUS "USER_PARAM_N0 = " ${USER_PARAM_N0} " (default)")
endif ()

if (USER_PARAM_SQ_THETA)
  message (STATUS "USER_PARAM_SQ_THETA = " ${USER_PARAM_SQ_THETA})
else ()
  set(USER_PARAM_SQ_THETA 0.49)
  message (STATUS "USER_PARAM_SQ_THETA = " ${USER_PARAM_SQ_THETA} " (default)")
endif ()

if (USER_PARAM_DEL)
  message (STATUS "USER_PARAM_DEL = " ${USER_PARAM_DEL})
else ()
  set(USER_PARAM_DEL 0.02)
  message (STATUS "USER_PARAM_DEL = " ${USER_PARAM_DEL} " (default)")
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
add_executable(ctreecode KITC_MRS_Combine_Data.cpp)
