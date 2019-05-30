#This file shows how to link to a PUMI
#installation using CMake
#it represents a simple 'CMakeLists.txt'
#file for a new project

cmake_minimum_required(VERSION 3.0.0)

project(MeshCreate VERSION 1.0.0 LANGUAGES CXX)
option(HAVE_SIMMETRIX "If Pumi was built with Simmetrix support" OFF)

# Starting here are the critical lines:

# Allow the user to indicate where they installed SCOREC
# via "-DSCOREC_PREFIX=/home/somewhere" when calling `cmake`
set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")

# If SCOREC_PREFIX was specified, only link to that directory,
# i.e. don't link to another installation in /usr/lib by mistake
if (SCOREC_PREFIX)
  find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
# IF SCOREC_PREFIX was not specified, look in typical system directories,
# and also in CMAKE_PREFIX_PATH (environment variable)
  find_package(
      SCOREC #package name, has to be SCOREC
      2.1.0  #version. can be omitted, and will match any installed version
             #greater than or equal to this one, as long as the major number
             #is the same
      REQUIRED #indicate that SCOREC is really needed to compile
      CONFIG   #skip the 'MODULE' search system, save some time and confusion
      )

endif()


MESSAGE( STATUS "CMAKE_BINARY_DIR:     " ${CMAKE_BINARY_DIR})
MESSAGE( STATUS "Pumi configuration file location: " ${SCOREC_CONFIG})
MESSAGE( STATUS "installation prefix: " ${CMAKE_INSTALL_PREFIX})

# Figure out of Pumi was built with Simmetrix support
get_target_property(SCOREC_LIBS SCOREC::core INTERFACE_LINK_LIBRARIES)
#MESSAGE(STATUS ${SCOREC_LIBS})
if(SCOREC_LIBS MATCHES "gmi_sim")
  MESSAGE(STATUS "Found gmi_sim")
  set (HAVE_SIMMETRIX ON)
else()
  MESSAGE(STATUS "gmi_sim not found")
endif()

configure_file (
  "${PROJECT_SOURCE_DIR}/pumiInterface_config.h.in"
  "${PROJECT_BINARY_DIR}/pumiInterface_config.h"
)
include_directories("${PROJECT_BINARY_DIR}")

# write the path info to be read from Julia
file(WRITE "PUMI_RPATH" "PUMI_CMAKE_PATH ${SCOREC_CONFIG}\n")
file(APPEND "PUMI_RPATH" "PUMIINTERFACE_PREFIX ${CMAKE_INSTALL_PREFIX}\n")
file(APPEND "PUMI_RPATH" "HAVE_SIMMETRIX ${HAVE_SIMMETRIX}")

# set RPATH for binaries
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")



#this is just example code, do your own thing
add_executable(test_init test_init.cc)
if (HAVE_SIMMETRIX)
  add_executable(test_cad test_cad.cc)
  add_executable(convertSimMesh convertSimMesh.cc)
  add_executable(test_init2 test_init2.cc)
endif (HAVE_SIMMETRIX)

add_executable(printGeoEntities printGeoEntities.cc)
add_executable(uniform_refine uniform_refine.cc)
add_executable(test_meshadapt test_meshadapt.cc)
add_library(pumiInterface funcs1.cc a2.cc apfSBPShape.cc apfSBPShape3.cc dgSBPShape1.cc dgSBPShape2.cc dgSBPShape4.cc dgSBPShape5.cc dgSBPShape8.cc dgSBPShape8.cc dgSBP3Shape1.cc dgSBP3Shape2.cc dgSBP3Shape4.cc dgSBP3Shape5.cc triangulation.cc triangulationDG.cc submesh_create.cc)

#for any targets that use PUMI, just use this command
#to it to include the right directories and link to all
#the scorec libraries automatically.
#we recommend PUBLIC if the target is a library and
#PRIVATE if the target is an executable
#target_link_libraries(mylibrary PUBLIC SCOREC::core)
target_link_libraries(pumiInterface PUBLIC SCOREC::core)
target_link_libraries(test_init PRIVATE SCOREC::core)
if (HAVE_SIMMETRIX)
  target_link_libraries(test_init2 PRIVATE SCOREC::core)
  target_link_libraries(test_cad PRIVATE SCOREC::core)
  target_link_libraries(convertSimMesh PRIVATE SCOREC::core)
endif (HAVE_SIMMETRIX)
target_link_libraries(printGeoEntities PRIVATE SCOREC::core)
target_link_libraries(uniform_refine PRIVATE SCOREC::core)
target_link_libraries(test_meshadapt PRIVATE SCOREC::core)

install (TARGETS pumiInterface DESTINATION lib)
install (TARGETS test_init DESTINATION bin)
if(HAVE_SIMMETRIX)
  install (TARGETS test_init2 DESTINATION bin)
  install (TARGETS test_cad DESTINATION bin)
  install (TARGETS convertSimMesh DESTINATION bin)
endif (HAVE_SIMMETRIX)
install (TARGETS printGeoEntities DESTINATION bin)
install (TARGETS uniform_refine DESTINATION bin)
install (TARGETS test_meshadapt DESTINATION bin)

install (FILES funcs1.h a2.h "${PROJECT_BINARY_DIR}/pumiInterface_config.h" DESTINATION include)
