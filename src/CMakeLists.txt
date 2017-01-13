#This file shows how to link to a PUMI
#installation using CMake
#it represents a simple 'CMakeLists.txt'
#file for a new project

cmake_minimum_required(VERSION 3.0.0)

project(MeshCreate VERSION 1.0.0 LANGUAGES CXX)

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

#this is just example code, do your own thing
#add_executable(mylibrary mylibrary.cpp)
add_library(funcs1 funcs1.cc a2.cc adaptFuncsJ.cc apfSBPShape.cc apfSBPShape3.cc dgSBPShape1.cc dgSBPShape2.cc dgSBP3Shape1.cc dgSBP3Shape2.cc triangulation.cc triangulationDG.cc)

#for any targets that use PUMI, just use this command
#to it to include the right directories and link to all
#the scorec libraries automatically.
#we recommend PUBLIC if the target is a library and
#PRIVATE if the target is an executable
#target_link_libraries(mylibrary PUBLIC SCOREC::core)
target_link_libraries(funcs1 PUBLIC SCOREC::core)

install (TARGETS funcs1 DESTINATION lib)