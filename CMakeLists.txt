cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(drprobe_clt)
enable_language(Fortran)


# Define the executable name
SET(celslc_exec celslc)

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 



SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-nologo -O3  -132 -assume byterecl -threads -static-intel"  # Intel
                )


############################################################
# Fourier transform
find_package(FFTW3 COMPONENTS double single threads REQUIRED)
############################################################
include_directories(${FFTW3_INCLUDE_DIR})

include_directories("src/common")
include_directories("src/celsc")

#TARGET_LINK_LIBRARIES(${FOOEXE} ${FFTW3_LIBRARIES})

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src/celslc)


# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${SRC})

# Add the source files
SET(CELSLC_src 
	${SRC}/celslc.f90
)


# Modules
add_library(celslcprm OBJECT ${SRC}/celslcprm.f90)
set(OBJLIB_celslcprm_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(emsdata OBJECT ${SRC}/emsdata.f90)
set(OBJLIB_emsdata_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(3dpot OBJECT ${SRC}/3dpot.f90)
set(OBJLIB_3dpot_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(AberrationFunctions OBJECT ${SRC}/AberrationFunctions.f90)
set(OBJLIB_AberrationFunctions_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(CellSlicer OBJECT ${SRC}/CellSlicer.f90)
set(OBJLIB_CellSlicer_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(fitfeprm OBJECT ${SRC}/fitfeprm.f90)
set(OBJLIB_fitfeprm_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(fscatab OBJECT ${SRC}/fscatab.f90)
set(OBJLIB_fscatab_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(cifio OBJECT ${SRC}/cifio.f90)
set(OBJLIB_cifio_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(symops OBJECT ${SRC}/symops.f90)
set(OBJLIB_symops_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)

add_library(spacegroups OBJECT ${SRC}/spacegroups.f90)
set(OBJLIB_spacegroups_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR} PARENT_SCOPE)


# the executable is built from 3 source files
add_executable(
  ${celslc_exec}
  ${CELSLC_src}
$<TARGET_OBJECTS:celslcprm>
$<TARGET_OBJECTS:emsdata>
$<TARGET_OBJECTS:3dpot>
$<TARGET_OBJECTS:AberrationFunctions>
$<TARGET_OBJECTS:CellSlicer>
$<TARGET_OBJECTS:fitfeprm>
$<TARGET_OBJECTS:fscatab>
$<TARGET_OBJECTS:cifio>
$<TARGET_OBJECTS:symops>
$<TARGET_OBJECTS:spacegroups>
 )

target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_celslcprm_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_emsdata_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_3dpot_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_AberrationFunctions_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_CellSlicer_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_fitfeprm_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_fscatab_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_cifio_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_symops_INCLUDE_DIRS} )
target_include_directories(${celslc_exec} PRIVATE ${OBJLIB_spacegroups_INCLUDE_DIRS} )


