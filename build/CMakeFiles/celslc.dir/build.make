# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.10.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.10.0/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/filipeoliveira/Documents/drprobe_clt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/filipeoliveira/Documents/drprobe_clt/build

# Include any dependencies generated for this target.
include CMakeFiles/celslc.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/celslc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/celslc.dir/flags.make

CMakeFiles/celslc.dir/src/celslc/celslc.f90.o: CMakeFiles/celslc.dir/flags.make
CMakeFiles/celslc.dir/src/celslc/celslc.f90.o: ../src/celslc/celslc.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/celslc.dir/src/celslc/celslc.f90.o"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/celslc.f90 -o CMakeFiles/celslc.dir/src/celslc/celslc.f90.o

CMakeFiles/celslc.dir/src/celslc/celslc.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/celslc.dir/src/celslc/celslc.f90.i"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/celslc.f90 > CMakeFiles/celslc.dir/src/celslc/celslc.f90.i

CMakeFiles/celslc.dir/src/celslc/celslc.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/celslc.dir/src/celslc/celslc.f90.s"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/celslc.f90 -o CMakeFiles/celslc.dir/src/celslc/celslc.f90.s

CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.requires:

.PHONY : CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.requires

CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.provides: CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.requires
	$(MAKE) -f CMakeFiles/celslc.dir/build.make CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.provides.build
.PHONY : CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.provides

CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.provides.build: CMakeFiles/celslc.dir/src/celslc/celslc.f90.o


# Object files for target celslc
celslc_OBJECTS = \
"CMakeFiles/celslc.dir/src/celslc/celslc.f90.o"

# External object files for target celslc
celslc_EXTERNAL_OBJECTS = \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/celslcprm.dir/src/celslc/celslcprm.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/emsdata.dir/src/celslc/emsdata.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/AberrationFunctions.dir/src/celslc/AberrationFunctions.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/CellSlicer.dir/src/celslc/CellSlicer.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/fitfeprm.dir/src/celslc/fitfeprm.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/fscatab.dir/src/celslc/fscatab.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/cifio.dir/src/celslc/cifio.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/symops.dir/src/celslc/symops.f90.o" \
"/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/spacegroups.dir/src/celslc/spacegroups.f90.o"

celslc: CMakeFiles/celslc.dir/src/celslc/celslc.f90.o
celslc: CMakeFiles/celslcprm.dir/src/celslc/celslcprm.f90.o
celslc: CMakeFiles/emsdata.dir/src/celslc/emsdata.f90.o
celslc: CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o
celslc: CMakeFiles/AberrationFunctions.dir/src/celslc/AberrationFunctions.f90.o
celslc: CMakeFiles/CellSlicer.dir/src/celslc/CellSlicer.f90.o
celslc: CMakeFiles/fitfeprm.dir/src/celslc/fitfeprm.f90.o
celslc: CMakeFiles/fscatab.dir/src/celslc/fscatab.f90.o
celslc: CMakeFiles/cifio.dir/src/celslc/cifio.f90.o
celslc: CMakeFiles/symops.dir/src/celslc/symops.f90.o
celslc: CMakeFiles/spacegroups.dir/src/celslc/spacegroups.f90.o
celslc: CMakeFiles/celslc.dir/build.make
celslc: CMakeFiles/celslc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable celslc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/celslc.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/celslc.dir/build: celslc

.PHONY : CMakeFiles/celslc.dir/build

CMakeFiles/celslc.dir/requires: CMakeFiles/celslc.dir/src/celslc/celslc.f90.o.requires

.PHONY : CMakeFiles/celslc.dir/requires

CMakeFiles/celslc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/celslc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/celslc.dir/clean

CMakeFiles/celslc.dir/depend:
	cd /Users/filipeoliveira/Documents/drprobe_clt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/filipeoliveira/Documents/drprobe_clt /Users/filipeoliveira/Documents/drprobe_clt /Users/filipeoliveira/Documents/drprobe_clt/build /Users/filipeoliveira/Documents/drprobe_clt/build /Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/celslc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/celslc.dir/depend
