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
include CMakeFiles/3dpot.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/3dpot.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/3dpot.dir/flags.make

CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o: CMakeFiles/3dpot.dir/flags.make
CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o: ../src/celslc/3dpot.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/3dpot.f90 -o CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o

CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.i"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/3dpot.f90 > CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.i

CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.s"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/3dpot.f90 -o CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.s

CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.requires:

.PHONY : CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.requires

CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.provides: CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.requires
	$(MAKE) -f CMakeFiles/3dpot.dir/build.make CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.provides.build
.PHONY : CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.provides

CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.provides.build: CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o


3dpot: CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o
3dpot: CMakeFiles/3dpot.dir/build.make

.PHONY : 3dpot

# Rule to build all files generated by this target.
CMakeFiles/3dpot.dir/build: 3dpot

.PHONY : CMakeFiles/3dpot.dir/build

CMakeFiles/3dpot.dir/requires: CMakeFiles/3dpot.dir/src/celslc/3dpot.f90.o.requires

.PHONY : CMakeFiles/3dpot.dir/requires

CMakeFiles/3dpot.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/3dpot.dir/cmake_clean.cmake
.PHONY : CMakeFiles/3dpot.dir/clean

CMakeFiles/3dpot.dir/depend:
	cd /Users/filipeoliveira/Documents/drprobe_clt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/filipeoliveira/Documents/drprobe_clt /Users/filipeoliveira/Documents/drprobe_clt /Users/filipeoliveira/Documents/drprobe_clt/build /Users/filipeoliveira/Documents/drprobe_clt/build /Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/3dpot.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/3dpot.dir/depend

