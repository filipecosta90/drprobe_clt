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
include CMakeFiles/symops.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/symops.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/symops.dir/flags.make

CMakeFiles/symops.dir/src/celslc/symops.f90.o: CMakeFiles/symops.dir/flags.make
CMakeFiles/symops.dir/src/celslc/symops.f90.o: ../src/celslc/symops.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/symops.dir/src/celslc/symops.f90.o"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/symops.f90 -o CMakeFiles/symops.dir/src/celslc/symops.f90.o

CMakeFiles/symops.dir/src/celslc/symops.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/symops.dir/src/celslc/symops.f90.i"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/symops.f90 > CMakeFiles/symops.dir/src/celslc/symops.f90.i

CMakeFiles/symops.dir/src/celslc/symops.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/symops.dir/src/celslc/symops.f90.s"
	/usr/local/bin/ifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /Users/filipeoliveira/Documents/drprobe_clt/src/celslc/symops.f90 -o CMakeFiles/symops.dir/src/celslc/symops.f90.s

CMakeFiles/symops.dir/src/celslc/symops.f90.o.requires:

.PHONY : CMakeFiles/symops.dir/src/celslc/symops.f90.o.requires

CMakeFiles/symops.dir/src/celslc/symops.f90.o.provides: CMakeFiles/symops.dir/src/celslc/symops.f90.o.requires
	$(MAKE) -f CMakeFiles/symops.dir/build.make CMakeFiles/symops.dir/src/celslc/symops.f90.o.provides.build
.PHONY : CMakeFiles/symops.dir/src/celslc/symops.f90.o.provides

CMakeFiles/symops.dir/src/celslc/symops.f90.o.provides.build: CMakeFiles/symops.dir/src/celslc/symops.f90.o


symops: CMakeFiles/symops.dir/src/celslc/symops.f90.o
symops: CMakeFiles/symops.dir/build.make

.PHONY : symops

# Rule to build all files generated by this target.
CMakeFiles/symops.dir/build: symops

.PHONY : CMakeFiles/symops.dir/build

CMakeFiles/symops.dir/requires: CMakeFiles/symops.dir/src/celslc/symops.f90.o.requires

.PHONY : CMakeFiles/symops.dir/requires

CMakeFiles/symops.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/symops.dir/cmake_clean.cmake
.PHONY : CMakeFiles/symops.dir/clean

CMakeFiles/symops.dir/depend:
	cd /Users/filipeoliveira/Documents/drprobe_clt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/filipeoliveira/Documents/drprobe_clt /Users/filipeoliveira/Documents/drprobe_clt /Users/filipeoliveira/Documents/drprobe_clt/build /Users/filipeoliveira/Documents/drprobe_clt/build /Users/filipeoliveira/Documents/drprobe_clt/build/CMakeFiles/symops.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/symops.dir/depend

