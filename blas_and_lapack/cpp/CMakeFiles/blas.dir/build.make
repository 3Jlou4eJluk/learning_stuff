# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp

# Include any dependencies generated for this target.
include CMakeFiles/blas.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/blas.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/blas.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/blas.dir/flags.make

CMakeFiles/blas.dir/main.o: CMakeFiles/blas.dir/flags.make
CMakeFiles/blas.dir/main.o: main.cpp
CMakeFiles/blas.dir/main.o: CMakeFiles/blas.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/blas.dir/main.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/blas.dir/main.o -MF CMakeFiles/blas.dir/main.o.d -o CMakeFiles/blas.dir/main.o -c /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp/main.cpp

CMakeFiles/blas.dir/main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/blas.dir/main.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp/main.cpp > CMakeFiles/blas.dir/main.i

CMakeFiles/blas.dir/main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/blas.dir/main.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp/main.cpp -o CMakeFiles/blas.dir/main.s

# Object files for target blas
blas_OBJECTS = \
"CMakeFiles/blas.dir/main.o"

# External object files for target blas
blas_EXTERNAL_OBJECTS =

blas: CMakeFiles/blas.dir/main.o
blas: CMakeFiles/blas.dir/build.make
blas: CMakeFiles/blas.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable blas"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/blas.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/blas.dir/build: blas
.PHONY : CMakeFiles/blas.dir/build

CMakeFiles/blas.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/blas.dir/cmake_clean.cmake
.PHONY : CMakeFiles/blas.dir/clean

CMakeFiles/blas.dir/depend:
	cd /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp /run/media/mishanya/MyBestHDD/progs/prakcodesem5/blas_and_lapack/cpp/CMakeFiles/blas.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/blas.dir/depend

