# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /usr/local/lib/python2.7/dist-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /usr/local/lib/python2.7/dist-packages/cmake/data/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/buba/Dokumentumok/Doktori/1D_por

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/buba/Dokumentumok/Doktori/1D_por

# Include any dependencies generated for this target.
include CMakeFiles/1D_drift.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/1D_drift.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/1D_drift.dir/flags.make

CMakeFiles/1D_drift.dir/apps/main.cpp.o: CMakeFiles/1D_drift.dir/flags.make
CMakeFiles/1D_drift.dir/apps/main.cpp.o: apps/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/buba/Dokumentumok/Doktori/1D_por/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/1D_drift.dir/apps/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/1D_drift.dir/apps/main.cpp.o -c /home/buba/Dokumentumok/Doktori/1D_por/apps/main.cpp

CMakeFiles/1D_drift.dir/apps/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/1D_drift.dir/apps/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/buba/Dokumentumok/Doktori/1D_por/apps/main.cpp > CMakeFiles/1D_drift.dir/apps/main.cpp.i

CMakeFiles/1D_drift.dir/apps/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/1D_drift.dir/apps/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/buba/Dokumentumok/Doktori/1D_por/apps/main.cpp -o CMakeFiles/1D_drift.dir/apps/main.cpp.s

# Object files for target 1D_drift
1D_drift_OBJECTS = \
"CMakeFiles/1D_drift.dir/apps/main.cpp.o"

# External object files for target 1D_drift
1D_drift_EXTERNAL_OBJECTS =

1D_drift: CMakeFiles/1D_drift.dir/apps/main.cpp.o
1D_drift: CMakeFiles/1D_drift.dir/build.make
1D_drift: lib/libdustdrift.so
1D_drift: CMakeFiles/1D_drift.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/buba/Dokumentumok/Doktori/1D_por/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable 1D_drift"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/1D_drift.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/1D_drift.dir/build: 1D_drift

.PHONY : CMakeFiles/1D_drift.dir/build

CMakeFiles/1D_drift.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/1D_drift.dir/cmake_clean.cmake
.PHONY : CMakeFiles/1D_drift.dir/clean

CMakeFiles/1D_drift.dir/depend:
	cd /home/buba/Dokumentumok/Doktori/1D_por && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/buba/Dokumentumok/Doktori/1D_por /home/buba/Dokumentumok/Doktori/1D_por /home/buba/Dokumentumok/Doktori/1D_por /home/buba/Dokumentumok/Doktori/1D_por /home/buba/Dokumentumok/Doktori/1D_por/CMakeFiles/1D_drift.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/1D_drift.dir/depend

