# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /snap/cmake/1005/bin/cmake

# The command to remove a file.
RM = /snap/cmake/1005/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/vagrant/code/motiffinder/genmap

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/vagrant/code/motiffinder/genmap-build

# Include any dependencies generated for this target.
include src/CMakeFiles/genmap.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/genmap.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/genmap.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/genmap.dir/flags.make

src/CMakeFiles/genmap.dir/genmap.cpp.o: src/CMakeFiles/genmap.dir/flags.make
src/CMakeFiles/genmap.dir/genmap.cpp.o: /home/vagrant/code/motiffinder/genmap/src/genmap.cpp
src/CMakeFiles/genmap.dir/genmap.cpp.o: src/CMakeFiles/genmap.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/genmap.dir/genmap.cpp.o"
	cd /home/vagrant/code/motiffinder/genmap-build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/genmap.dir/genmap.cpp.o -MF CMakeFiles/genmap.dir/genmap.cpp.o.d -o CMakeFiles/genmap.dir/genmap.cpp.o -c /home/vagrant/code/motiffinder/genmap/src/genmap.cpp

src/CMakeFiles/genmap.dir/genmap.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/genmap.dir/genmap.cpp.i"
	cd /home/vagrant/code/motiffinder/genmap-build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/vagrant/code/motiffinder/genmap/src/genmap.cpp > CMakeFiles/genmap.dir/genmap.cpp.i

src/CMakeFiles/genmap.dir/genmap.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/genmap.dir/genmap.cpp.s"
	cd /home/vagrant/code/motiffinder/genmap-build/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/vagrant/code/motiffinder/genmap/src/genmap.cpp -o CMakeFiles/genmap.dir/genmap.cpp.s

# Object files for target genmap
genmap_OBJECTS = \
"CMakeFiles/genmap.dir/genmap.cpp.o"

# External object files for target genmap
genmap_EXTERNAL_OBJECTS =

bin/genmap: src/CMakeFiles/genmap.dir/genmap.cpp.o
bin/genmap: src/CMakeFiles/genmap.dir/build.make
bin/genmap: src/CMakeFiles/genmap.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/genmap"
	cd /home/vagrant/code/motiffinder/genmap-build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/genmap.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/genmap.dir/build: bin/genmap
.PHONY : src/CMakeFiles/genmap.dir/build

src/CMakeFiles/genmap.dir/clean:
	cd /home/vagrant/code/motiffinder/genmap-build/src && $(CMAKE_COMMAND) -P CMakeFiles/genmap.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/genmap.dir/clean

src/CMakeFiles/genmap.dir/depend:
	cd /home/vagrant/code/motiffinder/genmap-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vagrant/code/motiffinder/genmap /home/vagrant/code/motiffinder/genmap/src /home/vagrant/code/motiffinder/genmap-build /home/vagrant/code/motiffinder/genmap-build/src /home/vagrant/code/motiffinder/genmap-build/src/CMakeFiles/genmap.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/genmap.dir/depend
