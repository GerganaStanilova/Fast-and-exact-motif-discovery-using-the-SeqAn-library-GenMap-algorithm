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

# Utility rule file for googlebenchmark.

# Include any custom commands dependencies for this target.
include benchmarks/CMakeFiles/googlebenchmark.dir/compiler_depend.make

# Include the progress variables for this target.
include benchmarks/CMakeFiles/googlebenchmark.dir/progress.make

benchmarks/CMakeFiles/googlebenchmark: benchmarks/CMakeFiles/googlebenchmark-complete

benchmarks/CMakeFiles/googlebenchmark-complete: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-install
benchmarks/CMakeFiles/googlebenchmark-complete: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-mkdir
benchmarks/CMakeFiles/googlebenchmark-complete: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-download
benchmarks/CMakeFiles/googlebenchmark-complete: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-patch
benchmarks/CMakeFiles/googlebenchmark-complete: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-configure
benchmarks/CMakeFiles/googlebenchmark-complete: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-build
benchmarks/CMakeFiles/googlebenchmark-complete: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'googlebenchmark'"
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/CMakeFiles
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/CMakeFiles/googlebenchmark-complete
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-done

benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-build: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Performing build step for 'googlebenchmark'"
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-build && $(MAKE)
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-build && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-build

benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-configure: benchmarks/googlebenchmark/tmp/googlebenchmark-cfgcmd.txt
benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-configure: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Performing configure step for 'googlebenchmark'"
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-build && /snap/cmake/1005/bin/cmake -DBENCHMARK_ENABLE_TESTING=false -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/vagrant/code/motiffinder/genmap-build/benchmarks -DCMAKE_CXX_COMPILER=/usr/bin/c++ "-GUnix Makefiles" /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-build && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-configure

benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-download: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-gitinfo.txt
benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-download: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (git clone) for 'googlebenchmark'"
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src && /snap/cmake/1005/bin/cmake -P /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/tmp/googlebenchmark-gitclone.cmake
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-download

benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-install: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Performing install step for 'googlebenchmark'"
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-build && $(MAKE) install
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-build && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-install

benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Creating directories for 'googlebenchmark'"
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-build
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/tmp
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E make_directory /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-mkdir

benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-patch: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/vagrant/code/motiffinder/genmap-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "No patch step for 'googlebenchmark'"
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E echo_append
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && /snap/cmake/1005/bin/cmake -E touch /home/vagrant/code/motiffinder/genmap-build/benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-patch

googlebenchmark: benchmarks/CMakeFiles/googlebenchmark
googlebenchmark: benchmarks/CMakeFiles/googlebenchmark-complete
googlebenchmark: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-build
googlebenchmark: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-configure
googlebenchmark: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-download
googlebenchmark: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-install
googlebenchmark: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-mkdir
googlebenchmark: benchmarks/googlebenchmark/src/googlebenchmark-stamp/googlebenchmark-patch
googlebenchmark: benchmarks/CMakeFiles/googlebenchmark.dir/build.make
.PHONY : googlebenchmark

# Rule to build all files generated by this target.
benchmarks/CMakeFiles/googlebenchmark.dir/build: googlebenchmark
.PHONY : benchmarks/CMakeFiles/googlebenchmark.dir/build

benchmarks/CMakeFiles/googlebenchmark.dir/clean:
	cd /home/vagrant/code/motiffinder/genmap-build/benchmarks && $(CMAKE_COMMAND) -P CMakeFiles/googlebenchmark.dir/cmake_clean.cmake
.PHONY : benchmarks/CMakeFiles/googlebenchmark.dir/clean

benchmarks/CMakeFiles/googlebenchmark.dir/depend:
	cd /home/vagrant/code/motiffinder/genmap-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/vagrant/code/motiffinder/genmap /home/vagrant/code/motiffinder/genmap/benchmarks /home/vagrant/code/motiffinder/genmap-build /home/vagrant/code/motiffinder/genmap-build/benchmarks /home/vagrant/code/motiffinder/genmap-build/benchmarks/CMakeFiles/googlebenchmark.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : benchmarks/CMakeFiles/googlebenchmark.dir/depend

