# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/thomas/Code/SuperLUbench

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/thomas/Code/SuperLUbench

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/thomas/Code/SuperLUbench/CMakeFiles /home/thomas/Code/SuperLUbench/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/thomas/Code/SuperLUbench/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named SuperLU_Dist_bench

# Build rule for target.
SuperLU_Dist_bench: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 SuperLU_Dist_bench
.PHONY : SuperLU_Dist_bench

# fast build rule for target.
SuperLU_Dist_bench/fast:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/build
.PHONY : SuperLU_Dist_bench/fast

#=============================================================================
# Target rules for targets named SuperLU_bench

# Build rule for target.
SuperLU_bench: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 SuperLU_bench
.PHONY : SuperLU_bench

# fast build rule for target.
SuperLU_bench/fast:
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/build
.PHONY : SuperLU_bench/fast

#=============================================================================
# Target rules for targets named SuperLU_MT_bench

# Build rule for target.
SuperLU_MT_bench: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 SuperLU_MT_bench
.PHONY : SuperLU_MT_bench

# fast build rule for target.
SuperLU_MT_bench/fast:
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/build
.PHONY : SuperLU_MT_bench/fast

src/SuperLU/dlinsolx.o: src/SuperLU/dlinsolx.c.o

.PHONY : src/SuperLU/dlinsolx.o

# target to build an object file
src/SuperLU/dlinsolx.c.o:
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/SuperLU/dlinsolx.c.o
.PHONY : src/SuperLU/dlinsolx.c.o

src/SuperLU/dlinsolx.i: src/SuperLU/dlinsolx.c.i

.PHONY : src/SuperLU/dlinsolx.i

# target to preprocess a source file
src/SuperLU/dlinsolx.c.i:
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/SuperLU/dlinsolx.c.i
.PHONY : src/SuperLU/dlinsolx.c.i

src/SuperLU/dlinsolx.s: src/SuperLU/dlinsolx.c.s

.PHONY : src/SuperLU/dlinsolx.s

# target to generate assembly for a file
src/SuperLU/dlinsolx.c.s:
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/SuperLU/dlinsolx.c.s
.PHONY : src/SuperLU/dlinsolx.c.s

src/SuperLU_Dist/pddrive_ABglobal.o: src/SuperLU_Dist/pddrive_ABglobal.c.o

.PHONY : src/SuperLU_Dist/pddrive_ABglobal.o

# target to build an object file
src/SuperLU_Dist/pddrive_ABglobal.c.o:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/SuperLU_Dist/pddrive_ABglobal.c.o
.PHONY : src/SuperLU_Dist/pddrive_ABglobal.c.o

src/SuperLU_Dist/pddrive_ABglobal.i: src/SuperLU_Dist/pddrive_ABglobal.c.i

.PHONY : src/SuperLU_Dist/pddrive_ABglobal.i

# target to preprocess a source file
src/SuperLU_Dist/pddrive_ABglobal.c.i:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/SuperLU_Dist/pddrive_ABglobal.c.i
.PHONY : src/SuperLU_Dist/pddrive_ABglobal.c.i

src/SuperLU_Dist/pddrive_ABglobal.s: src/SuperLU_Dist/pddrive_ABglobal.c.s

.PHONY : src/SuperLU_Dist/pddrive_ABglobal.s

# target to generate assembly for a file
src/SuperLU_Dist/pddrive_ABglobal.c.s:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/SuperLU_Dist/pddrive_ABglobal.c.s
.PHONY : src/SuperLU_Dist/pddrive_ABglobal.c.s

src/SuperLU_MT/pdlinsolx.o: src/SuperLU_MT/pdlinsolx.c.o

.PHONY : src/SuperLU_MT/pdlinsolx.o

# target to build an object file
src/SuperLU_MT/pdlinsolx.c.o:
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/SuperLU_MT/pdlinsolx.c.o
.PHONY : src/SuperLU_MT/pdlinsolx.c.o

src/SuperLU_MT/pdlinsolx.i: src/SuperLU_MT/pdlinsolx.c.i

.PHONY : src/SuperLU_MT/pdlinsolx.i

# target to preprocess a source file
src/SuperLU_MT/pdlinsolx.c.i:
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/SuperLU_MT/pdlinsolx.c.i
.PHONY : src/SuperLU_MT/pdlinsolx.c.i

src/SuperLU_MT/pdlinsolx.s: src/SuperLU_MT/pdlinsolx.c.s

.PHONY : src/SuperLU_MT/pdlinsolx.s

# target to generate assembly for a file
src/SuperLU_MT/pdlinsolx.c.s:
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/SuperLU_MT/pdlinsolx.c.s
.PHONY : src/SuperLU_MT/pdlinsolx.c.s

src/Util.o: src/Util.c.o

.PHONY : src/Util.o

# target to build an object file
src/Util.c.o:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/Util.c.o
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/Util.c.o
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/Util.c.o
.PHONY : src/Util.c.o

src/Util.i: src/Util.c.i

.PHONY : src/Util.i

# target to preprocess a source file
src/Util.c.i:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/Util.c.i
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/Util.c.i
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/Util.c.i
.PHONY : src/Util.c.i

src/Util.s: src/Util.c.s

.PHONY : src/Util.s

# target to generate assembly for a file
src/Util.c.s:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/Util.c.s
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/Util.c.s
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/Util.c.s
.PHONY : src/Util.c.s

src/dreadMM.o: src/dreadMM.c.o

.PHONY : src/dreadMM.o

# target to build an object file
src/dreadMM.c.o:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/dreadMM.c.o
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/dreadMM.c.o
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/dreadMM.c.o
.PHONY : src/dreadMM.c.o

src/dreadMM.i: src/dreadMM.c.i

.PHONY : src/dreadMM.i

# target to preprocess a source file
src/dreadMM.c.i:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/dreadMM.c.i
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/dreadMM.c.i
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/dreadMM.c.i
.PHONY : src/dreadMM.c.i

src/dreadMM.s: src/dreadMM.c.s

.PHONY : src/dreadMM.s

# target to generate assembly for a file
src/dreadMM.c.s:
	$(MAKE) -f CMakeFiles/SuperLU_Dist_bench.dir/build.make CMakeFiles/SuperLU_Dist_bench.dir/src/dreadMM.c.s
	$(MAKE) -f CMakeFiles/SuperLU_bench.dir/build.make CMakeFiles/SuperLU_bench.dir/src/dreadMM.c.s
	$(MAKE) -f CMakeFiles/SuperLU_MT_bench.dir/build.make CMakeFiles/SuperLU_MT_bench.dir/src/dreadMM.c.s
.PHONY : src/dreadMM.c.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... SuperLU_Dist_bench"
	@echo "... SuperLU_bench"
	@echo "... SuperLU_MT_bench"
	@echo "... src/SuperLU/dlinsolx.o"
	@echo "... src/SuperLU/dlinsolx.i"
	@echo "... src/SuperLU/dlinsolx.s"
	@echo "... src/SuperLU_Dist/pddrive_ABglobal.o"
	@echo "... src/SuperLU_Dist/pddrive_ABglobal.i"
	@echo "... src/SuperLU_Dist/pddrive_ABglobal.s"
	@echo "... src/SuperLU_MT/pdlinsolx.o"
	@echo "... src/SuperLU_MT/pdlinsolx.i"
	@echo "... src/SuperLU_MT/pdlinsolx.s"
	@echo "... src/Util.o"
	@echo "... src/Util.i"
	@echo "... src/Util.s"
	@echo "... src/dreadMM.o"
	@echo "... src/dreadMM.i"
	@echo "... src/dreadMM.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

