# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.19

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2021.1.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2021.1.2\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\CLionProjects\SBoxDereviankoCPP

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\CLionProjects\SBoxDereviankoCPP\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/SBoxDereviankoCPP.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SBoxDereviankoCPP.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SBoxDereviankoCPP.dir/flags.make

CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.obj: CMakeFiles/SBoxDereviankoCPP.dir/flags.make
CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\CLionProjects\SBoxDereviankoCPP\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.obj"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\SBoxDereviankoCPP.dir\main.cpp.obj -c D:\CLionProjects\SBoxDereviankoCPP\main.cpp

CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.i"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\CLionProjects\SBoxDereviankoCPP\main.cpp > CMakeFiles\SBoxDereviankoCPP.dir\main.cpp.i

CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.s"
	C:\PROGRA~2\MINGW-~1\I686-8~1.0-P\mingw32\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\CLionProjects\SBoxDereviankoCPP\main.cpp -o CMakeFiles\SBoxDereviankoCPP.dir\main.cpp.s

# Object files for target SBoxDereviankoCPP
SBoxDereviankoCPP_OBJECTS = \
"CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.obj"

# External object files for target SBoxDereviankoCPP
SBoxDereviankoCPP_EXTERNAL_OBJECTS =

SBoxDereviankoCPP.exe: CMakeFiles/SBoxDereviankoCPP.dir/main.cpp.obj
SBoxDereviankoCPP.exe: CMakeFiles/SBoxDereviankoCPP.dir/build.make
SBoxDereviankoCPP.exe: CMakeFiles/SBoxDereviankoCPP.dir/linklibs.rsp
SBoxDereviankoCPP.exe: CMakeFiles/SBoxDereviankoCPP.dir/objects1.rsp
SBoxDereviankoCPP.exe: CMakeFiles/SBoxDereviankoCPP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\CLionProjects\SBoxDereviankoCPP\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable SBoxDereviankoCPP.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\SBoxDereviankoCPP.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SBoxDereviankoCPP.dir/build: SBoxDereviankoCPP.exe

.PHONY : CMakeFiles/SBoxDereviankoCPP.dir/build

CMakeFiles/SBoxDereviankoCPP.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\SBoxDereviankoCPP.dir\cmake_clean.cmake
.PHONY : CMakeFiles/SBoxDereviankoCPP.dir/clean

CMakeFiles/SBoxDereviankoCPP.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\CLionProjects\SBoxDereviankoCPP D:\CLionProjects\SBoxDereviankoCPP D:\CLionProjects\SBoxDereviankoCPP\cmake-build-debug D:\CLionProjects\SBoxDereviankoCPP\cmake-build-debug D:\CLionProjects\SBoxDereviankoCPP\cmake-build-debug\CMakeFiles\SBoxDereviankoCPP.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SBoxDereviankoCPP.dir/depend

