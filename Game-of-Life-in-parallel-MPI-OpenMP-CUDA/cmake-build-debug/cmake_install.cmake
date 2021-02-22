# Install script for directory: /Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/game_mpi.exe")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA" TYPE EXECUTABLE FILES "/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/cmake-build-debug/game_mpi.exe")
  if(EXISTS "$ENV{DESTDIR}/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/game_mpi.exe" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/game_mpi.exe")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/-L/usr/local/Cellar/libevent/2.1.12/lib"
      "$ENV{DESTDIR}/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/game_mpi.exe")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/game_mpi.exe")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/xiejing/Game-of-Life-in-parallel-MPI-OpenMP-CUDA/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
