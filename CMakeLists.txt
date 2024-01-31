################################################################################
cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(OpenGLMINI LANGUAGES CXX CUDA)
#######################bin folder############################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXTERNAL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

#######################bin folder############################################
#compiler stuff
IF(MSVC)
  ADD_DEFINITIONS (/D  _CRT_SECURE_NO_WARNINGS)  ## add macro
 # ADD_DEFINITIONS (/D  _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING )
  #ADD_DEFINITIONS (/D   GLEW_STATIC ) ## use static glew library
  #ADD_DEFINITIONS (/D   GLFW_INCLUDE_NONE ) ## use static glew library
ENDIF()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# package
find_package(OpenMP REQUIRED)
add_subdirectory(renderer)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakePredefinedTargetsFolder")

