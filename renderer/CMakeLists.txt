file(GLOB_RECURSE SRCS *.cpp *.c  *.cu *.h )
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#add_executable(SampleBrowser WIN32  ${SRCS})
add_executable(Rasterizer   ${SRCS})

target_include_directories(Rasterizer  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Rasterizer  PUBLIC inc)
target_link_libraries(Rasterizer
PUBLIC
OpenMP::OpenMP_CXX
)
target_compile_features(Rasterizer PRIVATE cxx_std_17)
