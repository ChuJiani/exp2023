file(GLOB_RECURSE sources CONFIGURE_DEPENDS src/*.cpp include/*.h)
add_executable(align ${sources})
target_include_directories(align PUBLIC include)
