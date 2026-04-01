cmake_minimum_required(VERSION 3.20)

add_executable(snuggle_test main.cpp)

target_include_directories(snuggle_test PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/libslic3r
)

# libslic3r propagates SLIC3R_GUI, so GPU code is always visible.
# Link GL libs for the GPU evaluator.
target_link_libraries(snuggle_test libslic3r)

if (TARGET GLEW::GLEW)
    target_link_libraries(snuggle_test GLEW::GLEW)
endif()
if (TARGET OpenGL::GL)
    target_link_libraries(snuggle_test OpenGL::GL)
endif()
