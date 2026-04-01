cmake_minimum_required(VERSION 3.20)

add_executable(snuggle_test main.cpp)

# GPU tests need SLIC3R_GUI + OpenGL
option(SNUGGLE_TEST_GPU "Build with GPU evaluator tests" OFF)

if (SNUGGLE_TEST_GPU)
    target_compile_definitions(snuggle_test PRIVATE SLIC3R_GUI)
    target_link_libraries(snuggle_test libslic3r GLEW::GLEW OpenGL::GL)
else()
    target_link_libraries(snuggle_test libslic3r)
endif()

target_include_directories(snuggle_test PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/libslic3r
)

if (WIN32)
    # Copy runtime DLLs next to the test executable
    add_custom_command(TARGET snuggle_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:snuggle_test>
            $<TARGET_FILE_DIR:snuggle_test>
        COMMAND_EXPAND_LISTS
    )
endif()
