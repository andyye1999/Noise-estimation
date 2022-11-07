cmake_minimum_required(VERSION 2.8)
project(OMLSA LANGUAGES C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET(CMAKE_BUILD_TYPE "Release")
add_executable(OMLSA
        expint.c
        expint.h
        main.c
        anc.c
        speech_enhance.h)
target_link_libraries(OMLSA m)
