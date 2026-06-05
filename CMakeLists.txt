cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(displacement-analysis VERSION 1.0.0 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)

volt_add_plugin(displacement-analysis
    SOURCES src/displacements_engine.cpp
    MAIN_SOURCE src/plugin.cpp
)
