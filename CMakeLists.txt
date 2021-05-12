cmake_minimum_required(VERSION 3.10)

# set the project name
project(thermodynamic_integration VERSION 0.1.0 LANGUAGES C)

# paths
include_directories(/opt/local/include)
link_directories(/opt/local/lib)

# external packages
find_package(GSL REQUIRED)

# set cflags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

add_subdirectory(src)

