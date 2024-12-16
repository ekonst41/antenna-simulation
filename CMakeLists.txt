cmake_minimum_required(VERSION 3.29)
project(bic_simulation)

set(CMAKE_CXX_STANDARD 20)

find_package(Python3 REQUIRED)

add_executable(bic_simulation main.cpp)

get_filename_component(PYTHON_SCRIPT_PATH "${CMAKE_SOURCE_DIR}/visualisation.py" ABSOLUTE)

add_custom_command(
  TARGET bic_simulation POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E echo "Running C++ binary..."
  COMMAND $<TARGET_FILE:bic_simulation> ${CMAKE_SOURCE_DIR}
  COMMAND ${Python3_EXECUTABLE} ${PYTHON_SCRIPT_PATH}
  COMMENT "Visualisation..."
)