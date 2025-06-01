cmake_minimum_required(VERSION 3.27 FATAL_ERROR)
project(vp-solver)

# Make Ensign available
include(FetchContent)
FetchContent_Declare(Ensign
  GIT_REPOSITORY https://github.com/leinkemmer/Ensign.git
  GIT_TAG "development"
  GIT_SHALLOW ON
  GIT_PROGRESS ON
  FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(Ensign)

find_package(OpenMP)

add_executable(vp-solver main_template.cpp)
target_link_libraries(vp-solver PRIVATE Ensign::Ensign)
target_link_libraries(vp-solver PRIVATE OpenMP::OpenMP_CXX)