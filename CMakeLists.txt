cmake_minimum_required(VERSION 3.15...3.26)
project(ext)

set(PYTHON "ON")
set(OSQP_BUILD_UNITTESTS "OFF")
set(OSQP_USE_LONG "OFF")
set(OSQP_CUSTOM_PRINTING "${CMAKE_CURRENT_SOURCE_DIR}/cmake/printing.h")
set(OSQP_CUSTOM_MEMORY "${CMAKE_CURRENT_SOURCE_DIR}/cmake/memory.h")
set(OSQP_CODEGEN_INSTALL_DIR "codegen/codegen_src" CACHE PATH "" FORCE)

if(APPLE)
    message(STATUS "Building for Apple arches: ${CMAKE_OSX_ARCHITECTURES}")
endif()

include(FetchContent)

# 03/05/24 - Use modern python discovery
set(PYBIND11_FINDPYTHON "ON")

# Try to find pybind11, if not found, fetch it
find_package(pybind11 CONFIG QUIET)
if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found, fetching from GitHub...")
    FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG v2.11.1
    )
    FetchContent_MakeAvailable(pybind11)
else()
    message(STATUS "Found pybind11: ${pybind11_DIR}")
endif()

# find_package(pybind11 CONFIG REQUIRED)

# Find Python before configuring OSQP
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# 03/05/24 - Workaround because OSQP CMakeLists.txt is using old variable names
set(PYTHON_FOUND "ON")
set(PYTHON_INCLUDE_DIRS ${Python_INCLUDE_DIRS})

message(STATUS "Fetching/configuring OSQP")
list(APPEND CMAKE_MESSAGE_INDENT "  ")
FetchContent_Declare(
  osqp
#   GIT_REPOSITORY https://github.com/osqp/osqp.git
#   GIT_REPOSITORY https://github.com/lb3825/osqp.git
    SOURCE_DIR "C:/Users/baice/Documents/GitHub/osqp"
#   GIT_TAG v1.0.0
#   GIT_TAG restarts
#   GIT_TAG plot
)
list(POP_BACK CMAKE_MESSAGE_INDENT)
FetchContent_MakeAvailable(osqp)

# Define the module name if not already set
if(NOT DEFINED OSQP_EXT_MODULE_NAME)
    set(OSQP_EXT_MODULE_NAME "osqp")
endif()

if(${OSQP_ALGEBRA_BACKEND} STREQUAL "builtin")
    set(OSQP_EXT_MODULE_NAME "ext_builtin")
elseif(${OSQP_ALGEBRA_BACKEND} STREQUAL "mkl")
    set(OSQP_EXT_MODULE_NAME "osqp_mkl")
elseif(${OSQP_ALGEBRA_BACKEND} STREQUAL "cuda")
    set(OSQP_EXT_MODULE_NAME "osqp_cuda")
else()
    set(OSQP_EXT_MODULE_NAME "osqp")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/bindings.cpp.in
               ${CMAKE_CURRENT_SOURCE_DIR}/src/bindings.cpp)
pybind11_add_module(${OSQP_EXT_MODULE_NAME} src/bindings.cpp)
install(TARGETS ${OSQP_EXT_MODULE_NAME} DESTINATION . COMPONENT python)

# TODO: We shouldn't have to do this once the interfaces are set up correctly
if(${OSQP_ALGEBRA_BACKEND} STREQUAL "builtin")
    target_link_libraries(ext_builtin PUBLIC pybind11::module osqpstatic)
elseif(${OSQP_ALGEBRA_BACKEND} STREQUAL "mkl")
    if(APPLE)
        target_link_libraries(osqp_mkl PUBLIC pybind11::module osqpstatic)
    else()
        target_link_libraries(osqp_mkl PUBLIC pybind11::module osqpstatic $<LINK_ONLY:MKL::MKL>)
    endif()
elseif(${OSQP_ALGEBRA_BACKEND} STREQUAL "cuda")
    enable_language(CUDA)
    find_package(CUDA)
    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_directories(osqp_cuda PUBLIC ${CMAKE_CUDA_HOST_IMPLICIT_LINK_DIRECTORIES})
    target_link_libraries(osqp_cuda PUBLIC pybind11::module osqpstatic cublas cusparse)
endif()
