# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-src")
  file(MAKE_DIRECTORY "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-src")
endif()
file(MAKE_DIRECTORY
  "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-build"
  "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-subbuild/qdldl-populate-prefix"
  "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-subbuild/qdldl-populate-prefix/tmp"
  "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-subbuild/qdldl-populate-prefix/src/qdldl-populate-stamp"
  "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-subbuild/qdldl-populate-prefix/src"
  "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-subbuild/qdldl-populate-prefix/src/qdldl-populate-stamp"
)

set(configSubDirs Debug)
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-subbuild/qdldl-populate-prefix/src/qdldl-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "C:/Users/baice/Documents/GitHub/osqp-python/_deps/qdldl-subbuild/qdldl-populate-prefix/src/qdldl-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
