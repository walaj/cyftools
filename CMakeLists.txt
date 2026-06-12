cmake_minimum_required(VERSION 3.16)
project(cyftools
  VERSION 0.1.0
  DESCRIPTION "Command-line toolkit for multiplex-imaging single-cell tables (CYF format)"
  LANGUAGES CXX)

#
# cyftools — modern CMake build.
#
# The format core plus the streaming and spatial commands need only a C++17
# compiler and header-only libraries (cereal + nanoflann + delaunator are
# vendored/submoduled). Image output, plotting, neighbor embedding, and LDA are
# optional features behind the CYFTOOLS_WITH_* options below; turning them all
# off yields a build that needs no system libraries at all.
#
# Quick start:
#   git submodule update --init
#   cmake -S . -B build && cmake --build build -j
#   ctest --test-dir build              # run the format self-test
#
# Minimal core build (no system deps):
#   cmake -S . -B build -DCYFTOOLS_WITH_TIFF=OFF -DCYFTOOLS_WITH_CAIRO=OFF
#

# ---------------------------------------------------------------------------
# Standard / defaults
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Options (each maps to a HAVE_* compile definition the sources already gate on)
# ---------------------------------------------------------------------------
option(CYFTOOLS_64BIT        "Use 64-bit flag registers"                ON)
option(CYFTOOLS_WITH_OPENMP  "Parallelism via OpenMP"                   ON)
option(CYFTOOLS_WITH_TIFF    "TIFF image I/O (convolve, png)"           ON)
option(CYFTOOLS_WITH_CAIRO   "Cairo PDF/PNG plotting"                   ON)
option(CYFTOOLS_WITH_LDA     "LDA topic modeling (needs supervised-lda)" OFF)
option(CYFTOOLS_BUILD_TESTS  "Build the cyf format self-test"           ON)

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(EXT ${CMAKE_CURRENT_SOURCE_DIR}/external)

# ---------------------------------------------------------------------------
# Sources. Every file is compiled unconditionally — the heavy/optional code in
# each is guarded internally by the HAVE_* macros, so a file with its feature
# disabled simply compiles to little or nothing.
# ---------------------------------------------------------------------------
set(CYF_SOURCES
  ${SRC}/cyftools.cpp
  ${SRC}/cell_table.cpp
  ${SRC}/cell_table_cluster.cpp     # DBSCAN on nanoflann (was cell_table_mlpack)
  ${SRC}/cell_table_delaunay.cpp
  ${SRC}/cell_table_graph.cpp
  ${SRC}/cell_table_lda.cpp
  ${SRC}/cell_header.cpp
  ${SRC}/cell_graph.cpp
  ${SRC}/cell_flag.cpp
  ${SRC}/cell_utils.cpp
  ${SRC}/cell_processor.cpp
  ${SRC}/cell_row.cpp
  ${SRC}/cell_selector.cpp
  ${SRC}/cell_synth.cpp
  ${SRC}/cyf_io.cpp
  ${SRC}/bgzf.cpp
  ${SRC}/polygon.cpp
  ${SRC}/color_map.cpp
  ${SRC}/tiff_reader.cpp
  ${SRC}/tiff_writer.cpp
  ${SRC}/tiff_header.cpp
  ${SRC}/tiff_utils.cpp
  ${SRC}/tiff_ifd.cpp
  ${SRC}/tiff_image.cpp
  ${SRC}/tiff_cp.cpp
)

add_executable(cyftools ${CYF_SOURCES})
target_include_directories(cyftools PRIVATE
  ${SRC}                       # nanoflann.hpp + delaunator.hpp live here
  ${EXT}/cereal/include)       # cereal: still a core dep until the CYF migration completes

if(CYFTOOLS_64BIT)
  target_compile_definitions(cyftools PRIVATE USE_64_BIT)
endif()

# threads: the CSV reader (csv.h) uses a background std::thread worker. On GCC/Clang
# pass -pthread at BOTH compile and link (a bare -lpthread can deadlock std::thread);
# elsewhere use the Threads package.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(cyftools PRIVATE -pthread)
  target_link_options(cyftools PRIVATE -pthread)
else()
  find_package(Threads REQUIRED)
  target_link_libraries(cyftools PRIVATE Threads::Threads)
endif()

# zlib: the BGZF (.bcyf) compression container is built on it. Tiny, universal.
find_package(ZLIB REQUIRED)
target_link_libraries(cyftools PRIVATE ZLIB::ZLIB)

# ---------------------------------------------------------------------------
# OpenMP
# ---------------------------------------------------------------------------
if(CYFTOOLS_WITH_OPENMP)
  find_package(OpenMP)
  set(_omp_handled FALSE)
  if(NOT OpenMP_CXX_FOUND AND APPLE)
    # Apple's clang ships no OpenMP runtime — wire up Homebrew's libomp if present.
    find_program(_brew brew)
    if(_brew)
      execute_process(COMMAND ${_brew} --prefix libomp
                      OUTPUT_VARIABLE _libomp OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    endif()
    if(_libomp AND EXISTS "${_libomp}/include/omp.h")
      message(STATUS "OpenMP: using Homebrew libomp at ${_libomp}")
      target_compile_options(cyftools PRIVATE -Xpreprocessor -fopenmp -I${_libomp}/include)
      target_link_libraries(cyftools PRIVATE ${_libomp}/lib/libomp.dylib)
      target_compile_definitions(cyftools PRIVATE HAVE_OMP)
      set(_omp_handled TRUE)
    endif()
  endif()
  if(OpenMP_CXX_FOUND)
    target_link_libraries(cyftools PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(cyftools PRIVATE HAVE_OMP)
  elseif(NOT _omp_handled)
    message(WARNING "OpenMP not found; building single-threaded "
                    "(on macOS try 'brew install libomp', or set -DCYFTOOLS_WITH_OPENMP=OFF)")
  endif()
endif()

# ---------------------------------------------------------------------------
# knncolle: header-only KNN search used by the island / dist / tls commands.
# (Eigen-free; the UMAP feature and its umappp/CppIrlba/CppKmeans deps were removed.)
# ---------------------------------------------------------------------------
if(NOT EXISTS "${EXT}/knncolle/include/knncolle/knncolle.hpp")
  message(FATAL_ERROR
    "The knncolle submodule is not checked out.\n"
    "  Fix:  git submodule update --init")
endif()
target_include_directories(cyftools PRIVATE ${EXT}/knncolle/include)
target_compile_definitions(cyftools PRIVATE
  HAVE_KNNCOLLE KNNCOLLE_NO_KMKNN KNNCOLLE_NO_HNSW KNNCOLLE_NO_ANNOY)

# ---------------------------------------------------------------------------
# TIFF + JPEG (image I/O)
# ---------------------------------------------------------------------------
if(CYFTOOLS_WITH_TIFF)
  find_package(TIFF REQUIRED)
  find_package(JPEG REQUIRED)
  target_link_libraries(cyftools PRIVATE TIFF::TIFF JPEG::JPEG)
  target_compile_definitions(cyftools PRIVATE HAVE_TIFFLIB)
endif()

# ---------------------------------------------------------------------------
# Cairo (PDF/PNG plotting) — located via pkg-config
# ---------------------------------------------------------------------------
if(CYFTOOLS_WITH_CAIRO)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(CAIRO REQUIRED IMPORTED_TARGET cairo)
  target_link_libraries(cyftools PRIVATE PkgConfig::CAIRO)
  target_compile_definitions(cyftools PRIVATE HAVE_CAIRO)
endif()

# ---------------------------------------------------------------------------
# LDA (supervised-lda / ldaplusplus) — opt-in
# ---------------------------------------------------------------------------
if(CYFTOOLS_WITH_LDA)
  set(LDA_ROOT "$ENV{HOME}/git/supervised-lda" CACHE PATH "supervised-lda install root")
  find_path(LDA_INCLUDE_DIR ldaplusplus/LDA.hpp HINTS ${LDA_ROOT}/include)
  find_library(LDA_LIBRARY NAMES ldaplusplus HINTS ${LDA_ROOT}/lib ${LDA_ROOT}/build)
  if(NOT LDA_INCLUDE_DIR OR NOT LDA_LIBRARY)
    message(FATAL_ERROR "CYFTOOLS_WITH_LDA=ON but supervised-lda not found under ${LDA_ROOT}")
  endif()
  find_package(Eigen3 REQUIRED)   # LDA needs Eigen
  target_link_libraries(cyftools PRIVATE Eigen3::Eigen)
  target_include_directories(cyftools PRIVATE ${LDA_INCLUDE_DIR})
  target_link_libraries(cyftools PRIVATE ${LDA_LIBRARY})
  target_compile_definitions(cyftools PRIVATE HAVE_LDAPLUSPLUS)
endif()

install(TARGETS cyftools RUNTIME DESTINATION bin)

# ---------------------------------------------------------------------------
# Format self-test (dependency-light: only cereal). `ctest` runs it.
# ---------------------------------------------------------------------------
if(CYFTOOLS_BUILD_TESTS)
  enable_testing()
  add_executable(cyf_test ${SRC}/cyf_test.cpp ${SRC}/cyf_io.cpp ${SRC}/cell_header.cpp ${SRC}/bgzf.cpp)
  target_include_directories(cyf_test PRIVATE ${SRC} ${EXT}/cereal/include)
  target_link_libraries(cyf_test PRIVATE ZLIB::ZLIB)
  if(CYFTOOLS_64BIT)
    target_compile_definitions(cyf_test PRIVATE USE_64_BIT)
  endif()
  add_test(NAME cyf_test COMMAND cyf_test)
endif()

# ---------------------------------------------------------------------------
# Config summary
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "cyftools ${PROJECT_VERSION} — build configuration:")
message(STATUS "  build type     : ${CMAKE_BUILD_TYPE}")
message(STATUS "  OpenMP         : ${CYFTOOLS_WITH_OPENMP}")
message(STATUS "  TIFF/JPEG      : ${CYFTOOLS_WITH_TIFF}")
message(STATUS "  Cairo          : ${CYFTOOLS_WITH_CAIRO}")
message(STATUS "  LDA            : ${CYFTOOLS_WITH_LDA}")
message(STATUS "  self-test      : ${CYFTOOLS_BUILD_TESTS}")
message(STATUS "")
