cmake_minimum_required(VERSION 3.22)
project(ESMA VERSION 0.1 LANGUAGES CXX)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
message(STATUS "CXX is ${CXX}")
option(ESMA_ENABLE_TRAJECTORIES "Record and write ambulance trajectories" ON)
option(ESMA_TRAVEL_STREET "Use StreetTravel (OSRM) as Travel type" OFF)
option(ESMA_TRAVEL_CARTESIAN "Use CartesianTravel as Travel type" OFF)

# Dependências
find_package(Boost 1.74 REQUIRED COMPONENTS unit_test_framework program_options filesystem thread system iostreams chrono date_time regex)
find_package(GUROBI REQUIRED)
find_package(fmt REQUIRED)
find_package(spdlog REQUIRED)
find_package(cpptrace REQUIRED)
# OSRM: use static import at /usr/local/lib/libosrm.a if present, otherwise try
# a config package. Required for StreetTravel.
if(NOT osrm_DIR AND EXISTS "/usr/local/lib/libosrm.a")
  add_library(osrm STATIC IMPORTED)
  set_target_properties(osrm PROPERTIES
    IMPORTED_LOCATION "/usr/local/lib/libosrm.a"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include;/usr/local/include/osrm")
  set(osrm_FOUND TRUE)
else()
  find_package(osrm CONFIG REQUIRED)
endif()

set(ESMA_ENABLE_GPERF_PROFILER OFF)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(ESMA_ENABLE_GPERF_PROFILER ON)
endif()

if(ESMA_ENABLE_GPERF_PROFILER)
  find_library(GPERFTOOLS_PROFILER_LIB NAMES profiler)
  if(GPERFTOOLS_PROFILER_LIB)
    message(STATUS "Found gperftools profiler: ${GPERFTOOLS_PROFILER_LIB}")
  else()
    message(WARNING "gperftools profiler library (libprofiler) not found; CPU profiling will be disabled")
  endif()
endif()

# Include dirs
include_directories(${GUROBI_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(src)

message("Gurobi libs " ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY})

# Travel library
add_library(travel src/esma/travel.cpp)
target_compile_options(travel PRIVATE -Wall -Wextra -Werror -fmax-errors=1)
target_link_libraries(travel PRIVATE ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} osrm Boost::filesystem Boost::thread Boost::system Boost::iostreams Boost::chrono Boost::date_time Boost::regex pthread dl m)
if(ESMA_TRAVEL_STREET AND ESMA_TRAVEL_CARTESIAN)
  message(FATAL_ERROR "Only one of ESMA_TRAVEL_STREET or ESMA_TRAVEL_CARTESIAN can be ON.")
endif()
if(ESMA_TRAVEL_STREET)
  target_compile_definitions(travel PRIVATE ESMA_TRAVEL_STREET=1)
elseif(ESMA_TRAVEL_CARTESIAN)
  target_compile_definitions(travel PRIVATE ESMA_TRAVEL_CARTESIAN=1)
endif()

# Testes
add_executable(test_travel test/test_travel.cpp)
target_link_libraries(test_travel PRIVATE travel fmt::fmt cpptrace::cpptrace spdlog::spdlog ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} osrm Boost::filesystem Boost::thread Boost::system Boost::iostreams Boost::chrono Boost::date_time Boost::regex pthread dl m)

# Se estiver no MinGW (Windows), precisa do ws2_32
target_link_libraries(test_travel PRIVATE $<$<BOOL:${MINGW}>:ws2_32>)

target_compile_options(test_travel PRIVATE -Wall -Wextra -Werror -fmax-errors=1)
if(ESMA_ENABLE_TRAJECTORIES)
  target_compile_definitions(test_travel PRIVATE ESMA_ENABLE_TRAJECTORIES=1)
else()
  target_compile_definitions(test_travel PRIVATE ESMA_ENABLE_TRAJECTORIES=0)
endif()
if(ESMA_TRAVEL_STREET)
  target_compile_definitions(test_travel PRIVATE ESMA_TRAVEL_STREET=1)
elseif(ESMA_TRAVEL_CARTESIAN)
  target_compile_definitions(test_travel PRIVATE ESMA_TRAVEL_CARTESIAN=1)
endif()

# Programa principal
add_executable(esma src/main.cpp src/esma/param.cpp)
target_link_libraries(esma PRIVATE travel fmt::fmt cpptrace::cpptrace spdlog::spdlog ${Boost_LIBRARIES} ${GUROBI_LIBRARY} ${GUROBI_CXX_LIBRARY} osrm Boost::filesystem Boost::thread Boost::system Boost::iostreams Boost::chrono Boost::date_time Boost::regex pthread dl m)
if(ESMA_ENABLE_GPERF_PROFILER AND GPERFTOOLS_PROFILER_LIB)
  target_link_libraries(esma PRIVATE ${GPERFTOOLS_PROFILER_LIB})
endif()
target_link_libraries(esma PRIVATE $<$<BOOL:${MINGW}>:ws2_32>)

target_compile_options(esma PRIVATE -Wall -Wextra -Werror -fmax-errors=1)
if(ESMA_ENABLE_TRAJECTORIES)
  target_compile_definitions(esma PRIVATE ESMA_ENABLE_TRAJECTORIES=1)
else()
  target_compile_definitions(esma PRIVATE ESMA_ENABLE_TRAJECTORIES=0)
endif()
if(ESMA_TRAVEL_STREET)
  target_compile_definitions(esma PRIVATE ESMA_TRAVEL_STREET=1)
elseif(ESMA_TRAVEL_CARTESIAN)
  target_compile_definitions(esma PRIVATE ESMA_TRAVEL_CARTESIAN=1)
endif()

# Testes automatizados
enable_testing()
add_test(NAME ESMA_TEST COMMAND test_travel)

add_custom_target(run_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    DEPENDS test_travel
)

add_dependencies(run_tests test_travel)
