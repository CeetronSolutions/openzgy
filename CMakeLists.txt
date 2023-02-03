cmake_minimum_required(VERSION 3.15)

project(openzgy)

add_definitions(-DOPENZGY_STATIC -DHAVE_ZFP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-parameter -Wno-reorder -Wno-parentheses -Wno-switch -Werror"  )

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-delete-abstract-non-virtual-dtor -Wno-undefined-var-template -Wno-invalid-source-encoding -Wno-enum-compare -Wno-call-to-pure-virtual-from-ctor-dtor -Wno-unused-variable -Wno-unused-private-field -Wno-unused-lambda-capture -Wno-delete-non-abstract-non-virtual-dtor -Wno-braced-scalar-init -Wno-tautological-constant-out-of-range-compare -Wmissing-braces")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-restrict"  )
  endif()
endif()

if (MSVC)
    # Disable some annoying warnings (relative to warning level 3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4251")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  message(STATUS "Enabling OpenMP support")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OPENMP")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}" )
endif()

# Compiles static library
include(src/CMakeLists.txt)

add_library(${PROJECT_NAME} OBJECT ${HEADER_FILES} ${SOURCE_FILES} )

target_include_directories(${PROJECT_NAME}
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/zfp/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include/openzgy
  ${CMAKE_CURRENT_SOURCE_DIR}/include/zgyaccess
  ${CMAKE_CURRENT_SOURCE_DIR}/src/openzgy
)

set(BUILD_TESTING "OFF" CACHE BOOL "No testing enabled" FORCE)
add_subdirectory(zfp)

target_link_libraries(${PROJECT_NAME} zfp)
