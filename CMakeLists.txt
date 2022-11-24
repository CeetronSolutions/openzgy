cmake_minimum_required(VERSION 3.5)

project(openzgy)

add_definitions(-DOPENZGY_DLL -DHAVE_ZFP)

if (MSVC)
    # Disable some annoying warnings (relative to warning level 3)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4251")
endif()

# Compiles static library
include(src/CMakeLists.txt)

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES} )

target_include_directories(${PROJECT_NAME}
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/zfp/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

add_subdirectory(zfp)

target_link_libraries(${PROJECT_NAME} zfp)
