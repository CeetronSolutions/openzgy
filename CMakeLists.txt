cmake_minimum_required(VERSION 3.5)

project(openzgy)

# Compiles static library
include(src/CMakeLists.txt)

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES} )

target_include_directories(${PROJECT_NAME}
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# target_link_libraries(${PROJECT_NAME} ${QT_LIBRARIES})
