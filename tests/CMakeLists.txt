
# OpenZGY wrapper tests

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.11.0
)

set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
set(BUILD_GTEST ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE )

FetchContent_MakeAvailable(googletest)

set(THREADS_PREFER_PTHREAD_FLAG ON)

# Add executables first
add_executable(openzgy-tests reader_tests.cpp geometry_tests.cpp main.cpp)

# location of test data
CONFIGURE_FILE( ${CMAKE_CURRENT_LIST_DIR}/testdatafolder.h.cmake
                ${CMAKE_BINARY_DIR}/Generated/testdatafolder.h
)

target_include_directories(openzgy-tests PUBLIC ../include/ ${CMAKE_BINARY_DIR}/Generated)

if(MSVC)
  target_compile_options(openzgy-tests PRIVATE /W4 /WX)
else()
  target_compile_options(openzgy-tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()


target_compile_features(openzgy-tests PRIVATE cxx_std_20)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if (MSVC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
endif()

target_link_libraries(openzgy-tests PUBLIC openzgy gtest gtest_main )

add_test(NAME openzgy-tests COMMAND openzgy-tests)