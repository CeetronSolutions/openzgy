cmake_minimum_required(VERSION 3.5)

include(src/impl/CMakeLists.txt)

set(HEADER_FILES ${HEADER_FILES}
    src/api.h
    src/declspec.h
    src/exception.h
    src/iocontext.h
    src/safewriter.h
)

set(SOURCE_FILES ${SOURCE_FILES}
    src/api.cpp
    src/exception.cpp
    src/iocontext.cpp
)
