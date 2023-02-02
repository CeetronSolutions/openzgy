cmake_minimum_required(VERSION 3.5)

include(src/openzgy/impl/CMakeLists.txt)

set(HEADER_FILES ${HEADER_FILES}
    include/openzgy/exception.h
    src/openzgy/iocontext.h
    src/openzgy/safewriter.h
)

set(SOURCE_FILES ${SOURCE_FILES}
    src/openzgy/api.cpp
    src/openzgy/exception.cpp
    src/openzgy/iocontext.cpp
)
