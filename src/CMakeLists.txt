cmake_minimum_required(VERSION 3.5)

include(impl/CMakeLists.txt)

set(HEADER_FILES ${HEADER_FILES}
    api.h
    declspec.h
    exception.h
    iocontext.h
    safewriter.h
)

set(SOURCE_FILES ${SOURCE_FILES}
    api.cpp
    exception.cpp
    iocontext.cpp
)
