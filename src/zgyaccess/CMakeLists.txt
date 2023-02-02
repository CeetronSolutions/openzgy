cmake_minimum_required(VERSION 3.5)

set(HEADER_FILES ${HEADER_FILES}
	include/zgyaccess/zgyreader.h
	include/zgyaccess/seismicslice.h
)

set(SOURCE_FILES ${SOURCE_FILES}
	src/zgyaccess/zgyreader.cpp
	src/zgyaccess/seismicslice.cpp
)
