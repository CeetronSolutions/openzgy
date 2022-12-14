cmake_minimum_required(VERSION 3.5)

set(HEADER_FILES ${HEADER_FILES}
    src/impl/arrayops.h
    src/impl/bulk.h
    src/impl/compression.h
    src/impl/cornerpoints.h
    src/impl/databuffer.h
    src/impl/enum.h
    src/impl/environment.h
    src/impl/fancy_timers.h
    src/impl/file.h
    src/impl/file_consolidate.h
    src/impl/file_parallelizer.h
    src/impl/file_performance.h
    src/impl/file_relay.h
    src/impl/file_smallcache.h
    src/impl/genlod.h
    src/impl/guid.h
    src/impl/histogrambuilder.h
    src/impl/histogramdata.h
    src/impl/iltf2d.h
    src/impl/locker.h
    src/impl/lodalgo.h
    src/impl/lodsampling.h
    src/impl/logger.h
    src/impl/lookuptable.h
    src/impl/meta.h
    src/impl/minmaxscan.h
    src/impl/mtguard.h
    src/impl/perflogger.h
    src/impl/roundandclip.h
    src/impl/sdinterface.h
    src/impl/statisticdata.h
    src/impl/structaccess.h
    src/impl/subtiling.h
    src/impl/timer.h
    src/impl/transform.h
    src/impl/types.h
)

set(SOURCE_FILES ${SOURCE_FILES}
    src/impl/bulk.cpp
    src/impl/compression.cpp
    src/impl/compress_null.cpp
    src/impl/compress_zfp.cpp
    src/impl/cornerpoints.cpp
    src/impl/databuffer.cpp
    src/impl/environment.cpp
    src/impl/fancy_timers.cpp
    src/impl/file.cpp
    src/impl/file_consolidate.cpp
    src/impl/file_local.cpp
    src/impl/file_parallelizer.cpp
    src/impl/file_performance.cpp
    src/impl/file_relay.cpp
    src/impl/file_sd.cpp
    src/impl/file_smallcache.cpp
    src/impl/file_windows.cpp
    src/impl/genlod.cpp
    src/impl/guid.cpp
    src/impl/histogrambuilder.cpp
    src/impl/histogramdata.cpp
    src/impl/iltf2d.cpp
    src/impl/locker.cpp
    src/impl/lodalgo.cpp
    src/impl/lodsampling.cpp
    src/impl/logger.cpp
    src/impl/lookuptable.cpp
    src/impl/meta.cpp
    src/impl/minmaxscan.cpp
    src/impl/mtguard.cpp
    src/impl/perflogger.cpp
    src/impl/sdasyncio.cpp
    src/impl/sdinterface.cpp
    src/impl/statisticdata.cpp
    src/impl/structaccess.cpp
    src/impl/subtiling.cpp
    src/impl/timer.cpp
    src/impl/transform.cpp
    src/impl/types.cpp
)
