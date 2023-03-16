# openzgy
A wrapper around the native C++ part of the OpenZGY library from https://community.opengroup.org/osdu/platform/domain-data-mgmt-services/seismic/open-zgy

- Allows you to open the repo directly in Visual Studio 2022 using the built-in CMake support.
- Allows building the native C++ library using CMake
- Uses the ZFP library as a GIT submodule
- Uses the Open ZGY library as a GIT submodule
- Supports both Linux and Windows builds, using C++ 20
- Supports both clang and gcc on Linux (no tests built with clang)
- Allows you to build a linux version using Visual Studio 2022 and WSL2 support

Provides a simplified C++ API that is used by the ResInsight software ( https://resinsight.org/ ) to access seismic data:
- Open/close ZGY files on local or network disk
- Access file meta information and data histogram
- Read inline/crossline/z slices
- Read individual z traces

The API is provided by the ZGYAccess::ZGYReader class found in include/zgyreader.h

