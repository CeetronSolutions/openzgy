# Python bindings for the C++ version of OpenZGY

This is an extension enabling Python scripts to read ZGY files using the compiled C++ version of OpenZGY instead of the native Python version.

The main purpose is to compare the two implementations. Both the performance and the test result.

The API for the core part of the library is similar to both the old C++ wrapper around the closed-source ZGY library and to the pure Python implementation of OpenZGY.

Caveat: If you have obtained a pre-compiled wheel for this package, the naming convention of wheels means you can immediately see if it targets the wrong version of Python. But you will not see if the gcc version matches.

There is currently no Windows support. Adding this is doable but will further increase the maintenance cost.
