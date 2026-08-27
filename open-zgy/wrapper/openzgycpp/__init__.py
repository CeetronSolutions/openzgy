import os
from contextlib import contextmanager
from glob import glob

# Helpers to load dependencies of the compiled extension.
# This is WAY too complicated.

@contextmanager
def _temporary_setenv(updates):
    """
    In a "with" context, temporarily set an environment variable
    that will automatically be reset when the "with" goes out of scope.
    Why isn't there a standard function for doing this? Or is there?
    """
    saved = dict()
    for key, value in updates.items():
        if value != os.environ.get(key):
            saved[key] = os.environ.get(key)
            if value is None:
                print("Temporarily delete", key)
                del os.environ[key]
            else:
                print("Temporarily set", key, "=", value)
                os.environ[key] = value
    try:
        yield
    finally:
        for key, value in saved.items():
            if value is None:
                print("Restore (delete)", key)
                del os.environ[key]
            else:
                print("Restore", key, "=", value)
                os.environ[key] = value

@contextmanager
def _temporary_addenv(key, value, sep, prepend):
    """
    In a "with" context, temporarily modify an environment variable
    that contains colon- or semicolon separated values. The variable
    will automatically be reset when the "with" goes out of scope.
    """
    old = os.environ.get(key)
    if not old:
        new = value
    elif (sep + value + sep) not in (sep + old + sep):
        print('"{0}" is not in "{1}"'.format((sep + value + sep), (sep + old + sep)))
        new = value + sep + old if prepend else old + sep + value
    else:
        new = old
    with _temporary_setenv({key: new}):
        yield

@contextmanager
def _temporary_searchpath(folder):
    """
    In a "with" context, temporarily set a folder that will be searched
    to find dependencies of the compiled extension. Or any other dll for
    that matter. The setting will automatically be reset when the "with"
    goes out of scope.

    In Python for Windows 3.8 and later, searching for dependencies
    needed by a compiled extension only looks in the same directory
    as the extension and in a few system folders. The normal Windows
    search rules are disabled. Older versions of Python will look
    among other places in %PATH%.
    """
    if not folder:
        print("import openzgy.wrapper with default search path")
        yield
    elif os.name == 'nt':
        if hasattr(os, 'add_dll_directory'):
            with os.add_dll_directory(folder):
                print("import openzgy.wrapper with", folder, "in search path")
                yield
        else:
            with _temporary_addenv("PATH", folder, sep=';', prepend=True):
                print("import openzgy.wrapper with", folder, "in %PATH%")
                yield
    else:
        print("import openzgy.wrapper with", folder, "in $LD_LIBRARY_PATH")
        with _temporary_addenv("LD_LIBRARY_PATH", folder, sep=':', prepend=True):
            yield

def _find_openzgy_library_nt():
    libname = 'OpenZGY.dll'
    search  = [os.getenv("OPENZGY_SDK", None),
               os.path.join(os.getenv("ProgramFiles", r"C:\Program Files"), "OpenZGY"),
               r"C:\OpenZGY"]
    relpath = r"native\x64\Release"
    lib = os.path.join(os.path.dirname(__file__), libname)
    lib = lib.replace('/', '\\')
    if os.path.exists(lib):
        print("Found", lib, "bundled with openzgycpp")
        return None
    for folder in search:
        if folder:
            lib = os.path.join(folder, relpath, libname)
            lib = lib.replace('/', '\\')
            if os.path.exists(lib):
                print("Was found:", lib)
                return os.path.dirname(lib)
            print("Not found:", lib)
    print("Cannot find", libname, "anywhere.")
    print("If import of openzgycpp.wrapper fails, this is probably why.")
    return None

def _find_openzgy_library_linux():
    libname = "libopenzgy.so*"
    search  = [os.getenv("OPENZGY_SDK", "/usr/local/openzgy")]
    relpath = "native/*-gcc*"
    lib = os.path.join(os.path.dirname(__file__),
                       "..", "openzgy-sdk", "lib", libname)
    if len(glob(lib)) > 0:
        print("Found", lib, "bundled with openzgycpp")
        return None
    for folder in search:
        if folder:
            libdir = glob(os.path.join(folder, relpath))
            if len(libdir) == 1:
                print("Was found:", libdir[0])
                return libdir[0]
            elif len(libdir) > 1:
                # Problems reliably deciding what the current architecture is.
                raise RuntimeError("More than one OpenZGY SDK installed. Don't know which to pick:", *libdir)
            print("Not found in:", os.path.join(folder, relpath))
    print("Cannot find", libname, "anywhere.")
    print("If import of openzgycpp.wrapper fails, this is probably why.")
    return None

def _find_openzgy_library():
    if os.name == 'nt':
        return _find_openzgy_library_nt()
    else:
        return _find_openzgy_library_linux()

"""
Make all wrapper methods visible at the top level.

Handle the case where the SDK is installed separately.
If bundled with the Python code the SDK dlls will be found relative to,
or in, the folder where the compiled extension is. _find_openzgy_library
will return None in that case which makes _temporary_searchpath a no-op.

On Windows the separate SDK works. At least with Python 3.8 and later.
Because changing the search path is allowed while the process is running.
On Linux the value of LD_LIBRARY_PATH is read only at startup. So this
code won't work out of the box. There are several workarounds. All smelly.

 - When installing an sdist on linux with a separate SDK, the SDK path can
   be captured at install time and stored as the runpath of the compiled
   wrapper. Installing a wheel that way is not possible because there is no
   post-install hook. In that case the SDK needs to be installed in, or linked
   from, the default location /usr/local/openzgy. Or maybe just don't provide
   any wheels for the separate SDK case.

 - Or, re-exec the application from inside __init__.py to pick up the new
   search path. Which then remains in force until the application exits.
   The latter is no big issue but the former stinks.

-  Or, require the user to either set LD_LIBRARY_PATH before starting the
   application, or place the SDK where ldconfig will find it.

-  Or, put back the ctypes patch from the git history. ctypes can preload
   a shared object using an absolute path. This kludge would be tolerable
   on Windows but gets uglier on Linux because the .so files are typically
   known by more that one name.
"""
with _temporary_searchpath(_find_openzgy_library()):
    from .wrapper import *

# Might skip this so the code can be better unit tested.
del _find_openzgy_library
del _find_openzgy_library_linux
del _find_openzgy_library_nt
del _temporary_searchpath
del _temporary_addenv
del _temporary_setenv
