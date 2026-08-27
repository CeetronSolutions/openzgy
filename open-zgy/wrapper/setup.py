#!/usr/bin/env python

"""
Build and install OpenZGY bindings for Python.

Usage (Linux):
    /bin/ln -s $(BUILDROOT)/build/deploy/native/$(MYPLATFORM) openzgy-sdk/lib
    /bin/ln -s $(BUILDROOT)/native/src openzgy-sdk/include
    python3 setup.py sdist bdist_wheel

Usage if all-in-one is used (Windows):
    /bin/ln -s $(BUILDROOT)/build/deploy/native/$(MYPLATFORM) openzgy-sdk
    /bin/ln -s $(BUILDROOT)/native/src openzgy-sdk/include
    /bin/cp -t openzgy-sdk openzgycpp/__init__.py
"""

from setuptools import setup, Extension
import os
import sys
from glob import glob

########################################################################
# The bindings are provided as a compiled Python extension.
# setup.py only compiles the wrapper itself. The binary OpenZGY
# SDK is treated as data.
#
# Building a "bdist_wheel" results in a wheels specific to both
# a Python version and a Linux distro or Windows compiler.
#
# Creting an "sdist" is possible and will make the installation
# work with any Python version >= 3.6. But this will not be a
# true source distribution because of the bundled binary sdk.
# So the "sdist" is still specific to distro and compiler.
#
# An altenative would be to have the OpenZGY SDK installed by
# other means (dnf, yum, apt, msi, etc. or by extracting a tar)
# before pip install of the wrapper. This approach is not yet
# implented but would be fairly simple to do.
#
# An altenative would be to have the "sdist" build the entire
# OpenZGY. This is not feasible if Seismic Store is enabled,
# but might barely be doable for the core OpenZGY. Probably not,
# though. ZFP might be a hassle both technically and legally.
#
# Look in the git history, in particular the entry before
# "Choose the setup strategy to use", for experiments with
# other ways of setting this up.
########################################################################

version   = "0.2." + (os.getenv("AZURE_BUILDID", "dev0") or "dev0")
default_sdk = "C:\\OpenZGY" if os.name == 'nt' else "/usr/local/openzgy"
openzgy_sdk = os.getenv('OPENZGY_SDK', default_sdk) or default_sdk
linuxdistro = os.getenv('LINUXDISTRO', 'generic-linux')
print("Wrapper", version, "setup for", linuxdistro, "with SDK", openzgy_sdk, file=sys.stderr)

def platform():
    """
    The native SDK uses the platform name in paths etc.
    Here is a minor kludge to get the correct name.
    Note that os.path.islink() doesn't work on windows
    (wrap os.readlink in a try/catch instead) but
    currently there is just a single platform there.

    WARNING: Using this function will cause setup.py
    to behave differently depending on where it is
    run from. This is very much frowned upon.
    """
    if os.name == 'nt':
        return path.join('x64', 'Release')
    elif os.path.islink('openzgy-sdk/lib'):
        # Probably building sdist or wheel on Linux
        return os.path.basename(os.readlink('openzgy-sdk/lib'))
    elif os.path.islink('openzgy-sdk'):
        # Probably building sdist or wheel with Windows strategy
        return os.path.basename(os.readlink('openzgy-sdk'))
    else:
        installed_lib = openzgy_sdk + '/native/*-gcc*'
        globbed_lib = glob(installed_lib)
        if len(globbed_lib) == 1:
            # Probably a user installing or converting an sdist
            return os.path.basename(globbed_lib[0])
        elif len(globbed_lib) > 1:
            print("*** ERROR: Found multiple", installed_lib, file=sys.stderr)
            print("*** Please remove those you don't need.", file=sys.stderr)
            return 'unknown'
        else:
            print("*** ERROR: Cannot find", installed_lib, file=sys.stderr)
            print("*** Did you forget to set $OPENZGY_SDK ?", file=sys.stderr)
            return 'unknown'

def vanilla():
    """
    Decide whether to build or install a self contained package
    (regular) or one that expects the native OpenZGY sdk to
    already be installed on the target (vanilla).

    WARNING: Using this function will cause setup.py
    to behave differently depending on where it is
    run from. This is very much frowned upon.
    """
    if not os.path.exists('openzgy-sdk'):
        # Installing or converting a vanilla sdist.
        return True
    if os.getenv("OPENZGY_VANILLA", ""):
        # Used to force building vanilla sdist and/or wheel.
        return True
    return False

if os.name == 'nt':
    installed_root = openzgy_sdk.replace("\\", "/")
    installed_inc = installed_root + '/native/include/openzgy'
    installed_lib = installed_root + '/native/x64/Release'
    module1 = Extension(
        'openzgycpp.wrapper',
        libraries    = ['openzgy' ],
        sources      = ['wrappermodule.cpp'],
        include_dirs = ['openzgy-sdk/include', installed_inc],
        library_dirs = ['openzgy-sdk', installed_lib],
        # MAJOR kludge, %LIB% or /LIBPATH missing Python.dll ?
        #library_dirs = ['openzgy-sdk', installed_lib, r'c:\Program Files\Python38\libs'],
    )

    # Everything in a single package. Windows has no runtime_library_dirs.
    # This is a bit messy because the source tree needs to be modified so
    # it looks like the actual source (__init__.py) and the prebuilt binaries
    # are in the same folder. Handle this copying outside setup.py.
    data_args = dict(
        packages     = ['openzgycpp'],
        package_dir  = {'openzgycpp': 'openzgy-sdk'},
        package_data = {'openzgycpp': ["*.dll", "*.lib", "include/*.h"]},
    )
    # Or, assume native OpenZGY already installed in a well known location.
    if vanilla():
        data_args =  dict(packages = ['openzgycpp'])

else:
    installed_inc = openzgy_sdk + '/native/include/openzgy'
    installed_lib = openzgy_sdk + '/native/' + platform()
    module1 = Extension(
        'openzgycpp.wrapper',
        libraries    = ['openzgy' ],
        sources      = ['wrappermodule.cpp'],
        include_dirs = ['openzgy-sdk/include', installed_inc],
        library_dirs = ['openzgy-sdk/lib', 'openzgy-sdk', installed_lib],
        # This does add confusion in the vanilla a.k.a. generic case
        # because the wheel captures openzgy_sdk when built while the
        # sdist captures it on install. In the non-vanilla wheel and sdist
        # the .so is found relative to the origin in, so this issue is N/A.
        runtime_library_dirs = ['$ORIGIN/../openzgy-sdk/lib', installed_lib],
        #runtime_library_dirs = ['$ORIGIN/../openzgy-sdk/lib'],
        extra_compile_args   = ['-fopenmp', '-std=c++11', '-O2'] + (['-D_GLIBCXX_USE_CXX11_ABI=0'] if linuxdistro == 'omega8' else []),
        extra_link_args      = ['-fopenmp', '-Wl,-z,origin'])
    # Two Python packages. One for the Python part plus the
    # wrapper binary and one for the wrapper's dependencies.
    data_args =  dict(
        packages     = ['openzgycpp', 'openzgy-sdk'],
        package_dir  = {'openzgy-sdk': 'openzgy-sdk'},
        package_data = {'openzgy-sdk': ["lib/[!.]*.so*", "include/*.h"]})
    # Or, assume native OpenZGY already installed in a well known location.
    if vanilla():
        data_args =  dict(packages = ['openzgycpp'])

setup (name = 'OpenZgyBindings',
       license='Apache',
       classifiers=[
           "Programming Language :: C++ :: 11",
           "License :: OSI Approved :: Apache Software License",
       ],
       version      = version,
       description  = 'OpenZGY bindings for Python',
       author       = 'Schlumberger',
       #author_email = 'support@slb.com',
       zip_safe     = False,
       ext_modules  = [module1],
       **data_args)

# There is a warning about mixing debug and release on Windows.
# Which is odd, because a manual check says there isn't a problem.
import sys, sysconfig
print('Check for mixed Debug and Release:',
      'PyDEBUG', sysconfig.get_config_var('Py_DEBUG'),
      'gettotalrefcount', hasattr(sys, 'gettotalrefcount'),
      'interpreter', sys.implementation.name,
      'executable', sys.executable)

#
# See https://packaging.python.org/extensions
# See https://docs.python.org/3.4/distutils/index.html
# https://packaging.python.org/guides/tool-recommendations/
# https://setuptools.readthedocs.io/en/latest/

#
# Build and test instructions Windows.
# Do NOT run the tests with CWD=wrapper.
#    <start Visual Studio command shell>
#    cd .../wrapper
#    rd /s/q dist build test-venv
#    python setup.py bdist_wheel
#    virtualenv test-venv
#    test-venv\Scripts\activate.bat
#    pip install dist\OpenZgyBindings-0.2.dev0-cp38-cp38-win_amd64.whl
#    cd test
#    python test_black.py local
#

# Copyright 2017-2022, Schlumberger
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
