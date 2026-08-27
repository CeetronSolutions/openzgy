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

# Will export the following variables:
#    HOSTDISTRO LINUXDISTRO (os name), BUILDROOT TESTRUNDIR (paths)
#    PLATFORM CONFIG (consumed) MYPLATFORM MYCONFIG (produced)
SELF  := $(lastword $(MAKEFILE_LIST))
include ../scripts/Make.defs

MYNAME  = OpenZgyBindings
MYNAMEI = [Oo]pen[Zz]gy[Bb]indings
WHEELDIR= $(if $(OPENZGY_VANILLA),without-dependencies,$(MYPLATFORM))
WHEELEND= $(if $(OPENZGY_VANILLA),$(MYNAMEI)*.tar.gz,$(MYNAMEI)*.whl)
MYWHEEL = $(BUILDROOT)/build/deploy/wrapper/$(WHEELDIR)/$(WHEELEND)
MYSTAMP = $(BUILDROOT)/build/deploy/wrapper/$(WHEELDIR)/$(MYNAME).$(PYTHON).timestamp
DSO     = build/lib.linux-x86_64-*/openzgycpp.cpython-*.so
BIN_DIR = $(BUILDROOT)/build/deploy/native/$(MYPLATFORM)

INPUTS  =
INPUTS += Makefile setup.py requirements.txt
INPUTS += wrappermodule.cpp
INPUTS += $(BIN_DIR)/libopenzgy.so

# This Makefile does not support parallel execution with -j N.
# There is practically nothing at this level that can be done
# in parallel. For that reason, no effort has been make to ensure
# that the rules are correct as far as parallelization goes.
.NOTPARALLEL:

# This is to avoid polluting the source directory with compiled files.
# It is somewhat pointless if bdist_wheel still places intermediate
# files in source. Working on that.
export PYTHONDONTWRITEBYTECODE=1

#
# Here are the garbage files.
#
#   build/bdist.linux-x86_64/
#   build/lib.linux-x86_64-3.6/
#   build/temp.linux-x86_64-3.6/
#   OpenZgyBindings.egg-info/
#   dist/
#
# Suggested contents of setup.cfg
# This MIGHT work but needs more testing. Especially because
# tweaking that big mess which is setuptools might be a bad idea.
#
# [build]
# build_base = ../build/temp/wrapper/BUILDBASE
# debug = false
#
# [build_ext]
# debug = false
#
# [bdist_wheel]
# bdist_dir = ../build/temp/wrapper/BDIST
# dist_dir = ../build/temp/wrapper/DIST
# keep_temp = true
#
# [egg_info]
# egg_base = ../build/temp/wrapper
#

ifeq ($(strip $(LINUXDISTRO)),clang)
CC           = clang
CXX          = clang++
endif

all: generic specific

all_build: generic_build specific_build

# Note that OPENZGY_SDK is only set while running tests, if at all.
# This means that when installing generic/wheel this will always
# need to have the SDK installed in /usr/local/openzgy.
# The tests in the "generic" target will only test generic/sdist.

ifeq ($(strip $(PYTHONLIST)),)
PYTHONLIST := $(shell cd /usr/include; ls -d python3*)
endif
PYTHON ?= /bin/false

$(info Python interpreters: $(PYTHONLIST); current is $(PYTHON))

generic:
	@echo -e "\nBuilding 'generic' assuming OpenZGY SDK is present.\n"
	for PYTHON in $(PYTHONLIST); do \
		echo $@ PYTHON=$${PYTHON}; \
		env OPENZGY_VANILLA=1 OPENZGY_SDK='' $(MAKE) PYTHON=$${PYTHON} clean build-one; \
		env OPENZGY_VANILLA=1 OPENZGY_SDK=$(BUILDROOT)/build/deploy $(MAKE) PYTHON=$${PYTHON} test; \
	done

generic_build:
	@echo -e "\nBuilding 'generic' assuming OpenZGY SDK is present.\n"
	for PYTHON in $(PYTHONLIST); do \
		echo $@ PYTHON=$${PYTHON}; \
		env OPENZGY_VANILLA=1 OPENZGY_SDK='' $(MAKE) PYTHON=$${PYTHON} clean build-one; \
	done

specific:
	@echo -e "\nBuilding '$(MYPLATFORM)'\n"
	for PYTHON in $(PYTHONLIST); do \
		echo $@ PYTHON=$${PYTHON}; \
		env OPENZGY_VANILLA= $(MAKE) PYTHON=$${PYTHON} clean build-one test; \
	done

specific_build:
	@echo -e "\nBuilding '$(MYPLATFORM)'\n"
	for PYTHON in $(PYTHONLIST); do \
		echo $@ PYTHON=$${PYTHON}; \
		env OPENZGY_VANILLA= $(MAKE) PYTHON=$${PYTHON} clean build-one; \
	done

# Only for HAVE_SD and HAVE_ZFP variables.
include  ../native/src/Make.packages

# Caveat: Source tree needs to be partly writable to hold temp files.
# Final deliverable is ok; it can be specified by --dist-dir.
# bdist_wheel will create writable ./build and ./OpenZgyBindings.egg-info.
# The bdist_wheel --bdist-dir option is not much help; it defaults
# to ./build/bdist.linux-x86_64, and ./build gets created anyway.

# The line for building the wrapper in debug mode is commented out
# because, unlike Windows, mixing debug and release is not a problem.

build-one: $(MYSTAMP)

$(MYSTAMP): $(INPUTS)
	@echo "Cleaning to force a full rebuild."
	$(MAKE) -f $(SELF) clean
	$(MAKE) -f $(SELF) v.init
	/bin/mkdir -p $(dir $(MYSTAMP)) $(dir $(MYWHEEL))
	/bin/mkdir openzgy-sdk
	/bin/ln -s $(BIN_DIR) openzgy-sdk/lib
	/bin/ln -s $(BUILDROOT)/native/src openzgy-sdk/include
	/bin/rm -f setup.cfg
	@echo "Start the full rebuild now."
	@: Send the next line to setup.cfg if a debug build is desired.
	echo -e '[build]\ndebug = 1\n[build_ext]\ndebug = 1' > /dev/null
	(. $(ENVDIR)/bin/activate; python3 setup.py --no-user-cfg sdist bdist_wheel)
	/bin/cp -t $(dir $(MYWHEEL)) requirements.txt
	/bin/cp -t $(dir $(MYWHEEL)) dist/*
	touch $(MYSTAMP)
	@echo "Build of $@ is complete."

.PHONY: all all_build generic generic_build specific specific_build build test

#############################################################################
###   Virtualenv   ##########################################################
#############################################################################

ENVDIR = test-venv
PIP    = python3 -m pip --disable-pip-version-check -q
PIP    = pip3 --disable-pip-version-check -q

.PHONY: venv init install forceinstall uninstall

venv: $(ENVDIR)/bin/activate
	@/bin/true

# Create the $(ENVDIR) virtualenv if not done yet.
# Install third party packages but not the packages we are building.
# Note that if the virtualenv exists it is assumed to be up to date.
$(ENVDIR)/bin/activate: requirements.txt
	$(RM) -rf $(ENVDIR)
	$(PYTHON) -m venv $(ENVDIR)
	$(MAKE) -f $(SELF) v.init

init: requirements.txt
	$(if $(VIRTUAL_ENV),,$(error $@ must be run in a virtualenv))
	$(PIP) install  --upgrade pip
	$(PIP) install  --upgrade setuptools wheel
	$(PIP) install -r requirements.txt

install:
	$(if $(VIRTUAL_ENV),,$(error $@ must be run in a virtualenv))
	-$(PIP) uninstall -y $(MYNAME)
	$(PIP) install $(MYWHEEL)

# Install while skipping the virtualenv check. Useful in docker images.
forceinstall:
	-$(PIP) uninstall -y $(MYNAME)
	$(PIP) install $(MYWHEEL)

uninstall:
	-pip3 uninstall $(MYNAME)

# Any target v.FOO means "run FOO inside the virtualenv".
# Note that v.FOO shouldn't be declared phony, even if it arguably is.
v.%: $(ENVDIR)/bin/activate
	@echo "Building $(@:v.%=%) inside virtualenv"
	@(. $(ENVDIR)/bin/activate; $(MAKE) -f $(SELF) $(@:v.%=%))

EXAMPLE:
	$(if $(VIRTUAL_ENV),,$(error $@ must be run in a virtualenv))

#############################################################################
#############################################################################
#############################################################################

# If I want builds to run outside the virtualenv, declare it like this
# instead of having install and test depending on build.
# Use "$MYSTAMP" not "build" because the latter target does
# some work even if fully up to date.
v.forceinstall v.install v.realtest: $(MYSTAMP)

SAUTHENV=$(shell test ! -x ../../private/refreshtoken.sh || ../../private/refreshtoken.sh)
test: v.install
	@$(MAKE) $(TESTARG) $(SAUTHENV) v.realtest

verbose: $(MYSTAMP)
	unzip -l $(MYWHEEL)
	readelf -d $(DSO) | egrep 'NEEDED|RUNPATH'
	ldd $(DSO)

realtest:
	(cd test; python3 test.py)
	(cd test; python3 test_pickle.py)
	(cd test; python3 test_black.py local $(if $(strip $(HAVE_SD)),sd,,))

clean:
	$(RM) -rf MANIFEST OpenZgyBindings-1.0.dist-info OpenZgyBindings.egg-info build unittest1.zgy unittest2.zgy unittest3.zgy __pycache__ openzgycpp/__pycache__ test/__pycache__
	$(RM) -rf test-venv
	$(RM) -rf dist
	$(RM) -f openzgy-sdk/lib openzgy-sdk/include setup.cfg
	test ! -e openzgy-sdk || /bin/rmdir openzgy-sdk

clobber: clean
	$(RM) -rf $(BUILDROOT)/build/deploy/wrapper/$(MYPLATFORM)
	$(RM) -rf $(BUILDROOT)/build/deploy/wrapper/without-dependencies
	$(RM)     $(BUILDROOT)/build/deploy/wrapper/requirements.txt

.PHONY: verbose test realtest pushtest clean clobber
