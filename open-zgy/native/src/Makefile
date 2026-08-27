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
include ../../scripts/Make.defs

# Uncomment to disable parallel execution with -j N.
# The dependencies might not be 100% correct. My excuse for not fixing
# this properly is that the test target, which just runs a single
# executable, takes much longer than the compliation anyway. The
# automated build pipeline should not use -j. Developers building
# locally can do so at their own risk. Your mileage may vary regarding
# whether the -j option really does speed things up.
#.NOTPARALLEL:

.PHONY: all all-no-doxygen
all: build doxygen doxyint test
all-no-doxygen: build test

# Targets for integrating Seismic Store and ZFP.
include Make.packages

OPTIMIZE     = -O3
CXXFLAGS     =
CXX_ABI      =
ifeq ($(strip $(LINUXDISTRO)),omega8)
CXX_ABI      = -D_GLIBCXX_USE_CXX11_ABI=0
CXXFLAGS    += $(CXX_ABI)
endif
ifeq ($(strip $(LINUXDISTRO)),clang)
CC           = clang
CXX          = clang++
CXXFLAGS    += -std=c++20 -fopenmp=libomp
else
CXXFLAGS    += -std=c++11 -fopenmp
CXXFLAGS    += -Werror=maybe-uninitialized
endif
CXXFLAGS    += -Wall -Wno-deprecated -Wno-unknown-pragmas
CXXFLAGS    += -Werror=return-type
CXXFLAGS    += -Werror=uninitialized
CXXFLAGS    += -fPIC -DLINUX -g
CXXFLAGS    += $(OPTIMIZE)
#CXXFLAGS    += -fno-strict-aliasing -fno-tree-vectorize

ORIGIN       = -Wl,-z,origin,-rpath,\$$ORIGIN
ORIGIN      += -Wl,-rpath-link,$(BIN_DIR) # binutils < 2.28
TESTDATA     = ../../testdata
TESTDATARW   = $(BUILDROOT)/build/testdata

LIB =
LIB += $(ORIGIN)
ifeq ($(strip $(LINUXDISTRO)),clang)
LIB += -fuse-ld=lld -latomic
endif

CXXFLAGS += $(SD_CXXFLAGS) $(ZFP_CXXFLAGS)

NATIVE  = $(if $(COVERAGE),coverage,native)
OBJ_DIR = $(BUILDROOT)/build/temp/$(NATIVE)
BIN_DIR = $(BUILDROOT)/build/deploy/$(NATIVE)/$(MYPLATFORM)
INC_DIR = $(BUILDROOT)/build/deploy/$(NATIVE)/include/openzgy
DOX_DIR = $(BUILDROOT)/build/deploy/$(NATIVE)
LOG_DIR = $(BUILDROOT)/build/deploy/$(NATIVE)

# Input:       COVERAGE {MAIN,TEST,TOOL}_SRC
# Requires:    OPTIMIZE=
# Modifies:    CXXFLAGS LDFLAGS
# Adds target: "coverage"
include Make.coverage

LIBDSO_SYMLINK = $(BIN_DIR)/libopenzgy.so
LIBDSO_LINK2   = $(LIBDSO_SYMLINK).0.2
LIBDSO         = $(LIBDSO_SYMLINK).0.2.$(or $(AZURE_BUILDID),unknown)

PUBLISHED_INC = $(sort $(wildcard *.h))

MAIN_SRC = $(sort $(wildcard *.cpp impl/*.cpp capi/*.cpp))
TEST_SRC = $(sort $(wildcard test/*.cpp))
TOOL_SRC = $(sort $(wildcard tools/*.cpp))

MAIN_OBJ = $(addprefix $(OBJ_DIR)/, $(MAIN_SRC:.cpp=.o))
TEST_OBJ = $(addprefix $(OBJ_DIR)/, $(TEST_SRC:.cpp=.o))
TOOL_OBJ = $(addprefix $(OBJ_DIR)/, $(TOOL_SRC:.cpp=.o))

ALL_SRC = $(MAIN_SRC) $(TEST_SRC) $(TOOL_SRC)
ALL_OBJ = $(MAIN_OBJ) $(TEST_OBJ) $(TOOL_OBJ)
ALL_OBJ_DIR = $(sort $(BIN_DIR) $(patsubst %/,%,$(dir $(ALL_OBJ))))

$(info Default compiler and linker command lines:)
$(info $(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $$@ $$<)
$(info $(CXX) $(CXXFLAGS) $(ORIGIN) -o $$@ $$^)

.PHONY: all info build tools headers test xtest valgrind nm doxygen doxyint

info:
	@echo SRC = $(MAIN_SRC) -- $(TEST_SRC) -- $(TOOL_SRC)
	@echo OBJ = $(MAIN_OBJ) -- $(TEST_OBJ) -- $(TOOL_OBJ)
	@echo ALL_OBJ_DIR = $(ALL_OBJ_DIR)

build: $(CONFIG_FILE) $(LIBDSO) $(BIN_DIR)/test_all headers tools

headers: $(addprefix $(INC_DIR)/,$(PUBLISHED_INC))

$(INC_DIR)/%.h :: %.h | $(INC_DIR)
	/bin/cp -a -t $(INC_DIR) $<

# --no-exec is a temporary kludge. Without it, an Azure enabled build
# will crash while forking.
test: $(BIN_DIR)/test_all data | $(TESTRUNDIR)
	@echo
	$(BIN_DIR)/test_all --color=never --no-exec

XTEST=xxx
xtest: $(BIN_DIR)/test_all data
	@echo
	$(BIN_DIR)/test_all --no-exec --color=never `$(BIN_DIR)/test_all --no-exec --color=never --list | egrep '$(XTEST)'`

# Note, you might want to build with -O0 for better results.
# Note, you might want --no-exec, or maybe not.
# Note, temporarily use --gen-suppressions=all as needed.
valgrind: $(BIN_DIR)/test_all data | $(TESTRUNDIR)
	@echo
	-valgrind --log-file=$(LOG_DIR)/valgrind.txt --tool=memcheck --num-callers=12 --leak-check=full --suppressions=test/suppressions.txt $(BIN_DIR)/test_all --color=never --no-exec

tools: $(BIN_DIR)/zgycopyc $(BIN_DIR)/zgydumpc $(BIN_DIR)/zgydiffc $(BIN_DIR)/example
	@echo Tools have been built.

nm: $(LIBDSO)
	nm -gpC $(LIBDSO) | grep ' [TVW] ' | cut -d\  -f3- | sort -u

data: $(TESTDATARW)/Empty-v1.zgy $(TESTDATARW)/Empty-v3.zgy $(TESTDATARW)/Fancy-int8.zgy $(TESTDATARW)/EmptyOldFile.zgy

LATEX := $(wildcard /usr/bin/latex)
PDFTK := $(wildcard /usr/bin/pdftk)
# Uncommenting this will remove the "work in progress" background page.
PDFTK :=

doxygen:
	$(BUILDROOT)/scripts/build-doxygen.sh apidoc

doxyint:
	$(BUILDROOT)/scripts/build-doxygen.sh intdoc

$(TESTRUNDIR):
	mkdir -p $@

$(TESTDATARW)/Empty-v1.zgy: $(TESTDATA)/Empty-v1.zgy.bz2
	mkdir -p $(@D)
	bzip2 -dkc $^ > $@

$(TESTDATARW)/Empty-v3.zgy: $(TESTDATA)/Empty-v3.zgy.bz2
	mkdir -p $(@D)
	bzip2 -dkc $^ > $@

$(TESTDATARW)/Fancy-int8.zgy: $(TESTDATA)/Fancy-int8.zgy.bz2
	mkdir -p $(@D)
	bzip2 -dkc $^ > $@

$(TESTDATARW)/EmptyOldFile.zgy: $(TESTDATA)/EmptyOldFile.zgy.bz2
	mkdir -p $(@D)
	bzip2 -dkc $^ > $@

$(MAIN_OBJ) $(TEST_OBJ) $(TOOL_OBJ): $(OBJ_DIR)/%.o: %.cpp Makefile | $(ALL_OBJ_DIR)
	@echo Compile $<
	@$(CXX) -c -o $@ $(CPPFLAGS) $(CXXFLAGS) $<

$(ALLOBJ) $(LIBDSO) $(PRODUCT) $(TESTAPPS): Makefile | $(ALL_OBJ_DIR) $(BIN_DIR)

# TODO-Medium: Instead of a dependency from $(LIBDSO) to seismic store
# it would be much better to have a separate seismic store plug-in.
#
# sofixup is now done in the SDAPI build. So it isn't needed here.
#	test -h $(BIN_DIR)/libcurl.so || mv -f $(BIN_DIR)/libcurl.so $(BIN_DIR)/libcurl.so.vcpkg
#	$(BUILDROOT)/scripts/sofixup -l libcurl.so libcurl.so.4 $(BIN_DIR)/lib*.so*
#
# Note the use of $$(nm ...) in the rules instead of $(shell nm ...).
# The latter might run when the rule is parsed, not when that line
# of the rule gets executed.

$(LIBDSO): $(MAIN_OBJ) $(SD_SENTINEL) $(ZFP_SENTINEL) version.txt | $(BIN_DIR)
	$(CXX) -o $@ $(CXXFLAGS) -shared $(MAIN_OBJ) $(LIB) $(SD_LIBRARY) $(ZFP_LIBRARY) -Wl,-soname -Wl,$(notdir $(LIBDSO_LINK2)) -Wl,--no-undefined -Wl,-rpath '-Wl,$$ORIGIN'  -Wl,--version-script=version.txt
	$(RM) $(LIBDSO_SYMLINK)  $(LIBDSO_LINK2)
	/bin/ln -f -s $(notdir $(LIBDSO)) $(LIBDSO_LINK2)
	/bin/ln -f -s $(notdir $(LIBDSO_LINK2)) $(LIBDSO_SYMLINK)
	$(RM) $(BIN_DIR)/libsdapi.so* $(BIN_DIR)/libzfp.so*
ifneq ($(strip $(SD_LIBRARY)),)
	/bin/cp -f -a -t $(BIN_DIR) $(dir $(strip $(SD_LIBRARY)))/*
	ldd $(BIN_DIR)/libsdapi.so
	/bin/mkdir -p $(BUILDROOT)/build/deploy/sdapi/$(PLATFORM)
	test ! -r $(SD_ZIPFILE) || /bin/cp -p -t $(BUILDROOT)/build/deploy/sdapi/$(PLATFORM) $(SD_ZIPFILE)
endif
ifneq ($(strip $(ZFP_LIBRARY)),)
	/bin/cp -a -t $(BIN_DIR) $(strip $(ZFP_LIBRARY))
endif
	echo "OPENZGY:" $$($(BUILDROOT)/scripts/gethead.sh $(BUILDROOT)) > $(BIN_DIR)/../../version.txt
	test ! -r $(SD_BUILDDIR)/version.txt || cat $(SD_BUILDDIR)/version.txt >> $(BIN_DIR)/../../version.txt
	echo >> $(BIN_DIR)/../../version.txt
	echo "DYNAMIC $@" >> $(BIN_DIR)/../../version.txt
	ldd $@ >> $(BIN_DIR)/../../version.txt
	@echo "Exported C bindings:" $$(nm -gp $(LIBDSO) | grep ' T oz_' | wc -l)

# SD_LIBRARY is only needed because of test/sdutils.cpp doing direct access
# to SDAPI and bypassing OpenZGY entirely.
$(BIN_DIR)/test_all: $(TEST_OBJ) $(LIBDSO) $(SD_SENTINEL) $(ZFP_SENTINEL)
	$(CXX) -o $@ $(CXXFLAGS) $(ORIGIN) $(TEST_OBJ) $(LIBDSO) $(SD_LIBRARY)

#$(BIN_DIR)/zgycopyc: $(OBJ_DIR)/tools/zgycopyc.o $(OBJ_DIR)/test/mock.o $(LIBDSO)
#	$(CXX) -o $@ $(CXXFLAGS) $(ORIGIN) $^ -fopenmp
#	@echo Tool $(@F) has been built.
#
# zgycopyc is linked statically with OpenZGY so that gprof can work.
# It still links dynamically with sdapi and zfp (might change).
# DO NOT COPY/PASTE this rule unless you are writing a tool specifically
# for calling some methods in $(LIBDSO) that you need to profile.
# The declared dependency to $(LIBDSO) is wrong, but I need almost
# everything that $(LIBDSO) needs. So this is simpler to maintain.
# Note that I link with all the object files, not some libopenzgy.a
# which would need several hacks to work. Because there are object
# files that register with various factories when the code is loaded.
# Those won't be pulled in from the archive.
$(BIN_DIR)/zgycopyc: $(OBJ_DIR)/tools/zgycopyc.o $(OBJ_DIR)/tools/readwritemirror.o $(OBJ_DIR)/tools/readlodcrop.o $(OBJ_DIR)/test/mock.o $(LIBDSO)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJ_DIR)/tools/zgycopyc.o $(OBJ_DIR)/tools/readwritemirror.o $(OBJ_DIR)/tools/readlodcrop.o $(OBJ_DIR)/test/mock.o $(MAIN_OBJ) $(SD_LIBRARY) $(ZFP_LIBRARY) $(LIB)
	@echo Tool $(@F) has been built.

$(BIN_DIR)/zgydumpc $(BIN_DIR)/zgydiffc: $(BIN_DIR)/%: $(OBJ_DIR)/tools/%.o $(LIBDSO)
	@echo "Link    $@"
	@$(CXX) -o $@ $(CXXFLAGS) $(ORIGIN) $^

$(BIN_DIR)/example: example.h $(LIBDSO)
	g++ -o $@ $(CXX_ABI) -I $(INC_DIR)/.. -x c++ -Wall -std=c++11 $< -L $(BIN_DIR) $(ORIGIN) -lopenzgy
	(message=`$@ 2>&1; /bin/true`; \
	 match=`expr "$$message" : "Usage.*infile outfile"`; \
	 test $$match -gt 0)

FRC:
.PHONY: FRC

$(ALL_OBJ_DIR) $(INC_DIR):
	mkdir -p $@

.PHONY: clean clobber
clean: packages-clean
	$(RM) -rf $(OBJ_DIR)
	$(RM) $(TESTDATARW)/Empty-v1.zgy $(RM) $(TESTDATARW)/Empty-v3.zgy $(TESTDATARW)/Fancy-int8.zgy $(TESTDATARW)/EmptyOldFile.zgy

clobber: clean
	$(RM) -rf $(BIN_DIR)
	$(RM) -rf $(BUILDROOT)/build/deploy/sdapi
	$(RM) -rf $(BUILDROOT)/build/deploy/version.txt
	$(RM) -f $(DOX_DIR)/apidoc.tgz $(DOX_DIR)/apidoc.pdf
	$(RM) -f $(DOX_DIR)/intdoc.tgz $(DOX_DIR)/intdoc.pdf
	$(RM) -f $(addprefix $(INC_DIR)/,$(PUBLISHED_INC))

# The automatic dependency generator doesn't catch SDAPI / ZFP header files.
$(OBJ_DIR)/impl/file_sd.o: $(SD_SENTINEL)
$(OBJ_DIR)/test/test_file.o: $(SD_SENTINEL)
$(OBJ_DIR)/test/test_sdutils.o: $(SD_SENTINEL)
$(OBJ_DIR)/impl/compress_zfp.o: $(ZFP_SENTINEL)

# Note: If gcc finds that foo.o depends on inc.h, I will make foo.o, */foo.o
# depend on inc.h. Yes this is messy but it works and you probably wouldn't
# have noticed it without my confession.
# If a build is only housekeeping or documentation, then no compiler
# generated depenencies are needed.
REALGOALS = $(filter-out clean clobber doxygen doxyint,$(MAKECMDGOALS))
$(info Real goals are "$(REALGOALS)")
ifneq (,$(strip $(REALGOALS)))
.PHONY: depends
depends: $(OBJ_DIR)/depends.make
$(OBJ_DIR)/depends.make: $(ALL_SRC) $(SD_SENTINEL) $(ZFP_SENTINEL) | $(OBJ_DIR)
	@echo Build $(OBJ_DIR)/depends.make
	@$(CXX) -MM $(CXXFLAGS) $(ALL_SRC) > $@
	sed -i -e 's; [a-zA-Z0-9_-][a-zA-Z0-9_-]*/\.\./; ;g' $@
	sed -i -e 's;^[a-zA-Z0-9_-][a-zA-Z0-9_-]*\.o;$$(wildcard $$(OBJ_DIR)/& $$(OBJ_DIR)/*/&);' $@
include $(OBJ_DIR)/depends.make
endif
