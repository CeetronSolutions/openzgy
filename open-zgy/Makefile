# Copyright 2017-2021, Schlumberger
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
include scripts/Make.defs

TAG = build-openzgy-$(LINUXDISTRO)

OPENZGY_TOKEN ?= $(firstword $(strip $(shell sdutil auth idtoken || true)))
export OPENZGY_TOKEN := $(OPENZGY_TOKEN)
NEWTOKEN := $(shell test ! -x private/refreshtoken.sh || private/refreshtoken.sh)
$(if $(NEWTOKEN),$(info Refreshed the SAuth token used for unit tests),$(info Makefile did not refresh SAuth token))
$(if $(NEWTOKEN),$(eval export override $(NEWTOKEN)))
undefine NEWTOKEN

.PHONY: native wrapper build clean clobber docker-build docker-buildtest docker-test docker noop doxygen

# This Makefile does not support parallel execution with -j N.
#
# Make would conclude that it may choose to build native and wrapper
# in parallel. But the "wrapper" target will build
# (but not test) the "native" library using a recursive make. While
# the "native" target does a similar recursive build, with tests.
# There are probably others. It is possible to restructure the rules
# so Make knows about those dependencies. But that is hardly useful
# because most rules in this file would end up being serialized
# anyway.
#
# The next declaration causes any "-j N" option to be ingored in this
# Makefile but it will still be honored by recursive makes. Unless
# those sub-makefiles choose to do the same.

.NOTPARALLEL:

all: native wrapper

noop:
	@echo no-op

native: | $(TESTRUNDIR)
	$(MAKE) -C native/src all-no-doxygen

wrapper: | $(TESTRUNDIR)
	$(MAKE) -C native/src build
	$(MAKE) -C wrapper generic specific

build: | $(TESTRUNDIR)
	$(MAKE) -C native/src data build
	$(MAKE) -C wrapper generic_build specific_build venv

doxygen: | $(TESTRUNDIR)
	$(MAKE) -C native/src doxygen doxyint

# Collect the source files needed to run all the tests and deploy them
# as a single tar file. This is useful if you want to run the tests
# in an environment where you only have the deployed files available.
# Don't include folders in the tar file, only regular files.
# The gitlab runner doesn't allow chmod or touch on a folder.
.PHONY: testscripts
testscripts:
	tar cvf build/deploy/testscripts.tar \
	    scripts/pip3.sh \
	    scripts/python3pip.sh \
	    native/src/runtests.sh \
	    wrapper/runtests.sh \
	    wrapper/test/test.py \
	    wrapper/test/test_black.py \
	    wrapper/test/test_pickle.py \
	    wrapper/test/test_utils.py \
	    wrapper/test/docker/test_inside.sh \
	    wrapper/test/docker/azure_outside.sh \
	    native/src/test/suppressions.txt \
	    testdata/Empty-v1.zgy.bz2 \
	    testdata/Empty-v3.zgy.bz2 \
	    testdata/Fancy-int8.zgy.bz2 \
	    testdata/EmptyOldFile.zgy.bz2 \
	    $(wildcard private/refreshtoken.sh*) \
	    $(wildcard private/grabtoken.sh*) \
	    $(wildcard private/showtoken.py*)

# Note, this recipe will be duplicated in the Azure build script
# instead of Azure running "make", which would give an additional
# dependency for the Azure agent. Also, the rule below assumes
# the developer has seismic-service-bin in the work area, or a symlink.
# The azure task will instead fetch it every time.
# The Azure task cannot use docker run -v ... because the agent itself
# might be running as a docker container.

# TODO-Low, cannot run more than one build of the exact same env
# because the docker names would collide.

DOCKERFILE   := $(wildcard scripts/Dockerfile-$(LINUXDISTRO))
DOCKERTARGET := $(if $(findstring xdoc, x$(strip $(LINUXDISTRO))),doxygen,build testscripts)

docker-build:
	$(if $(DOCKERFILE),,$(error LINUXDISTRO=$(LINUXDISTRO) is not supported))
	docker tag $(TAG) $(TAG):old || true
	mkdir -p seismic-service-bin private
	docker build --pull -t $(TAG) -f $(DOCKERFILE) .
	rmdir --ignore-fail-on-non-empty seismic-service-bin private
	docker rmi -f $(TAG):old $(TAG):test || true
	$(RM) cid.txt
	docker run --cidfile cid.txt $(TAG) make -j 8 $(DOCKERTARGET)
	docker cp $$(cat cid.txt):/home/me/oz/build/deploy/. - | gzip -9 > openzgy-install.tgz
	docker rm $$(cat cid.txt)
	$(RM) cid.txt
	docker build -t $(TAG):test -f $(DOCKERFILE)-test .

# Only valid right after docker-build, useful if experimenting
# with changes to $(DOCKERFILE)-test and nothing else.
docker-buildtest:
	docker build -t $(TAG):test -f $(DOCKERFILE)-test .

docker-test:
	@echo -e '\n*** Test native started at' $$(date)
	docker run --rm -e OPENZGY_TOKEN -e OPENZGY_SDURL -e OPENZGY_SDAPIKEY $(TAG):test native/src/runtests.sh
	@echo -e '\n*** Test wrapper started at' $$(date)
	docker run --rm -e OPENZGY_TOKEN -e OPENZGY_SDURL -e OPENZGY_SDAPIKEY $(TAG):test wrapper/runtests.sh
	@echo -e '*** Tests done at' $$(date) '\n'

docker: docker-build docker-test
	@true

.PHONY: flexnet
flexnet:
	zip -r OpenZGY.zip README.md wrapper native scripts Makefile

clean:
	$(MAKE) -C wrapper $@
	$(MAKE) -C native/src $@

clobber:
	$(MAKE) -C wrapper $@
	$(MAKE) -C native/src $@
	$(RM) -rf $(TESTRUNDIR)
	$(RM) pure.tar native.tar wrap.tar deploy.tgz OpenZGY.zip
	$(RM) build/deploy/testscripts.tar
	@echo Checking that ./build is empty after clobbering
	test ! -d build || find build -type f -print
	test ! -d build || test `find build -type f -print | wc -l` -eq 0

$(TESTRUNDIR):
	mkdir -p $(TESTRUNDIR)

.PHONY: todo-native todo eyesore
SELECT  = $(addprefix --include=,\*.py \*.cpp \*.h \*.txt)
SELECT += --exclude-dir test-venv
COUNT1  = sed -e 's/^.*\(TODO[^ ,.:]*\).*/\1/' | sort | uniq -c
COUNT2  = sed -e 's/^.*\(TODO[^ ,.:]*\).*/(sum)/' | sort | uniq -c

todo-native:
	@echo
	@echo Native
	@grep -r $(SELECT) TODO native/src wrapper | $(COUNT1)
	@grep -r $(SELECT) TODO native/src wrapper | $(COUNT2)

todo-plans:
	@echo
	@echo Plans
	@grep -r $(SELECT) TODO plans.txt | $(COUNT1)
	@grep -r $(SELECT) TODO plans.txt | $(COUNT2)

NEED_EYESORE = *.* native/src wrapper
SELECT2 = $(SELECT) --include Makefile --include \*.defs
eyesore:
	@echo
	grep -rl $(SELECT2) 'http://www.apache.org/licenses/LICENSE-2.0' $(NEED_EYESORE) | sort -o eyesore_1
	grep -rl $(SELECT2) 'WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND' $(NEED_EYESORE) | sort -o eyesore_2
	grep -rl $(SELECT2) '.' $(NEED_EYESORE) | sort -o eyesore_3
	comm -12 eyesore_1 eyesore_2 > eyesore_4
	comm -23 eyesore_3 eyesore_4 > eyesore_missing
	$(RM) eyesore_1 eyesore_2 eyesore_3 eyesore_4
	@echo "The following files do not have a copyright notice:"
	@echo
	@cat eyesore_missing

todo: todo-plans todo-native
	@/bin/true
