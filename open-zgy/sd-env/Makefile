
# Point to where you cloned SDAPI, and where you want its binaries to go.

CWD       := $(shell pwd)
WORKAREA   = $(abspath $(CWD)/../seismic-store-cpp-lib)
DISTROS    = rocky8 focal jammy noble
MANUAL_ALL = $(addprefix manual-,$(DISTROS))
SDAPI_ALL = $(addprefix sdapi-,$(DISTROS))
RUN_ALL = $(addprefix run-,$(DISTROS))

.PHONY: help manual run clean clobber
#.PHONY $(MANUAL_ALL) $(RUN_ALL)
# .PHONY doesn't play well with pattern rules.
RUNARGS = -it --user 800
RUNARGS += --cap-add sys_ptrace
RUNARGS += -v $(WORKAREA):/home/me/sdapi:ro
RUNARGS += -v $(CWD)/build-template.sh:/home/me/build.sh
RUNCRED = -e OPENZGY_TOKEN -e OPENZGY_SDAPIKEY -e OPENZGY_SDURL

PULL=--pull
#PULL=

# This Makefile does not support parallel execution with -j N
# The only really time consuming part is to build the Azure
# prerequisites, and the task for that will enable parallel itself.
.NOTPARALLEL:

help:
	@#
	@# NOTE: The help text in Makefile should be kept in sync with README.md
	@#
	@echo "Make targets for unattended builds. This is the normal usage."
	@echo
	@echo "Output is stored as ../sdapi-pkgs/DIST-sdapi_linux64.tar.gz"
	@echo "It needs to be copied to ../pkg/sdapi_linux64.tar.gz before"
	@echo "building OpenZGY for DIST. The files in ../pkg are not named"
	@echo "after the distro, so you can only have one set in that folder."
	@echo
	@echo "    make sdapi      -- Build SDAPI only"
	@echo "    make sdapi-DIST -- Build SDAPI only, for DIST only"
	@echo
	@echo "Make targets for manual builds, typically run in this order:"
	@echo
	@echo "    make clobber    -- Removes docker instances but not images"
	@echo "    make manual     -- Rebuilds images used for manual builds"
	@echo "                      Source mapped from $(WORKAREA)"
	@echo "    make run        -- Create instances but don't start them"
	@echo "    make xterm      -- Start instances in individual windows"
	@echo "    ./build.sh      -- to be executed inside each xterm"
	@echo
	@echo 'To manually build one distro only, use e.g. manual-focal and'
	@echo 'run-focal. Then start with "docker start -ai sd-focal"'
	@echo 'instead of "make xterm".'
	@echo
	@echo 'Note that manual builds expect to find the SDAPI source at'
	@echo "$(WORKAREA)"
	@echo 'in the host. The source will be mapped as a volume into the'
	@echo 'container. Any patching needs to be done outside. "make sdapi"'
	@echo 'will instead copy the source and will handle patching itself.'
	@echo

manual: $(MANUAL_ALL)

sdapi: $(SDAPI_ALL)

run: $(RUN_ALL)

FRC:

manual-%: FRC
	docker tag sd-env-$*:manual sd-env-$*:manual-old || true
	docker build $(PULL) -f Dockerfile-$* --target sdapi-manual -t sd-env-$*:manual ..
	docker rmi sd-env-$*:manual-old || true

sdapi-%: FRC
	docker tag sd-env-$*:sdapi sd-env-$*:sdapi-old || true
	docker build $(PULL) -f Dockerfile-$* --target sdapi-build -t sd-env-$*:sdapi ..
	docker rmi sd-env-$*:sdapi-old || true
	( \
	  cid=$$(docker create sd-env-$*:sdapi); \
	  tarball=sdapi_linux64.tar.gz; \
	  mkdir -p ../sdapi-pkgs; \
	  docker cp $${cid}:/home/me/$${tarball} ../sdapi-pkgs/$*_$${tarball}; \
	  docker rm $${cid} \
	)

run-%: FRC
	docker rm --force sd-$* 2>/dev/null || true
	docker create $(RUNARGS) --name sd-$*  sd-env-$*:manual

.PHONY: xterm
xterm:
	xterm -sl 256  -T 'rocky8 (gcc 8.5.0)'       -geometry 200x24+100+200 -e docker start -ai sd-rocky8 &
	@sleep 2
	xterm -sl 256  -T 'Ubuntu focal (gcc 9.3.0)' -geometry 200x24+150+300 -e docker start -ai sd-focal &
	@sleep 2
	xterm -sl 256  -T 'Ububtu jammy (11.4)'      -geometry 200x24+200+400 -e docker start -ai sd-jammy &
	@sleep 2
	xterm -sl 256  -T 'Ububtu noble (gcc 13.2)'  -geometry 200x24+250+500 -e docker start -ai sd-noble &

# Caveat: this only works on local builds where the timestamp can be preseved.
# Caveat: Make sure all images used in FROM lines are listed here.
# Pull all the base images, but not more frequently than every 24 hours.
# If this is called before any container build then the --pull option
# can be removed from "docker build". Docker is now rate limiting pulls
# to 100 per 6 hours, and apparently they also count pulls that do nothing
# because the image is already up to date.

.PHONY: pull
pull:
	touch -d yesterday pull-needed.timestamp
	$(MAKE) pull-done.timestamp

pull-done.timestamp: pull-needed.timestamp
	docker pull rockylinux/rockylinux:8
	docker pull ubuntu:focal
	docker pull ubuntu:jammy
	docker pull ubuntu:noble
	touch pull-done.timestamp

clobber clean:
	docker rm --force $(patsubst %,sd-%,$(DISTROS)) || true
