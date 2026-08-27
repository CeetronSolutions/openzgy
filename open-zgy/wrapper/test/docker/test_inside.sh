#!/bin/bash
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

ALT="${1:?}" # global/user/venv/virtualenv
TYPE="${2:?}" # wheel/sdist/generic_wheel/generic_sdist
ARGS_BLACK="local" # local and/or sd
PYTHON="${3:?}"    # pythonX.Y, do NOT use full path.
USAGE='Usage: $0 [global|user|venv|virtualenv] [wheel|sdist] pythonX.Y'
DEPLOY=/home/me/oz/build/deploy
TESTDIR=/home/me/oz/wrapper/test

# AD-HOC tests. Need to be tweaked before being run.
#
# Check several of the ways OpenZgyBindings can be installed.
# These tests should be run inside a Docker container.
# Start a new container for each test.
# Currently the tests are expected to run manually.
#
# On CentOS 8: quay.io/centos/centos:stream8 -> quicktest
#     dnf install python3 python3-devel gcc-c++ sudo virtualenv
#     pip3 install numpy
#     useradd -m -u 1000 -G wheel --shell /bin/bash -p '' me
#
# On Ubuntu: ubuntu:focal -> quickfocal
#     apt update
#     apt install python3 python3-pip g++ sudo virtualenv
#     pip3 install numpy
#     useradd -m -u 1000 -G sudo --shell /bin/bash -p '' me
#

# Convert the Python version number to the version embedded in a wheel's name.
function wheelversion {
    python=${1:?}
    minor=$(${python} --version | sed -e 's;^Python 3.\([0-9][0-9]*\)\..*$;\1;')
    case "${minor}" in
        6) echo cp36m;;
        [0-9]|[0-9][0-9]) echo cp3${minor}-cp3${minor};;
        *) echo fix-this-error;;
    esac
}

echo "running $@"

if [ "${ALT}" = "virtualenv" -a ! -e /usr/bin/virtualenv -a -e /etc/redhat-release ]
then
    # /usr/bin/virtualenv is broken in rocky8 and Dockerfile-rocky8-test
    # will not install it. See the Dockerfile for a discusion.
    echo
    echo "Skipped because virtualenv is not installed."
    echo
    #sudo ${PYTHON} -m pip install virtualenv
    exit 0
fi

if [ ! -r /.dockerenv ]; then
    echo >&2 "This test can only be run inside Docker"
    exit 1
fi
if [ ! -d ${DEPLOY} -o ! -d ${TESTDIR} ]
then
    echo >&2 "${DEPLOY} and/or ${TESTDIR} is missing. Did you forget a -v in 'docker run'?"
    exit 1
fi

# Nothing should be pre-installed when we are testing the installers.
# TODO: Instead, use a Dockerfile that only contains the installer tarball.
# Note: Uninstall from ${PYTHON} would have been sufficient here.
# Note: Avoid "grep -q" because python complains about broken pipes.
PYTHONS=$(/bin/ls -1 /usr/bin/python3[.0-9]* | grep '/python3[.0-9]*$')
for pp in ${PYTHON} ${PYTHONS}
do
    if ${pp} -m pip list | grep OpenZgyBindings > /dev/null
    then
        echo "Uninstall OpenZgyBindings from ${pp}"
        sudo PIP_ROOT_USER_ACTION=ignore ${pp} -m pip uninstall -y OpenZgyBindings
    else
        echo "OpenZgyBindings was not pre-installed in ${pp}"
    fi
done

# Only oz/build/deploy/wrapper oz/wrapper/test is needed for the wheel
# and sdist installer. generic_wheel needs oz/build/deploy/native/{arch}/*.so*
# as well. And generic_sdist needs most of the oz folder. To improve the
# tests, it is possible (but maybe not worth the trouble) to restrict ~/oz
# to precisely what is needed. CAREFUL! If "os" has been bind-mounted into
# this container, removing unwanted files this could DESTROY the work area!
#if [ ! -e oz-unreadable ]
#then
#    pushd /home/me
#    case "$TYPE" in
#        generic_wheel)
#            tar -cf oz.tar oz/build/deploy/wrapper oz/wrapper/test oz/build/deploy/native/*-gcc*;;
#        generic_sdist)
#            tar -cf oz.tar oz;;
#        *)
#            tar -cf oz.tar oz/build/deploy/wrapper oz/wrapper/test;;
#    esac
#    mv oz oz-unreadable
#    chmod 0 oz-unreadable
#    tar -xf /home/me/oz.tar
#    popd
#fi

# generic_wheel and generic_sdist require that the SDK is pre-installed.
# generic_wheel can only pick it up from the default location.
# Unless made available using ldconfig or $LD_LIBRARY_PATH.
case "$TYPE" in
    generic_wheel)
	# Either of these should work.
	#sudo /bin/ln -s -f ${DEPLOY} /usr/local/openzgy;;
	export LD_LIBRARY_PATH=$(echo ${DEPLOY}/native/*-gcc*);;
    generic_sdist)
	export OPENZGY_SDK=${DEPLOY};;
esac

tag=$(wheelversion ${PYTHON})
echo "INFO: SDK is '${OPENZGY_SDK}'"
echo "INFO: Python is $(${PYTHON} --version) (wheel ${tag})"

case "$TYPE" in
    wheel) BINDINGS=$(/bin/ls -1d ${DEPLOY}/wrapper/*-gcc*/[Oo]pen[Zz]gy[Bb]indings-*-linux_x86_64.whl | grep ${tag});;
    sdist) BINDINGS=$(/bin/ls -1d  ${DEPLOY}/wrapper/*-gcc*/[Oo]pen[Zz]gy[Bb]indings-*.tar.gz);;
    generic_wheel) BINDINGS=$(/bin/ls -1d ${DEPLOY}/wrapper/without-dependencies/[Oo]pen[Zz]gy[Bb]indings-*-linux_x86_64.whl | grep ${tag});;
    generic_sdist) BINDINGS=$(/bin/ls -1d  ${DEPLOY}/wrapper/without-dependencies/[Oo]pen[Zz]gy[Bb]indings-*.tar.gz);;
    *) echo >&2 "$USAGE"; exit 2;;
esac

if [ -z "${BINDINGS}" -o ! -r "${BINDINGS}" ]
then
    echo >&2 "Cannot find \"${BINDINGS}\""
    echo >&2 "Here is the content of \"${DEPLOY}/wrapper\""
    find ${DEPLOY}/wrapper
    exit 1
fi

# Very rough test on package size to verify we got the right type.
# wheel and sdist contain OpenZGY binaries and may need > 30 MB.
# A prebuilt generic wheel has just one binary. 260 kB on CentOS.
# A generic sdist has just cpp and py files. 68 kB at last count.
size=$(stat -c %s ${BINDINGS})
case "$TYPE" in
    wheel|sdist)
	if [ "${size}" -lt 3000000 ]; then
	    echo >&2 "Size of ${BINDINGS} is < 3 MB. This is probably an error."
	    exit 1
	fi;;
    generic_wheel|generic_sdist)
	if [ "${size}" -gt 3000000 ]; then
	    echo >&2 "Size of ${BINDINGS} is > 3 MB. This is probably an error."
	    exit 1
	fi;;
    *) echo >&2 "Internal error in switch."
       exit 1;;
esac

cd /tmp
/bin/rm -rf vv

case "$ALT" in
    global)
	echo PIP_ROOT_USER_ACTION=ignore OPENZGY_SDK=${OPENZGY_SDK} ${PYTHON} -m pip install $pyopt $BINDINGS
	sudo PIP_ROOT_USER_ACTION=ignore OPENZGY_SDK=${OPENZGY_SDK} ${PYTHON} -m pip install $pyopt $BINDINGS
	;;

    user|local)
        # In noble, even with --user, --break-system-packages option is needed.
        # Unless it is user root that installed in --user mode. Go figure.
        # Note, used to call the "pip3" command instead of "python3 -m pip"
        # for no special reason. But, while the pip3.10 symlink exists e.g.
        # in jammy, pip3.11 does not. So the command can't pick up the right
        # version.
        echo PIP_BREAK_SYSTEM_PACKAGES=1 ${PYTHON} -m pip install --user $BINDINGS
        PIP_BREAK_SYSTEM_PACKAGES=1 ${PYTHON} -m pip install --user $BINDINGS
        ;;

    venv)
	${PYTHON} -m venv vv
	. vv/bin/activate
	pip3 install --upgrade pip
	echo pip3 install $BINDINGS
	pip3 install $BINDINGS
	pip3 install -q numpy
    ;;

    virtualenv)
	virtualenv -p ${PYTHON} vv
	. vv/bin/activate
	echo pip3 install $BINDINGS
	pip3 install $BINDINGS
	pip3 install -q numpy
    ;;

    *) echo >&2 "$USAGE"; exit 2;;
esac

echo

${PYTHON} -c 'import openzgycpp; print(openzgycpp.__file__); print(openzgycpp.wrapper.__file__)'
find / -name libopenzgy.so -print 2>/dev/null

${PYTHON} ${TESTDIR}/test_black.py "${ARGS}" && ${PYTHON} ${TESTDIR}/test.py
RC=$?
echo "Test case ${ALT} ${TYPE} returning ${RC}"
exit ${RC}
