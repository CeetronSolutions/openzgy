#!/bin/bash

# This script is an alternative to running "make test" in this folder.
# Running the script is preferable when inside a build pipeline,
# as it gives better separation between build and test. It also tries
# to use only the deployed results. "make test" will actually rebuild
# any missing artifacts which means it is useless for ensuring that
# everything has been deployed properly.

set -e

cd $(dirname $0)/..

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

# Local tricks in my environment to get a valid SAuth token.
eval $(test ! -x private/grabtoken.sh || private/grabtoken.sh -v)

# I'd rather not depend on scripts/getdistro and g++ -dumpfullversion here,
# because the tests ought to be runnable even without the compiler present.
# So, just look at the deploy folder and use whatever architecture is there.
# The caller of the script can set $MYPLATFORM manually if desired.

if [ "x$MYPLATFORM" = "x" ]; then
    # All should use the same distro, only the python version differs.
    # So, just pick the first one.
    BINDING0=$(/bin/ls -1 build/deploy/wrapper/*-gcc*/OpenZgyBindings-*.whl | head -1)
    MYPLATFORM=$(basename $(dirname ${BINDING0}))
    if [ -d "build/deploy/wrapper/${MYPLATFORM}" ]; then
	echo "MYPLATFORM=${MYPLATFORM} (from deploy)"
    else
	echo >&2 "ERROR, cannot determine distro name."
	exit 3
    fi
else
    echo "MYPLATFORM=${MYPLATFORM} (explicit)"
fi

/bin/rm -rf build/run
/bin/mkdir -p build/run
cd build/run

#PIP="python3 -m pip --disable-pip-version-check -q"
PIP="../../scripts/pip3.sh --disable-pip-version-check -q"

rm -f -rf test-venv
# TODO: Explicit choice of python version.
python3 -m venv test-venv
. test-venv/bin/activate
${PIP} install --upgrade pip
${PIP} install --upgrade setuptools wheel
${PIP} install -r ../../build/deploy/wrapper/${MYPLATFORM}/requirements.txt
tag=$(wheelversion python3)
for wheel in ../../build/deploy/wrapper/${MYPLATFORM}/OpenZgyBindings-*.whl
do
    if echo ${wheel} | /bin/grep -q ${tag}
    then
        echo "Try to install $(basename ${wheel}) in $(python --version)"
        ${PIP} install ${wheel} || echo "Failure noted and ignored."
    else
        echo "Skip install of ${wheel}. Wrong version."
    fi
done
TESTLIST=local
if [ -r ../../build/deploy/libsdapi.so.3 -a "x${OPENZGY_TOKEN}" != "x" ]
then
    TESTLIST="$TESTLIST sd"
fi

echo "Python environment for running tests (${TESTLIST}):"
${PIP} list

if [ $# -eq 1 -a "x$1" = "xvalgrind" ]
then
    set +e
    # Half-hearted attempt to get python to play nice with valgrind.
    # YMMV. Consider building a special version of Python instead.
    export PYTHONMALLOC=malloc
    mkdir -p ../../build/deploy/wrapper
    echo "valgrind ../../wrapper/test/test.py ${TESTLIST}"
    valgrind \
	--log-file=../../build/deploy/wrapper/valgrind-1.txt \
	--tool=memcheck \
	--num-callers=12 \
	--leak-check=full \
	--suppressions=../../native/src/test/suppressions.txt \
	python3 ../../wrapper/test/test.py
    echo "valgrind ../../wrapper/test/test_black.py ${TESTLIST}"
    valgrind \
	--log-file=../../build/deploy/wrapper/valgrind-2.txt \
	--tool=memcheck \
	--num-callers=12 \
	--leak-check=full \
	--suppressions=../../native/src/test/suppressions.txt \
	python3 ../../wrapper/test/test_black.py
    set -e
else
    echo python3 ../../wrapper/test/test.py
         python3 ../../wrapper/test/test.py
    echo python3 ../../wrapper/test/test_pickle.py
         python3 ../../wrapper/test/test_pickle.py
    echo python3 ../../wrapper/test/test_black.py ${TESTLIST}
         python3 ../../wrapper/test/test_black.py ${TESTLIST}
fi
