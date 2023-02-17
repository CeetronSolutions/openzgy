#!/bin/bash

# This script is an alternative to running "make test" in this folder.
# Running the script is preferable when inside a build pipeline,
# as it gives better separation between build and test. It also tries
# to use only the deployed results. "make test" will actually rebuild
# any missing artifacts which means it is useless for ensuring that
# everything has been deployed properly.

set -e

# Workaround for SDAPI bug that causes 5 minutes of retries on something
# that definitely shouldn't. This may case negative tests to run way
# too long.
export OPENZGY_SD_BACKOFF=4

cd $(dirname $0)/../..

# Local tricks in my environment to get a valid SAuth token.
eval $(test ! -x private/grabtoken.sh || private/grabtoken.sh -v)

# Selectively enable or disable tests depending on what options were included.
# Note that some tests handle this themselves. E.g. test_all.
# TODO, it would be cleaner to always use the config.sh.
source build/deploy/config.sh

# I'd rather not depend on scripts/getdistro and g++ -dumpfullversion here,
# because the tests ought to be runnable even without the compiler present.
# So, just look at the deploy folder and use whatever architecture is there.
# The caller of the script can set $MYPLATFORM manually if desired.

if [ "x$MYPLATFORM" = "x" ]; then
    MYPLATFORM=$(basename $(dirname build/deploy/native/*/test_all) | awk '{print $1}')
fi

ALL_TESTDATA="
Empty-v1.zgy
Empty-v3.zgy
EmptyOldFile.zgy
Fancy-int8.zgy
"

mkdir -p   build/testdata
for file in $ALL_TESTDATA
do
    bzip2 -dkc testdata/${file}.bz2   > build/testdata/${file}
done

/bin/rm -rf build/run
/bin/mkdir -p build/run
cd build/run

TESTAPP=../../build/deploy/native/${MYPLATFORM}/test_all
TESTCMP=../../build/deploy/native/${MYPLATFORM}/zgydiffc
TESTMATRIX=../../build/deploy/native/${MYPLATFORM}/testmatrix.sh
TESTCONST=../../build/deploy/native/${MYPLATFORM}/testconstant.sh

# cutest can take a list of tests to run, but doesn't allow wildcards.
# Working around that limitation.
TESTLIST=
if [ $# -eq 1 -a "x$1" != "" -a "x$1" != "xvalgrind" ]; then
    TESTLIST=$(${TESTAPP} -color=never --no-exec --list | egrep "$1")
fi

if [ $# -eq 1 -a "x$1" = "xvalgrind" ]
then
    set +e
    mkdir -p ../../build/deploy/native
    echo "NOT RUN" > ../../build/deploy/native/valgrind.txt
    echo "Testing OpenZGY/C++ for ${MYPLATFORM}"
    echo "valgrind ${TESTAPP} --color=never --no-exec ${TESTLIST}"
    set -x
    valgrind \
	--log-file=../../build/deploy/native/valgrind.txt \
	--tool=memcheck \
	--num-callers=12 \
	--leak-check=full \
	--suppressions=../../native/src/test/suppressions.txt \
	${TESTAPP} --color=never --no-exec ${TESTLIST} || true
    pwd
    find ../../build
    set -e
else
    echo "Testing OpenZGY/C++ for ${MYPLATFORM}"
    echo "${TESTAPP} --color=never --no-exec ${TESTLIST}"
    ${TESTAPP} --color=never --no-exec ${TESTLIST}
fi
