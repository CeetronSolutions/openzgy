#!/bin/bash
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

# Build the software and report code coverage.
# Coverage information can be viewed in a browser.
# Usage: Run this script as a replacement for
# "make build testscripts".

# If running the script manually:
#docker run -it build-openzgy-rocky8

# if running the script in batch, fine tuning the script without
# wanting to build a new base image every time:
#docker run -v ~/oz/scripts/build-coverage.sh:/home/me/oz/scripts/build-coverage-tmp.sh:ro build-openzgy-rocky8 /bin/bash scripts/build-coverage-tmp.sh

# Normally the script is found in ${BASELINE}/scripts/
# If running the script manually the path may need to be adjusted.
BASELINE=$(realpath $(dirname ${0})/..)
echo "Building coverage report for ${BASELINE}"

# Install additional software needed for coverage report.
sudo dnf install -y perl perl-XML-Simple perl-XML-Dumper gd-devel perl-GD

# Build the software, creating .gcno files for the C++ parts.
cd ${BASELINE}
make OPTIMIZE="-g -fprofile-arcs -ftest-coverage" clobber build testscripts

cd ${BASELINE}/wrapper
rm -rf gcov-venv; python3.8 -m venv gcov-venv; source gcov-venv/bin/activate
pip3 install --upgrade pip
pip3 install -r requirements.txt
ls -l ../build/deploy/wrapper/rocky8-gcc85 || true
pip3 install ../build/deploy/wrapper/rocky8-gcc85/OpenZgyBindings-*-cp38-cp38-linux_x86_64.whl
deactivate

# If running manually to verify this script, make a checkpoint.
# docker commit `docker ps -lq` lcov1

# Make sure functions in untouched code are shown (with zero counts).
cd ${BASELINE}/native/src
lcov --directory ../../build/temp/native --base-directory . --zerocounters
lcov --directory ../../build/temp/native --base-directory . --capture --initial -o initial.info --include '*native/src/*' --exclude '*native/src/test/*' --exclude '*native/src/tools/*'

cd ${BASELINE}
# Local tricks in my environment to get a valid SAuth token.
eval $(test ! -x private/grabtoken.sh || private/grabtoken.sh -v)

#############################################################################
###   Run the actual tests, creating .gcda files.   #########################
#############################################################################

cd ${BASELINE}/native/src
FAIL=0

# Exclude some tests that are problematic in conjunction with coverage.
#  - timer.overhead -- takes too long when coverage or valgrind is enabled.
ALL=$(${BASELINE}/build/deploy/native/rocky8-gcc85/test_all --list |
          egrep -v '^Unit tests:|~|timer.overhead')
echo "Native tests:" ${ALL} | fold -s

${BASELINE}/build/deploy/native/rocky8-gcc85/test_all --no-exec ${ALL}
test $? -eq 0 -o "${FAIL}" -ne -0 || FAIL=1

# Run all tests again, with all logging code enabled but output discarded.
env OPENZGY_VERBOSE=42 OPENZGY_RECORD_LOGFILE=/dev/null \
    ${BASELINE}/build/deploy/native/rocky8-gcc85/test_all --no-exec ${ALL}
test $? -eq 0 -o "${FAIL}" -ne -0 || FAIL=2

# TODO, cursory check of the recorder output.

# If I don't do this, coverage is reported as too low because of debug code.
# If I do it, coverage might be reported as too high because some attributes
# may be read that would otherwise not be accessed and presumably checked.

env OPENZGY_MEASURE_KB=256 OPENZGY_MEASURE_LOGFILE=/dev/null ${BASELINE}/build/deploy/native/rocky8-gcc85/test_all --no-exec .verbose meta.open_v1 meta.open_v3 file.localfilefactory > /dev/null 2>&1
test $? -eq 0 -o "${FAIL}" -ne -0 || FAIL=3

#############################################################################
###   Done running native tests.   ##########################################
#############################################################################

source ../../wrapper/gcov-venv/bin/activate

echo 'Activated for running Python tests'
python3 -m pip list

echo python3 ../../wrapper/test/test.py
python3 ../../wrapper/test/test.py
test $? -eq 0 -o "${FAIL}" -ne -0 || FAIL=4

echo python3 ../../wrapper/test/test_pickle.py
python3 ../../wrapper/test/test_pickle.py
test $? -eq 0 -o "${FAIL}" -ne -0 || FAIL=5

echo python3 ../../wrapper/test/test_black.py local
python3 ../../wrapper/test/test_black.py local
test $? -eq 0 -o "${FAIL}" -ne -0 || FAIL=6

deactivate

#############################################################################
###   Done running Python tests.   ##########################################
#############################################################################

# If running manually to verify this script, make a checkpoint.
# docker commit `docker ps -lq` lcov1

# Process the collected data into .gcov format.
cd ${BASELINE}/native/src
for file in *.cpp impl/*.cpp; do gcov --relative-only --object-dir ../../build/temp/native/$(dirname $file) $file; done

# Further process using lcov.
lcov --directory ../../build/temp/native --base-directory . --capture -o module.info --include '*native/src/*' --exclude '*native/src/test/*'

# Add in the zero counters to make sure we report everything
lcov -o pluszero.info -a initial.info -a module.info

# --include and --relative and --extract won't work. --remove might work.
lcov -o final.info -r pluszero.info '/usr/include/*' '*/native/src/test/*' '*/native/src/tools/*' | tee final.txt
tail -4 < final.txt > summary.txt

# And then to html
rm -rf coverage; genhtml --num-spaces 8 --title OpenZGY final.info -o coverage
tar zcf ../../build/deploy/native/coverage.tgz coverage

# If running manually, exit the container and retrieve the result.
#rm -rf coverage coverage.tgz; docker cp `docker ps -lq`:/home/me/oz/build/deploy/native/coverage.tgz .; tar zxf coverage.tgz; firefox coverage/index.html &

if [ ${FAIL} -ne 0 ]
then
    echo >&2 "FAILED: Test suite ${FAIL} and possibly others."
    # It is debatable whether to let the pipeline fail and not
    # deliver any coverage results at all. If this happens too
    # often, it migh be better to ignore the error.
    exit 1
fi
