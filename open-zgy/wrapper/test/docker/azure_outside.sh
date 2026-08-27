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

# This test script is used to try out different ways the Python wrapper
# might be installed: Native package manager vs. self contained wheel,
# sdist vs. wheel, global vs. user vs. virtualenv vs. venv.
# A small subset of the regular unit tests are run for each test case.
#
# The script can exercise all 16 combinations. If this is too much hassle
# then run the 4 combinations of (wheel, generic_sdist) and (user, venv).
# Those are the ones most likely to be used by end users.

#############################################################################
# THIS VERSION OF THE SCRIPT IS FOR USE BY THE AZURE DEVOPS PIPELINE.
# Use test_run.sh and test_outside.sh for running tests manually.
# And try to keep the scripts in sync where that makes sense.
#############################################################################

LINUXDISTRO="${1:?}"
SCRIPTDIR=$(realpath "$(dirname $0)")
BUILDROOT=$(realpath "${SCRIPTDIR}/../../..")
TAG="${2:?}"

echo "RUN $0 LINUXDISTRO=${LINUXDISTRO} TAG=${TAG}"
errors=0

# Run the tests for all versions of Python that are installed in the
# test image. Need to ask the image for a list. It had been simpler to
# have the loop(s) inside the test environment.

#PYTHONS=$(/bin/ls -1 /usr/bin/python3[.0-9]* | sed -ne 's;/usr/bin/\(python3[.0-9]*\)$;\1;p')
PYTHONS=$(docker run --rm ${TAG} /bin/bash -c '/bin/ls -1 /usr/bin/python3[.0-9]* | sed -ne "s;/usr/bin/\(python3[.0-9]*\)$;\1;p"')
echo "Python versions to test: ${PYTHONS}"

for PYTHON in ${PYTHONS}
do
  echo -e "\nRunning smoke tests for ${PYTHON}\n"
  for type in wheel sdist generic_wheel generic_sdist
  do
    for mode in global user venv virtualenv
    do
	echo
	echo '************************************************************'
	echo '***' $mode $type ${PYTHON}
	echo '************************************************************'
	echo .../test/docker/test_inside.sh $mode $type
	docker run --rm --cap-add sys_ptrace \
	       ${TAG} \
	       /home/me/oz/wrapper/test/docker/test_inside.sh $mode $type $PYTHON
	RC=$?
	if [ "${RC}" -ne 0 ]; then
	    errors=$(expr ${errors} + 1)
	fi
	echo "Test case ${mode} ${type} returned ${RC}. Current errors: ${errors}."
    done
  done
done
echo
echo "TOTAL ERROR COUNT: ${errors}"
echo
exit ${errors}
