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

############################################################
#   This file is used for manual tests. See test_run.sh.   #
#   Azure devops will use azure_outside.sh instead.        #
#   Try to keep this file in sync with azure_outside.sh.   #
############################################################

LINUXDISTRO="${1:-centos8}"
SCRIPTDIR=$(realpath "$(dirname $0)")
BUILDROOT=$(realpath "${SCRIPTDIR}/../../..")
TAG="${2:-openzgy-qtest-${LINUXDISTRO}}"

echo "RUN $0 LINUXDISTRO=${LINUXDISTRO} TAG=${TAG}"
errors=0

for type in wheel sdist generic_wheel generic_sdist
do
    for mode in global user venv virtualenv
    do
	echo
	echo '************************************************************'
	echo '***' $mode $type
	echo '************************************************************'
	echo .../test/docker/test_inside.sh $mode $type
	docker run --rm --cap-add sys_ptrace \
	       -v ${BUILDROOT}/build/deploy:/home/me/oz/build/deploy:ro \
	       -v ${BUILDROOT}/wrapper/test:/home/me/oz/wrapper/test:ro \
	       -e OPENZGY_SDURL -e OPENZGY_SDAPIKEY -e OPENZGY_TOKEN \
	       ${TAG} /home/me/oz/wrapper/test/docker/test_inside.sh $mode $type python3
	RC=$?
	if [ "${RC}" -ne 0 ]; then
	    errors=$(expr ${errors} + 1)
	fi
	echo "Test case ${mode} ${type} returned ${RC}. Current errors: ${errors}."
    done
done
echo
echo "TOTAL ERROR COUNT: ${errors}"
echo
exit ${errors}
