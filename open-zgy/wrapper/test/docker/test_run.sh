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

# See comments in test_outside.sh about the purpose of these tests.

#############################################################################
# Note that test_run.sh and the invoked test_setup.sh and test_outside.sh
# are meant for running tests manually, with the source code bind-mounted
# into a Docker container. These scripts are currently (Jan 2025) out of date!
# For the Azure build, use azure_outside.sh which uses the pre-existing
# Dockerfile-{distro}-test images with the source code copied in.
#############################################################################

LINUXDISTRO="${1:-centos8}"
SCRIPTDIR=$(realpath "$(dirname $0)")
BUILDROOT="${SCRIPTDIR}/../../.."
IID="/tmp/openzgy-qtest.$$"

${SCRIPTDIR}/test_setup.sh   "${LINUXDISTRO}" "${IID}"
${SCRIPTDIR}/test_outside.sh "${LINUXDISTRO}" $(cat "${IID}")
RC=$?

/bin/rm -f "${IID}"

exit $RC
