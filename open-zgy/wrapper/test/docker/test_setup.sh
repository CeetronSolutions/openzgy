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

LINUXDISTRO="${1:-centos8}"
TAG="openzgy-qtest-${LINUXDISTRO}"
IID="${2:-/tmp/${TAG}.iid}"

############################################################
#   This file is used for manual tests. See test_run.sh.   #
#   Azure devops will use azure_outside.sh instead.        #
############################################################

# The tag and the default location of iidfile are only relevant when
# runnig the tests manually. When used as an automated test, a name
# collisions might cause some docker image to not be deleted after use.
# That is harmless.

/bin/rm -f "${IID}"
docker tag "${TAG}" "${TAG}:old" || /bin/true

case "${LINUXDISTRO}"  in
centos8)
docker build --pull -t "${TAG}" --iidfile "${IID}" - << @EOF
FROM quay.io/centos/centos:stream8
RUN dnf install -y python3 python3-devel gcc-c++ sudo virtualenv
RUN pip3 install numpy
RUN useradd -m -u 1000 -G wheel --shell /bin/bash -p '' me
USER me
ENV HOME=/home/me
@EOF
;;

focal)
docker build --pull -t "${TAG}" --iidfile "${IID}" - << @EOF
FROM ubuntu:focal
RUN apt-get update && apt-get install -y python3 python3-pip g++ sudo virtualenv
RUN pip3 install numpy
RUN useradd -m -u 1000 -G sudo --shell /bin/bash -p '' me
USER me
ENV HOME=/home/me
@EOF
;;

*)
    echo >&2 "$0: distro '${LINUXDISTRO}' is unknown."
    exit 1;;
esac

docker rmi "${TAG}:old" > /dev/null 2>&1 || /bin/true
