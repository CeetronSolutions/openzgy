#!/usr/bin/env python3

"""

Example code using the Python wrapper around the C++ OpenZGY API.
Install the Python API as follows:

    $ python3 -m venv ~/vv
    $ source ~/vv/bin/activate
    $ python3 -m pip install --upgrade pip
    $ python3 -m pip install -r build/deploy/wrapper/*-gcc*/requirements.txt
    $ python3 -m pip install build/deploy/wrapper/*-gcc*/OpenZgyBindings-*.whl

Run the script:

    $ python3 wrapper/test/example.py somefile.zgy someoutput.zgy

"""

import numpy as np
import sys
from openzgycpp import ZgyReader, ZgyWriter
from openzgycpp import SampleDataType, UnitDimension, DecimationType
from test_utils import SDCredentials, ProgressWithDots

def roundup(n, align):
    return ((n + align - 1) // align) * align

def copyloop(r, w):
    """Copy data between two open files, one brick-column at a time."""
    size  = r.size
    bs    = r.bricksize
    data  = np.zeros((bs[0], bs[1], roundup(size[2], bs[2])), dtype=np.float32)
    total = ((size[0]+bs[0]-1) // bs[0]) * ((size[1]+bs[1]-1) // bs[1])
    done  = 0
    p1    = ProgressWithDots()
    p2    = ProgressWithDots()
    for ii in range(0, size[0], bs[0]):
        for jj in range(0, size[1], bs[1]):
            r.read((ii, jj, 0), data);
            w.write((ii,jj, 0), data);
            done += 1
            p1(done, total)
    w.finalize(progress=p2)

def copy(srcfilename, dstfilename):
    with ZgyReader(srcfilename, iocontext = SDCredentials()) as r:
        with ZgyWriter(dstfilename,
                       size        = r.size,
                       iocontext   = SDCredentials(),
                       zfp_compressor = 50,
                       bricksize   = r.bricksize,
                       datatype    = SampleDataType.float,
                       zunitdim    = r.zunitdim,
                       zunitname   = r.zunitname,
                       zunitfactor = r.zunitfactor,
                       hunitdim    = r.hunitdim,
                       hunitname   = r.hunitname,
                       hunitfactor = r.hunitfactor,
                       zstart      = r.zstart,
                       zinc        = r.zinc,
                       annotstart  = r.annotstart,
                       annotinc    = r.annotinc,
                       corners     = r.corners) as w:
            copyloop(r, w)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: example infile.zgy outfile.zgy")
        sys.exit(1)
    copy(sys.argv[1], sys.argv[2])

# Copyright 2024-2024, Schlumberger
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
