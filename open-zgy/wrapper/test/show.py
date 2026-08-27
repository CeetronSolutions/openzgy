#!/usr/bin/env python3

#print('Running' if __name__ == '__main__' else 'Importing', __file__)

import numpy as np
import sys
import os
import time
from PIL import Image

import openzgycpp as newzgy
from openzgycpp import SampleDataType, UnitDimension, DecimationType
from test_utils import SDCredentials, TempFileAutoDelete, LocalFileAutoDelete, CloudFileAutoDelete, HasSeismicStore, HasZFPCompression, ProgressWithDots
from viewzgy import savePNG, showFileInTk

def read_data_all_at_once(reader, lod, start, size):
    section = np.zeros(size, dtype=np.float32)
    reader.read(start, section, lod = lod)
    return section

def read_data_b_at_a_time(reader, lod, start, size):
    bs = np.array(reader.bricksize, dtype=np.int64)
    padsize = ((np.array(size, np.int64) + bs - 1) // bs) * bs
    brick = np.zeros(bs, dtype=np.float32)
    section = np.zeros(padsize, dtype=np.float32)
    for ii in range(0, size[0], bs[0]):
        for jj in range(0, size[1], bs[1]):
            for kk in range(0, size[2], bs[2]):
                reader.read((start[0]+ii, start[1]+jj, start[2]+kk),brick, lod = lod)
                section[ii:ii+bs[0], jj:jj+bs[1], kk:kk+bs[2]] = brick
    return section[:size[0],:size[1],:size[2]]

def timing_report(reader, lod, start, size, elapsed):
    bs = np.array(reader.bricksize if isinstance(reader, newzgy.ZgyReader) else (64, 64, 64), dtype=np.int64)
    padsize = ((np.array(size, np.int64) + bs - 1) // bs) * bs
    bandwidth = np.prod(padsize) / elapsed # should I use size or padsize?
    bandwidth /= (1024*1024)
    print("Elapsed {0:6.2f} seconds, bandwidth {1:6.2f} MVoxel/s reading {2} lod {3} size {4} start {5}".format(elapsed, bandwidth, reader.datatype, lod, tuple(size), tuple(start)))

_zgycloud_inited = False

def run(filename, *, lods = [0], direction = 0, slurp=True,
        readerfactory = newzgy.ZgyReader, outname = None, iocontext=None):
    """
    Read 64 traces or slices in the specified direction.
    Optionally save the first of these as PNG.
    """
    if iocontext is None and filename[:5] == "sd://":
        iocontext = SDCredentials()
    if False:
        print("Read", filename,
              ("64 inlines", "64 crosslines", "64 slices")[direction],
              "slurp" if slurp else "block at a time")
    allslices = []
    with readerfactory(filename, iocontext=iocontext) as reader:
        slicenumber = reader.size[direction] // 2
        #reader.dump()
        for lod in lods:
            step = 1<<lod
            start = [0, 0, 0]
            size = np.array(reader.size, dtype=np.int64) // (1 << lod)
            start[direction] = slicenumber >> lod
            size[direction] = 1
            if direction == 2:
                size[0] = min(size[0], 1024)
                size[1] = min(size[1], 1024)
            start = tuple(start)
            size = tuple(size)
            starttime = time.time()
            if slurp:
                section = read_data_all_at_once(reader, lod, start, size)
            else:
                section = read_data_b_at_a_time(reader, lod, start, size)
            #timing_report(reader, lod, start, size, time.time() - starttime)
            # Display only the first section or slice
            if outname:
                if direction == 0:
                    myslice = section[0,...]
                elif direction == 1:
                    myslice = section[:,0,:]
                else:
                    myslice = section[...,0]
                #savePNG(myslice, outname + "_lod" + str(lod) + ".png")
                allslices.append(myslice)

    if outname:
        s = allslices[0].shape
        w = np.sum([allslices[i].shape[0] for i in range(3)])
        combined = np.zeros((w, s[1]), dtype=allslices[0].dtype)
        combined[0:s[0],:] = allslices[0]
        ss = allslices[1].shape
        combined[s[0]:s[0]+ss[0], 0:ss[1]] = allslices[1]
        sss = allslices[2].shape
        combined[s[0]+ss[0]:s[0]+ss[0]+sss[0], 0:sss[1]] = allslices[2]
        combined = combined[::-1,::]
        savePNG(combined, outname + ".png")
        showFileInTk(outname + ".png", outname)


def Main(args):

    if not args:
        args = ["../build/testdata/Empty-v3.zgy",
                "../build/testdata/Empty-v1.zgy"]

    suffix = 1
    for filename in args:
        out = os.path.join(os.getenv("TESTRUNDIR", '.'), "cpp" + str(suffix))
        try:
            run(filename, lods=range(3), direction=0, slurp=True,
                readerfactory=newzgy.ZgyReader, outname = out)
        except ZgyMissingFeature as ex:
            print("{0}: {1}".format(filename, str(ex)))

        suffix += 1

if __name__ == "__main__":
    np.seterr(all='raise')
    Main(sys.argv[1:])

# Copyright 2017-2020, Schlumberger
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
