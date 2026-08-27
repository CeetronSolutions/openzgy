#!/usr/bin/env python3

"""
Test ZGY bindings for Python.
"""

import openzgycpp as zgy
import numpy as np
import os
from test_utils import LocalFileAutoDelete

def dump(message, obj, verbose = False):
    if message: print(message)
    for x in dir(obj):
        if x[0] != '_':
            if verbose:
                doc = '\n' + getattr(obj.__class__, x).__doc__
                doc = doc.replace('\n', '\n\t\t')
                value = getattr(obj, x)
                print('\t' + x, "=", value, doc)
            else:
                value = getattr(obj, x)
                if not callable(value):
                    print('\t' + x, "=", value)

def checkmeta(meta, datatype = zgy.SampleDataType.float32, complete = True):
    assert(meta.size == (30, 60, 100))
    assert(meta.datatype == datatype)
    # This works differently in OpenZGY than in ZGY-Public.
    # In OpenZGY the user specified datarange is ignored for float.
    # The range will be empty if no data exists yet, or the actual
    # value range if it does. Similarly, if cloning a float cube
    # into an integer one the actual value range of the float cube
    # will be used as the coding range.
    #assert(meta.datarange == (-1, 1))
    assert(meta.zunitdim == zgy.UnitDimension.time)
    assert(meta.zunitname == "ms")
    assert(abs(meta.zunitfactor - 0.001) < 1.0e-5)
    assert(meta.hunitdim == zgy.UnitDimension.length)
    assert(meta.hunitname == "ft")
    assert(abs(meta.hunitfactor - 0.3048) < 0.0001)
    assert(meta.zstart == 2500)
    assert(abs(meta.zinc - 4.1) < 0.0001)
    assert(meta.annotstart == (1234, 5678))
    assert(meta.annotinc == (5, 2))
    assert(meta.corners == ((1000, 1000), (3900, 1000), (1000, 6900), (3900, 6900)))
    assert(meta.nlods == 2 if complete else 1)

def createFileFromScratch(filename, unlock_gil=False):
    with zgy.ZgyWriter(filename,
                       size = (30, 60, 100),
                       datatype = zgy.SampleDataType.float32,
                       datarange = (-1, 1),
                       zunitdim = zgy.UnitDimension.time,
                       zunitname = "ms",
                       zunitfactor = 0.001,
                       hunitdim = zgy.UnitDimension.length,
                       hunitname = "ft",
                       hunitfactor = 0.3048,
                       zstart = 2500,
                       zinc = 4.1,
                       annotstart = (1234, 5678),
                       annotinc = (5, 2),
                       corners = ((1000, 1000), (3900, 1000), (1000, 6900), (3900, 6900)),
                       unlockgil=unlock_gil,
    )as demo:
        #dump("Unit test 1:", demo)
        checkmeta(demo, complete=False)
        data = np.random.standard_normal(demo.size).astype(np.float32)
        demo.write((0, 0, 0), data)

def readBackFirstFile(filename, unlock_gil=False):
    with zgy.ZgyReader(filename, unlockgil=unlock_gil) as demo:
        #dump("Unit test 1:", demo)
        checkmeta(demo)
        orig = ( 0, 0, 0)
        size = ( 1, demo.size[1], demo.size[2] )
        data = np.zeros(size, dtype=np.float32)
        demo.read(orig, data)
        #print(data)
        avg = np.average(data.flat)
        #print("Average :", avg)
        assert(avg >= -1 and avg <= 1)
        assert(np.min(data.flat) != np.max(data.flat))

def createFileAsClone(filename, clonefrom):
    with zgy.ZgyWriter(filename, clonefrom, datatype = zgy.SampleDataType.int16) as demo:
        #dump("Unit test 2:", demo)
        checkmeta(demo, zgy.SampleDataType.int16, complete=False)

def createFileAsClone8(filename, clonefrom, unlock_gil=False):
    with zgy.ZgyWriter(filename, clonefrom, datatype = zgy.SampleDataType.int8, datarange = (-128, +127), unlockgil=unlock_gil) as demo:
        #dump("Unit test 3:", demo)
        data = (np.random.standard_normal(demo.size) * 100).astype(np.int8)
        demo.write((0, 0, 0), data)
    with zgy.ZgyReader(filename) as demo:
        # Read back as float.
        data = np.zeros(demo.size, dtype=np.float32)
        demo.read((0,0,0), data, lod=0)
        #print(data)
        # This is an 8-bit file, so we can also read it as such.
        data = np.zeros(demo.size, dtype=np.int8)
        demo.read((0,0,0), data, lod=0)
        #print(data)

def readBackSecondFile(filename):
    # Open the file we just created and check meta data.
    # Note that we didn't actually write any data to this one.
    with zgy.ZgyReader(filename) as demo:
        #dump("Unit test 2:", demo)
        checkmeta(demo, zgy.SampleDataType.int16)
        orig = ( 0, 0, 0)
        size = ( 1, demo.size[1], demo.size[2] )
        try:
            data = np.zeros(size, dtype=np.float32)
            demo.read(orig, data)
            #print(data)
            # The latest version handles 'brick not found' itself.
            # So 'no exception' is a valid result.
            #assert("Should have gotten an exception" == "!")
        except zgy.ZgyError as ex:
            # Missing data is normally replaced with zero,
            # but if the requested range was completely outside
            # the valid data we get 'Brick not found'.
            #print("Got expected exception ", ex)
            assert(ex.args[0] == 16)
            assert(ex.args[1] == 'brick not found')

def checkConversions(filename):
    with zgy.ZgyReader(filename) as demo:
        #dump("", demo, True)
        a = demo.indexToAnnot((3, 7))
        i = demo.annotToIndex(a)
        #print(a, i)
        assert(a == (1249, 5692) and i == (3, 7))
        w = demo.indexToWorld((3, 7))
        i = demo.worldToIndex(w)
        #print(w, i)
        assert(w == (1300, 1700) and i == (3, 7))
        w = demo.annotToWorld(a)
        a = demo.worldToAnnot(w)
        #print(w, a)
        assert(w == (1300, 1700) and (1249, 5692))

#def readUsingAnotherAPI(filename):
#    with zgy.reader(filename) as demo:
#        dump("Unit test 3:", demo)
#        checkmeta(demo)

def readFromClosedFile(filename):
    with zgy.ZgyReader(filename) as reader:
        pass
        data = np.zeros((1,1,1), dtype=np.float32)
    try:
        reader.read((0, 0, 0), data)
        assert(False) # Should have raised an exception
    except zgy.ZgyError as ex:
        assert "Not open" in str(ex)

def dumpAPI(m):
    for name in sorted(dir(m)):
        if name[0] != "_":
            print(name)
            instance = getattr(zgy, name)
            if type(instance) == type:
                print("  " + "\n  ".join([e for e in sorted(dir(instance)) if e[0] != "_"]))

def dumpInstallLocation():
    """For debugging setup.py."""
    print(zgy.__file__)
    print(zgy.wrapper.__file__)
    if os.name != 'nt':
        os.system("ldd " + zgy.wrapper.__file__ + " | sed -ne '/libopenzgy/s/^.* => \\(.*\\) .*$/\\1/p'")
        os.system("readelf -d " + zgy.wrapper.__file__ + " | grep PATH")

dumpInstallLocation()
#dumpAPI(zgy)
#with zgy.ZgyReader("../../build/testdata/Empty-v3.zgy") as r:
#    print(r.meta)
#print(list(zgy.SampleDataType.__members__))
#print(list(zgy.UnitDimension.__members__))
#print(list(zgy.DecimationType.__members__))

with    LocalFileAutoDelete("unittest1.zgy") as lad1, \
        LocalFileAutoDelete("unittest2.zgy") as lad2, \
        LocalFileAutoDelete("unittest3.zgy") as lad3, \
        LocalFileAutoDelete("unittest4.zgy") as lad4, \
        LocalFileAutoDelete("unittest5.zgy") as lad5:

    createFileFromScratch(lad1.name)
    readBackFirstFile(lad1.name)
    createFileAsClone(lad2.name, lad1.name)
    readBackSecondFile(lad2.name)
    checkConversions(lad1.name)
    createFileAsClone8(lad3.name, lad1.name)
    #readUsingAnotherAPI(lad1.name)
    readFromClosedFile(lad1.name)

    # tests with unlocking gil
    createFileFromScratch(lad4.name, unlock_gil=True)
    readBackFirstFile(lad4.name, unlock_gil=True)
    createFileAsClone8(lad5.name, lad1.name, unlock_gil=True)

    if os.name == 'nt':
        # The script may crash with no message whatsoever
        # e.g. if mixing Debug and Release on windows.
        print("wrapper/test.py completed sucessfully")

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
