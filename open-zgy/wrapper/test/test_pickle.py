#!/usr/bin/env python3

"""
Test that types defined by this module and instances of them can be pickled.

  - Enums:
      - SampleDataType
      - UnitDimension
      - DecimationType
      - FinalizeAction

  - Named tuples:
      - Statistics
      - Histogram
      - FileStats

  - Exceptions
      - ZgyError
      - ZgyFormatError
      - ZgyCorruptedFile
      - ZgyUserError
      - ZgyInternalError
      - ZgyEndOfFile
      - ZgySegmentIsClosed
      - ZgyAborted
      - ZgyMissingFeature
      - ZgyIoError

  - Main types (don't test instances)
      - ZgyReader
      - ZgyWriter
      - ZgyUtils
"""

from enum import Enum
import pickle
import openzgycpp as zgy

nt_stats = zgy.SampleStatistics(999, 100.0, 1234.0, -100.0, +42.0)
nt_histo = zgy.SampleHistogram(16, 0.0, 1.0, list([1 for x in range(16)]))
nt_fstat = zgy.FileStats(*([0 for x in range(18)] + [1.0,False,(0,)]))
#print(repr(nt_stats))
#print(repr(nt_histo))
#print(repr(nt_fstat))

# List all custom Enum and NamedTuple types, both the type object itself
# and all known tags (for Enums) or example instances (for NamedTuple).
# Plus a few builtin types just to test the tester.

source = [
    str,
    int,
    float,

    "a string",
    42,
    3.1415826,

    zgy.SampleStatistics,
    zgy.SampleHistogram,
    zgy.FileStats,

    nt_stats,
    nt_histo,
    nt_fstat,

    zgy.SampleDataType,
    zgy.SampleDataType.unknown,
    zgy.SampleDataType.int8,
    zgy.SampleDataType.int16,
    zgy.SampleDataType.float,
    zgy.SampleDataType.float32,

    zgy.UnitDimension,
    zgy.UnitDimension.unknown,
    zgy.UnitDimension.time,
    zgy.UnitDimension.length,
    zgy.UnitDimension.arcangle,
    zgy.DecimationType.LowPass,

    zgy.DecimationType,
    zgy.DecimationType.WeightedAverage,
    zgy.DecimationType.Average,
    zgy.DecimationType.Median,
    zgy.DecimationType.Minimum,
    zgy.DecimationType.Maximum,
    zgy.DecimationType.Decimate,
    zgy.DecimationType.DecimateSkipNaN,
    zgy.DecimationType.AllZero,
    zgy.DecimationType.MostFrequent,
    zgy.DecimationType.MostFrequentNon0,
    zgy.DecimationType.AverageNon0,

    zgy.FinalizeAction,
    zgy.FinalizeAction.Delete,
    zgy.FinalizeAction.Keep,
    zgy.FinalizeAction.BuildIncremental,
    zgy.FinalizeAction.BuildFull,

    # Exceptions; only list the type object itself.
    # Instances are in source_ex.
    zgy.ZgyError,
    zgy.ZgyFormatError,
    zgy.ZgyCorruptedFile,
    zgy.ZgyUserError,
    zgy.ZgyInternalError,
    zgy.ZgyEndOfFile,
    zgy.ZgySegmentIsClosed,
    zgy.ZgyAborted,
    zgy.ZgyMissingFeature,
    zgy.ZgyIoError,

    # Complex types; only list the type object itself.
    zgy.ZgyReader,
    zgy.ZgyWriter,
    zgy.ZgyUtils,
]

source_ex = [
    zgy.ZgyError("x"),
    zgy.ZgyFormatError("x"),
    zgy.ZgyCorruptedFile("x"),
    zgy.ZgyUserError("x"),
    zgy.ZgyInternalError("x"),
    zgy.ZgyEndOfFile("x"),
    zgy.ZgySegmentIsClosed("x"),
    zgy.ZgyAborted("x"),
    zgy.ZgyMissingFeature("x"),
    zgy.ZgyIoError("x"),
]

# Create an enum the usual way.
class Color(Enum):
    RED = 1
    GREEN = 2
    BLUE = 5

# Create an enum using the functional API, as done in the module.
FColor = Enum(
    "FColor", [("RED",1),("GREEN",2),("BLUE",5)], module=__name__
)

def test_pickle_simple_enum():
    data = pickle.dumps([Color.BLUE])
    #print("Simple pickle", data)
    check = pickle.loads(data)
    #print("Simple round trip", check)
    assert check[0] == Color.BLUE

def test_pickle_functional_enum():
    data = pickle.dumps([FColor.BLUE])
    #print("Functional pickle", data)
    check = pickle.loads(data)
    #print("Functional round trip", check)
    assert check[0] == FColor.BLUE

def test_basic_enum_features():
    """
    Use one of the openzgycpp enums, verify that it works for
    converting enum->name, enum->value, name->enum, value->enum
    and is hashable and iterable.
    """
    assert zgy.SampleDataType.int16.name == "int16"
    assert zgy.SampleDataType.int16.value == 1002
    assert zgy.SampleDataType["int16"] == zgy.SampleDataType.int16
    assert zgy.SampleDataType(1002) == zgy.SampleDataType.int16
    assert zgy.SampleDataType.int16 in set([zgy.SampleDataType.int16])
    assert "int16" in list([ x.name for x in zgy.SampleDataType ])

def test_roundtrip(values, repr_only = False):
    values_repr = list([repr(x) for x in values])
    #print(*values_repr, sep="\n")
    data = pickle.dumps(values)
    #print(data)
    check = pickle.loads(data)
    check_repr = list([repr(x) for x in check])
    #print(*check_repr, sep="\n")
    assert check_repr == values_repr
    # Two Exception instances won't compare equal even when equivalent.
    # This is true both for built-in types and ZgyError etc.
    # So, need to relax the checking for those.
    if not repr_only:
        assert check == values

def declare_enum(T):
    """
    Purely for experiments, take an existing enum type and output code
    to create that same enum procedurally. When used with one of the
    enums defined in C++ the result should be the same but the enums
    won't be equal in the "is" sense.
    """
    a = [ '("{}",{})'.format(x.name, x.value) for x in T ]
    a = "[" + ", ".join(a) + "]"
    return 'my_{0} = Enum("{0}", {1}, module="{2}")'.format(
        T.__name__, a, T.__module__)

test_pickle_simple_enum()
test_pickle_functional_enum()
test_basic_enum_features()
test_roundtrip(source)
test_roundtrip(source_ex, True)

if False:
    print(declare_enum(zgy.SampleDataType))
    print(declare_enum(zgy.UnitDimension))
    print(declare_enum(zgy.DecimationType))
    print(declare_enum(zgy.FinalizeAction))
    print(declare_enum(Color))


