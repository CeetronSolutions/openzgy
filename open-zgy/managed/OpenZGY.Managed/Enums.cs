// Copyright 2017-2022, Schlumberger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

namespace OpenZGY.Managed
{
  public enum SampleDataType
  {
    // must be kept in sync with the declaration in native/src/api.h
    unknown = 1000,
    int8 = 1001,
    int16 = 1002,
    float32 = 1003,
  }

  public enum UnitDimension
  {
    // must be kept in sync with the declaration in native/src/api.h
    unknown = 2000,
    time = 2001,
    length = 2002,
    arcangle = 2003,
  }

  public enum DecimationType
  {
    // must be kept in sync with the declaration in native/src/api.h
    LowPass = 100, ///< \brief Lowpass Z / decimate XY
    WeightedAverage = 101, ///< \brief Weighted averaging (depends on global stats)
    Average = 102, ///< \brief Simple averaging
    Median = 103, ///< \brief Somewhat more expensive averaging
    Minimum = 104, ///< \brief Minimum value
    Maximum = 105, ///< \brief Maximum value
    //MinMax         = 106, ///< \brief Checkerboard of minimum and maximum values
    Decimate = 107, ///< \brief Simple decimation, use first sample
    DecimateSkipNaN = 108, ///< \brief Use first sample that is not NaN
    //DecimateRandom = 109, ///< \brief Random decimation using a fixed seed
    AllZero = 110, ///< \brief Just fill the LOD brick with zeroes
    //WhiteNoise     = 111, ///< \brief Fill with white noise
    MostFrequent = 112, ///< \brief The value that occurs most frequently
    MostFrequentNon0 = 113, ///< \brief The non-zero value that occurs most frequently
    AverageNon0 = 114, ///< \brief Average value, but treat 0 as NaN.
  }

  public enum FinalizeAction
  {
    // must be kept in sync with the declaration in native/src/api.h
    /**
     * Remove any existing information. If force=true this is unconditinal.
     * Otherwise it changes to "Keep" if the information already exists
     * and is not stale. With force=off, Delete and Keep actually do the same.
     */
    Delete = 3001,

    /**
     * Keep any existing information. If force=true this will even keep the
     * information if it is looks stale. Otherwise the mode changes to "Delete"
     * if OpenZGY cannot guarantee that the information is still correct.
     */
    Keep = 3002,

    /**
     * Do the minimum amount of work needed to bring the derived data up to
     * date. May cause the statistics and histogram to be less accurate.
     * Changes to "Keep" if the information is already up to date. changes
     * to "BuildFull" if information does not exist already or if the
     * file is compressed or if incremental cannot work for other reasons.
     * Setting force=true might do a little more work, to be defimed.
     */
    BuildIncremental = 3003,

    /**
     * Delete and re-create derived information from scratch. If force=true
     * this is unconditional. That may be useful if you want to change the
     * decimation algorithm. Otherwise change to "Keep" if the information
     * is already up to date. Never changes to BuildIncremental.
     */
    BuildFull = 3004,

    /**
     * As BuildFull, but do not collect or store histogram and statistics.
     * This might be useful to speed up writes that we know will only be
     * read by Petrel. Or other apps that don't strictly need them.
     * CAVEAT: This is incompatible with DecimationType::WeightedAverage.
     * So, indirectly this option will lower the quality of LoD data,
     */
    BuildNoHistogram = 3005,

    /**
     * Value used as default argument in finalize, inside close, etc.
     */
    BuildDefault = BuildFull,
  }
}
