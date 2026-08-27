// Copyright 2017-2023, SLB
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

#pragma once

#include "../declspec.h"

#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <cstdint>
#include <functional>

namespace InternalZGY {
  struct LevelAndForce;
}

namespace InternalZGY {
#if 0
}
#endif

/**
 * Provide methods to query and update brick status.
 * This helps figure out when low resolution data should be computed.
 *
 * The class is NOT THREADSAFE. Adding locking inside the class
 * probably isn't helpful, because the caller will need to hold a lock
 * that covers both a call to getWork() and processing the contents.
 * Marking additional bricks as dirty while genlod is running is
 * probably a bad idea.
 *
 * Keep a per-brick flag for  physical status (Clean, FullyWritten, etc.).
 * This is the only mutable information in the class.
 *
 * All methods in the class deal with brick positions, not trace/sample.
 * The class doesn't even know what the current brick size is.
 *
 * Provide a virtual per-brick status computed from neighbors in the next
 * lower level. Virtual status can also be requested from a region such as
 * a brick columnn, with PartlyReady being returned if the individual
 * bricks don't have the same status.
 *
 * The bulk layer will need one instance per open file, and must
 * update status whenever data is written. At least for level 0;
 * higher levels might be updated here (to be decided).
 *
 * Provide methods to iterate over the survey or parts of it, with
 * the purpose of generating low resolution data and possibly
 * collect statistics and histogram. Only the iteration is handled
 * by this class.
 *
 * Each brick has a "physical" state:
 *  - Clean (never written).
 *  - PartlyWritten (caller knows there might be more).
 *  - FullyWritten (no more are expected).
 *  - WasWritten (Low resolution data has been generated).
 *    Apart from debugging, WasWritten is treated as Clean.
 *
 * And a "virtual state", for lod1+, computed from lod-1:
 *  - NotNeeded (none of the 8 lod-1 bricks are partly or fully written).
 *  - PartlyReady (some lod-1 written. Can generate lowres now if closing).
 *  - FullyReady (the low resolution brick can be written immediately).
 *
 * The first implementation used the same enum for both. That might have
 * caused some confusion. Because they aren't really the same.
 *
 * It is possible to enforce the "write only once" recommendation by
 * looking at illegal state transitions for the "physical" state and
 * illegal combinations of "physical" and "virtual". But note that
 * write once is only a recommendation. At least for on-prem.
 * Violating it will just cause a few extra passes of low resolution
 * compute. So, turn on enforcing only for testing.
 *
 * The physical state can go Clean->Partly->Partly->Fully->WasWritten,
 * or directly from Clean to Fully written. Most of the other
 * transitions imply that the brick has now been written more than
 * once. Without having the excuse of the user doing a r/m/w.
 * The remainder should not be possible. E.g. there is no way
 * of setting the "physical" back to Clean without resetting the
 * entire tracker.
 *
 * When the "virtual" state (which only applies to low resolution bricks)
 * is reported as Partly or Fully ready, the "physical" state of that
 * same brick is expected to be Clean.
 *
 * Because:
 *   - WasWritten clearly indicates that this brick will be written twice.
 *   - FullyWritten means the lowres hasn't been generated yet, but it
 *     would have been allowed. So there is still something wrong.
 *   - PartlyWritten should not happen at all for low resolution bricks
 *     due to how plan B and C is implemented. No r/m/w when generating
 *     low resolution.
 *   - Clean is thus the only valid state.
 */
class OPENZGY_TEST_API TrackTouched
{
public:
  enum class BrickPhysFlag : std::uint8_t
    {Clean=0,
     PartlyWritten,
     FullyWritten,
     WasWritten};
  enum class BrickVirtFlag : std::uint8_t
    {NotNeeded=0,
     PartlyReady,
     FullyReady};

  /**
   * A single I,J,K position.
   * Depending on context, this can contain brick numbers or trace / samples.
   * And the position can be either relative to lod0 (full resolution)
   * coordinates or relative to the current lod.
   */
  typedef std::array<std::int64_t,3> index3_t;

  /**
   * A list of I,J,K positions.
   */
  typedef std::vector<index3_t> indexlist_t;

  /**
   * I,J,K position(s), including the lod they belong to.
   *
   * Also indicate the vertical size at this level; the latter is a
   * convenience for when the position is actually 2d and the user
   * ignores pos[2]. If processing is per-brick instead of per
   * brick-column (I don't know whether both will be needed) then
   * zsize will always be 1.
   */
  typedef struct task {index3_t pos; int lod; std::int64_t zsize;} task_t;
  typedef std::vector<task_t> tasklist_t;

  /**
   * Information about a cube's size at each lod level. Indexed by lod.
   * For each lod it has the size in bricks. The type just happens
   * to be the same as an indexlist_t but is logically different.
   * The contents depend only on the survey size, so it will not change.
   */
  typedef std::vector<std::array<std::int64_t,3>> lodsizes_t;

  /**
   * Information about a cube's lookup table. Indexed by lod.
   * For each lod it has the offset in the lookup table of the first entry.
   * An extra element is pushed onto the end with the total slot count.
   * so it is valid to access offsets_[nlods].
   * The contents depend only on the survey size, so it will not change.
   */
  typedef std::vector<std::int64_t> offsets_t;

  /**
   * Logger for debugging.
   */
  typedef std::function<bool(int, const std::string&)> logger_t;

private:
  // Information about the file that doesn not change after it is created.
  // lodsizes_ and brickoffsets_ are computed from the survey size.
  // If called from inside a ZgyInternalBulk instance, use the following:
  //     this->_metadata->ih().lodsizes();
  //     this->_metadata->ih().brickoffsets();
  const lodsizes_t lodsizes_;
  const offsets_t  brickoffsets_;
  const logger_t   loggerfn_;
  // Status of each full- and low resolution brick.
  std::vector<BrickPhysFlag> modified_bricks_;

public:
  TrackTouched(
       const lodsizes_t& lodsizes,
       const offsets_t&  brickoffsets,
       const logger_t& logger);
  TrackTouched(const TrackTouched&) = delete;
  TrackTouched(const TrackTouched&&) = delete;
  TrackTouched& operator=(const TrackTouched&) = delete;
  TrackTouched& operator=(const TrackTouched&&) = delete;
  TrackTouched(const TrackTouched&, bool dummy);

public:
  void setAllClean();
  bool isAllClean() const;
  void set1Written(const index3_t& brickpos, std::int32_t lod, bool rmw);
  void setCWritten(const indexlist_t& brickposlist, std::int32_t lod, bool rmw);
  void setWritten(const indexlist_t& brickposlist, std::int32_t lod, bool rmw);
  static std::int64_t countIO(const task_t& task);
  tasklist_t getWork(const LevelAndForce& lodmode, const indexlist_t& roi) const;
  tasklist_t getWorkAndClear(const LevelAndForce& lodmode, const indexlist_t& roi);
  void showDirty(bool lowlevel, int loglevel) const;

private:
  static bool doThisNow(BrickVirtFlag flag, int level, int force, int maxlevel);
  static indexlist_t removeDuplicates(const indexlist_t&);
  static indexlist_t getLowResBrickColumns(const indexlist_t&, bool column);
  void setPhysByLinearIndex(std::int64_t index, BrickPhysFlag flag);
  BrickVirtFlag getVirtByBrickPosition(const index3_t& brickpos, std::int32_t lod, bool rmw) const;
  BrickPhysFlag getPhysByLinearIndex(std::int64_t index, bool ignore_was_written) const;

  indexlist_t listOneLevel(int lod, int force, const indexlist_t& roi) const;
  tasklist_t listUpToLevel(const LevelAndForce& lodmode, const indexlist_t& roi);

  std::int64_t getBrickLookupIndex(index3_t pos, std::int64_t lod) const;

  std::vector<std::int64_t> getParentBrickIndices(
      std::int64_t i0, std::int64_t j0, std::int64_t k0,
      std::int64_t ni, std::int64_t nj, std::int64_t nk,
      std::int64_t lod) const;

  bool logger_(int priority, const std::string& ss = std::string()) const;
  bool logger_(int priority, const std::ios& ss) const;
};

} // namespace
