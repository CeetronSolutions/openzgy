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

#include "tracktouched.h"
#include "lookuptable.h"

#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <cstdint>
#include <functional>
#include <algorithm>
#ifndef _WIN32
#include <unistd.h>
#endif

namespace InternalZGY {
#if 0
}
#endif

TrackTouched::TrackTouched(
     const lodsizes_t& lodsizes,
     const offsets_t&  brickoffsets,
     const logger_t&   logger)
  : lodsizes_(lodsizes)
  , brickoffsets_(brickoffsets)
  , loggerfn_(logger)
  , modified_bricks_()
{
  logger_(10, std::stringstream().flush() << "TrackTouched created");
}

/**
 * The copy constructor has been made explicit by adding a dummy
 * parameter. Copy or assign by accident can easily create bugs.
 */
TrackTouched::TrackTouched(const TrackTouched& other, bool /*dummy*/)
  : lodsizes_(other.lodsizes_)
  , brickoffsets_(other.brickoffsets_)
  , loggerfn_(other.loggerfn_)
  , modified_bricks_(other.modified_bricks_)
{
  logger_(10, std::stringstream().flush() << "TrackTouched copied");
}

/**
 * Mark all bricks as clean, in the sense that no low resolution data
 * needs to be generated. Call this after writing the entire survey
 * with a constant value, followed by a trivial finalize that just
 * sets all the low resolution bricks to the same value.
 */
void
TrackTouched::setAllClean()
{
  modified_bricks_.clear();
}

/**
 * Return true if setAllClean has been called.
 * Return false if not, but it is still possible that everything is clean.
 */
bool
TrackTouched::isAllClean() const
{
  return modified_bricks_.empty();
}

/**
 * \brief Mark a brick as having been written.
 *
 * \details
 *
 * At lod 0, just this brick gets flagged. Usually as FullyWritten,
 * but if the caller did a read/modify/write then we must expect
 * additional writes. So the brick is then flagged PartlyWritten.
 * With no attempt to detect when there has beem enough r/m/w to
 * fill the entire brick.
 *
 * At higher lods, this brick gets flagged the same way. In addion,
 * the (usually 8) source brics get changed back to WasWritten (which
 * is mostly treated as Clean) since there is no longer a need to
 * generate this particular brick.
 *
 * The positions are in trace / sample at the specified lod.
 */
void
TrackTouched::set1Written(
     const index3_t& brickpos,
     std::int32_t lod,
     bool maybe_more)
{
  const std::int64_t bix = getBrickLookupIndex(brickpos, lod);
  setPhysByLinearIndex(bix, maybe_more ?
                       BrickPhysFlag::PartlyWritten :
                       BrickPhysFlag::FullyWritten);
  logger_(10, std::stringstream().flush()
          << "set1Written(lod=" << lod << ", linear=" << bix << ")");
  if (lod > 0) {
    std::vector<std::int64_t> parents =
      getParentBrickIndices
      (brickpos[0], brickpos[1], brickpos[2], 1, 1, 1, lod);
    for (std::int64_t parent : parents)
      setPhysByLinearIndex(parent, BrickPhysFlag::WasWritten);
  }
}

/**
 * Mark all bricks in one brick-column as written.
 * pos[2] is ignored.
 * See set1Written() for details.
 */
void
TrackTouched::setCWritten(
     const indexlist_t& brickposlist,
     std::int32_t lod,
     bool maybe_more)
{
  const std::int64_t vsize = this->lodsizes_[lod][2];
  for (const index3_t& brickpos : brickposlist)
    for (std::int64_t kk = 0; kk < vsize; ++kk)
      set1Written(index3_t{brickpos[0], brickpos[1], kk}, lod, maybe_more);
}

/**
 * Mark multiple bricks as being written.
 * See set1Written() for details.
 */
void
TrackTouched::setWritten(
     const indexlist_t& brickposlist,
     std::int32_t lod,
     bool maybe_more)
{
  for (const index3_t& pos : brickposlist)
    set1Written(pos, lod, maybe_more);
}

/**
 * \brief See if the specified brick or brick-column of lowres data is
 * ready to compute, based on the status of its parents at lod-1.
 *
 * \details
 *
 * The input position is in trace / sample numbers relative to lod.
 * Calling this with lod=0 does not make sense, as that lod is not computed.
 *
 * The return value will be one of the following:
 *
 *  - Clean: All parents are unchanged. No need to generate unless
 *    the application asked for a full rebuild.
 *    Lod and/or position outside the survey are also assumed to be clean.
 *  - PartlyWritten: Only some parents ready. Or some (or all) flagged Partial
 *    because a read/modify/write was done. There might be more writes.
 *    Defer generate until status changes to FullyWritten or the file
 *    is being closed.
 *  - FullyWritten: All parents written. Ready to generate now,
 *    but application might still wish to defer everything to the end.
 *
 * Note that the fuction only looks one level down. Caller must check level
 * 1 first, and generate bricks at that level before it can checks level 2.
 * Same for the other levels.
 *
 * TODO-WIP-BrickedAPI: Caveat: Files opened for update are assumed to be
 * all clean. But, depending on the parameters to the previous finalize()
 * it might not be. At the next format change a "lod dirty" flag should
 * be added. That also signifies that statistics and histogram might
 * be wrong.
 */
TrackTouched::BrickVirtFlag
TrackTouched::getVirtByBrickPosition(
     const index3_t& brickpos,
     std::int32_t lod,
     bool entire_column) const
{
  const std::int64_t vsize     = this->lodsizes_[lod][2];
  const std::int64_t vbeg      = entire_column ? 0 : brickpos[2];
  const std::int64_t vcount    = entire_column ? vsize : 1;

  std::vector<std::int64_t> parents =
    getParentBrickIndices
    (brickpos[0], brickpos[1], vbeg, 1, 1, vcount, lod);
  if (!parents.empty()) {
    const BrickPhysFlag first = getPhysByLinearIndex(parents.front(), true);
    for (std::int64_t parent : parents)
      if (getPhysByLinearIndex(parent, true) != first)
        return BrickVirtFlag::PartlyReady;
    // If we get here, all parents have the same phys flag.
    // Clean and WasWritten count as identical in this case.
    switch (first) {
    case BrickPhysFlag::Clean:
    case BrickPhysFlag::WasWritten:
      return BrickVirtFlag::NotNeeded;
    case BrickPhysFlag::PartlyWritten:
    default:
      return BrickVirtFlag::PartlyReady;
    case BrickPhysFlag::FullyWritten:
      return BrickVirtFlag::FullyReady;
    }
  }
  else {
    return BrickVirtFlag::NotNeeded; // All outside range.
  }

  // TODO-WIP-BrickedAPI: Enforce that bricks can only be written,
  // once, throwing an exception or invoking a callback
  // as soon as possible. Maybe long before the
  // actual write is done. At this point, if the virtual state
  // is partly or fully written then our own physical flag
  // should not be FullyWritten or WasWritten. Because that
  // would imply that we dodged a bullet (FullyWritten)
  // or will actually do a duplicate write some time later
  // (WasWritten). It gets a bit messy if entire_column is true.
  // Maybe not check here, and just wait for setPhysByLinearIndex to
  // to catch it.
}

/**
 * Choose whether a particular brick or region should be generated now,
 * based on how the "maxlevel" and "force" parameters.
 * "force" is currently an integer but might change to an enum later.
 *
 *   0 - No-op. Don't generate any lods at this point.
 *
 *   1 - Only generate lods if all 8 parent bricks (or fewer at survey
 *       edge) have been modified. The assumption is that the
 *       application will write each brick at most once. So, the
 *       application probably won't be writing these 8 bricks again,
 *       and we might as well finish the lod generation for this area.
 *       It is appropriate but not required to do this after every
 *       single write.
 *
 *   2 - Generate lods if any of the 8 parent bricks have been nodified.
 *       This is appropriate when flushing the file.
 *
 *   3 - Generate lods unconditionally. Use this to fix inconsistencies
 *       or to switch to using a different LOD generation algorithm.
 *
 * "maxlevel" has the following meaning:
 *
 *   -1: Write all LOD levels as soon as "force" lets us do so.
 *
 *    0: do not write anything here, wait for explicit flush.
 *
 *    1: write LOD1 bricks as soon as we have data, LOD2 and higher later.
 *       This avoids the issue of LOD2+ depending on global statistics.
 *
 *    2: Similarly, it is possible to write LOD1+LOD2 now. Et cetera.
 *
 * Low resolution bricks can be generated both after every write and in
 * finalize() / close() of the file. In the first case, use force=1
 * and any maxlevel. Including 0 to disable this feature. In the second
 * case use force 2 or 3, and maxlevel -1.
 */
bool /*static*/
TrackTouched::doThisNow(BrickVirtFlag flag, int level, int force, int maxlevel)
{
  if (maxlevel >= 0 && level > maxlevel)
    return false;
  switch (force) {
  default:
  case 0: return false;                               // Never
  case 1: return (flag == BrickVirtFlag::FullyReady); // All parents dirty
  case 2: return (flag != BrickVirtFlag::NotNeeded);  // Any parent dirty
  case 3: return true;                                // Always
  }
}

/**
 * Remove duplicates from a vector of ijk positions.
 */
TrackTouched::indexlist_t /*static*/
TrackTouched::removeDuplicates(const indexlist_t& list)
{
  indexlist_t copy(list);
  if (list.size() > 1) {
    std::sort(copy.begin(), copy.end());
    copy.erase(std::unique(copy.begin(), copy.end()), copy.end());
  }
  return copy;
}

/**
 * Given a vector of ijk positions, return these positions in the next
 * higher level. Remove duplicates, so passing 8 adjacent fully aligned
 * cubes may result in just one output position.
 *
 * If entire_column = true, the vertical position is ignored. So all
 * bricks in the same brick column will be duplicates of the first.
 *
 * The input uses brick (not sample) numbers relative to lod, and the
 * output is brick numbers relative to lod+1. At this low level
 * it would also have worked with trace / sample positons if they
 * are all aligned the same way.
 */
TrackTouched::indexlist_t /*static*/
TrackTouched::getLowResBrickColumns(const indexlist_t& list, bool entire_column)
{
  std::vector<std::array<std::int64_t,3>> result;
  for (const auto& pos : list)
    result.push_back(std::array<std::int64_t,3>
                     {pos[0]/2, pos[1]/2, entire_column ? 0 : pos[2]/2});
  return removeDuplicates(result);
}

/**
 * \brief Set the state of a single brick only.
 *
 * \details
 *
 * The brick's position should be already converted to a linear index
 * into the modified_bricks_ vector.
 *
 * \internal
 *
 * If nothing has been touched yet, and caller wants to set a brick as
 * clean, this is a no-op and we don't need to allocate space for the
 * lookup if it hasn't been done already.
 */
void
TrackTouched::setPhysByLinearIndex(std::int64_t index, BrickPhysFlag flag)
{
  BrickPhysFlag oldflag = BrickPhysFlag::Clean;
  if (!this->modified_bricks_.empty() || flag != BrickPhysFlag::Clean) {
    const std::int64_t max_index = this->brickoffsets_.back();
    if (index >= 0 && index < max_index) {
      if ((std::int64_t)this->modified_bricks_.size() < max_index)
        this->modified_bricks_.resize(max_index, BrickPhysFlag::Clean);
      oldflag = this->modified_bricks_[index];
      this->modified_bricks_[index] = flag;
    }
  }

  // Enforce that bricks can only be written once,
  // throwing an exception or invoking a callback
  // as soon as possible. Maybe long before the
  // actual write is done.
  // TODO-WIP-BrickedAPI: Don't throw exceptions in production.
  // A polite message in the log window might be acceptable.
  if ((int)flag < (int)oldflag ||
    (flag == BrickPhysFlag::FullyWritten && oldflag == flag))
  {
    // Some unit tests violate this recommendation on purpose.
    // There needs to be some way for the application to
    // control whether this is enforced.
    //throw std::runtime_error("Illegal phys state change. Bricks will be written more than once.");
  }
}

/**
 * \brief Get the state of a single brick only.
 *
 * \details
 *
 * The brick's position should be already converted to a linear index
 * into the modified_bricks_ vector.
 *
 * Any brick outside the valid range is assumed to be clean.
 *
 */
TrackTouched::BrickPhysFlag
TrackTouched::getPhysByLinearIndex(std::int64_t index, bool ignore_was_written) const
{
  if (index < 0 || index >= (std::int64_t)this->modified_bricks_.size())
    return BrickPhysFlag::Clean;
  else if (ignore_was_written && this->modified_bricks_[index] == BrickPhysFlag::WasWritten)
    return BrickPhysFlag::Clean;
  else
    return this->modified_bricks_[index];
}

/**
 * \brief
 *
 * List lods to be generated for a single level, possibly limited to
 * a list of bricks that have recently been modified.
 *
 * \details
 *
 * There is an assumption that the lod-1 layer has already been processed
 * and the dirty bits have been updated accordingly. If looping over all
 * lods, the caller needs to immediately call setCWritten() on the result.
 * This means that the state gets updated before the actual processing starts.
 * If a dry run is needed then a temporary will be needed.
 *
 * roi should contain brick numbers, not trace / sample numbers.
 * If roi is empty, process the entire survey.
 *
 * The ROI list refers to lod-1, i.e. the *source* of the lowres
 * bricks to be written. The first call to this function will
 * typically have lod=1 and a list of positions in lod0 (if
 * incremental) or an empty roi (if called from finalize).
 *
 * On the other hand, the positions in the "children" array,
 * which is what the generator sees, refer to lod.
 *
 * The source cubes must be already tagged with the correct state.
 *
 * Force describes how eager we are to generate the lods,
 * see GenerateLOD_brick for details.
 *
 * The return value is a list of bricks, relative to lod, that need
 * to be (a) generated and (b) used as the roi for further decimation.
 * Note that listing the next level up can only be done after the
 * bricks at this level have been generated. Or appear to be generated
 * by calling setCWritten(). It doesn't make sense to run the algorithm
 * on higher lod levels until whatever caused a brick to not be
 * generated has been fixed.
 *
 * \internal
 *
 * - Hard code the choice to do whole brick columns, and never r/m/w.
 * - Loop over roi or entire survey, generate brick-columns.
 *   Remove duplicates so there is only one call per brick-column.
 * - Optionally mark the low resolution bricks as written, and the
 *   source as clean. Note that clear=true is required if called from
 *   listUpToLevel() that returns ready bricks in all levels.
 * - Return the list of low resolution bricks that can be decimated
 *   further, subject to the same checks done at this level.
 *
 * Some alternatives, I haven't decided yet which to use.
 *
 * 1) generateOneLevel() invokes the actual call to decimate by calling
 *    a provided functor. One brick-column at a time. generateOneLevel()
 *    will also clear the dirty bits. Whatever is inside the functor
 *    will probably do the same. Doing it twice won't hurt.
 *    When finalizing, the total count is needed up front for the
 *    progress callback. Do this using a dry run in generateUpToLevel()
 *    where the functor just tallies the request. This must be done
 *    in a TEMPORARY CLONE.
 *
 * 2) [CURRENT]
 *    Remove the functor. Rename generate* to list*. listUpToLevel()
 *    is changed to return a list of bricks to be generated, including
 *    lod number, bricks in column, ipos, jpos. Caller then needs to
 *    compute and store each of these, in the order listed. The order
 *    listed will be to group by lod. Starting at the lower level and
 *    completing that before moving on to the next. The function will
 *    reset dirty bits. So the caller might need to run in a
 *    temporary; this depends on the logic inside the caller.
 *
 * 3) As (2), but do higher levels as soon as possible. Increasing the
 *    chance that the decimator will find data in the cache. But, the
 *    WeightedAverage algorithm might produce even less reproducible
 *    results.
 *
 * 4) Remove the functor. Remove generateUpToLevel() and implement
 *    that logic inside the caller instead. generateOneLevel() then
 *    doesn't need to clear the dirty bits up front. This is cleaner
 *    in the sense that dirty bits are cleared immediately after the
 *    data is written. It is less clean in the sense that some
 *    internal logic is moved out of the TrackTouched class.
 *    That logic is currently so subtle that this is a really bad idea.
 */
TrackTouched::indexlist_t
TrackTouched::listOneLevel(int lod, int force, const indexlist_t& roi) const
{
  indexlist_t children;
  indexlist_t next;
  if (roi.empty()) {
    const index3_t size = this->lodsizes_[lod];
    index3_t pos{0,0,0};
    for (pos[0]= 0; pos[0] < size[0]; ++pos[0]) {
      for (pos[1]= 0; pos[1] < size[1]; ++pos[1]) {
        children.push_back(pos);
      }
    }
  }
  else {
    children = getLowResBrickColumns(roi, /*columns=*/true);
  }
  for (const index3_t& pos : children) {
    BrickVirtFlag flag = getVirtByBrickPosition(pos, lod, /*columns=*/true);
    if (doThisNow(flag, lod, force, -1)) {
      next.push_back(pos);
    }
  }
  return next;
}

/**
 * \brief
 * List LODs ready to compute from level 1 up to and including "maxlevel".
 *
 * \details
 * See doThisNow() for an explanation of maxlevel and force.
 * If "roi" is empty, the entire survey is processed. Otherwise only
 * lowres data based on the indicated bricks, in LOD0 coordinates.
 *
 * The dirty bricks will be cleared as if the bricks
 * actually got generated.
 *
 * TODO-WIP-BrickedAPI: remember to test corner cases, surveys with
 * 0 and 1 lod levels.
 *
 * The function is NOT THREAD SAFE. The class itself can be made
 * thread safe by protecting its own data with a lock. But, when
 * called from genlod, a lock should be held until both obtaining
 * the list and processing it has been completed. Or just make sure
 * the processing is already single threaded at this point.
 */
TrackTouched::tasklist_t
TrackTouched::listUpToLevel(const LevelAndForce& lodmode, const indexlist_t& roi)
{
  int maxlevel = lodmode.level;
  int force = lodmode.force;
  tasklist_t result;
  const std::int32_t nlods = (std::int32_t)this->lodsizes_.size();
  if (maxlevel == 0 || force == 0 || (force != 3 && isAllClean()) || nlods == 1)
    return result;

  // brickoffsets_ has size nlods+1, prefer to use lodsizes_
  // but this is still an implementation detail.
  if (maxlevel < 0 || maxlevel >= nlods)
    maxlevel = nlods - 1;

  // Ordering of the loop is vital, as (re)generation of LOD N+1 depends
  // on dirty range for LOD N, so LOD N must be done before LOD N+1.

  if (force < 2) {
    // TODO-WIP-BrickedAPI: TEST:
    // This code might miss some bricks, but those should be caught by finalize.
    // The code assumes that if no work was found in level N, then there is
    // no need to check higher levels. More generaly, int assumes that when
    // the last brick of 8 parents is ready then the low resolution verisons
    // will be computed immediately. This might not be true if maxlevel was set.
    indexlist_t next = listOneLevel(/*lod=*/1, force, roi);
    setCWritten(next, /*lod=*/1, /*maybe_more=*/false);
    for (const index3_t& pos : next)
      result.push_back(task_t{pos, /*lod=*/1, this->lodsizes_[/*lod=*/1][2]});
    // At higher levels, an empty roi should do nothing, not do all.
    for (int lod = 2; lod <= maxlevel && !next.empty(); ++lod) {
      next = listOneLevel(lod, force, next);
      setCWritten(next, lod, /*maybe_more=*/false);
      for (const index3_t& pos : next)
        result.push_back(task_t{pos, lod, this->lodsizes_[lod][2]});
    }
  }
  else {
    // More thorough check, usually done in finalize with an empty ROI.
    // Could have used this branch unconditionally, at some extra cost
    // testing all lowres bricks above us in spite of them almost
    // certainly being not ready.
    indexlist_t roi_at_lod(roi);
    for (int lod = 1; lod <= maxlevel; ++lod) {
      indexlist_t next = listOneLevel(lod, force, roi_at_lod);
      logger_(2, "\nScanning lod " +
              std::to_string(lod) +
              " got " +
              std::to_string(next.size()) +
              " dirty");
      setCWritten(next, lod, /*maybe_more=*/false);
      for (const index3_t& pos : next)
        result.push_back(task_t{pos, lod, this->lodsizes_[lod][2]});
      // Adjust for the next (more decimated) level.
      // Note: When computing the size of a more decimated level, round up.
      // When finding a brick position in a more decimated level, round down
      for (index3_t& pos : roi_at_lod)
        for (int dim = 0; dim < 3; ++dim)
          pos[dim] /= 2;
      roi_at_lod = removeDuplicates(roi_at_lod);
    }
  }
  return result;
}

/**
 * In plan B, one task represents one brick-column to be written at
 * the current lod, and 4 brick-columns of twice the height to be
 * read. Reads outside the survey (because of odd sizes) are also
 * counted. The number of bricks in a brick-column at this lod is
 * helpfully provided in the task_t instance.
 *
 * The idea is to get as accurate progress information as possible.
 * A simpler approach would be to report brick-columnbs instead of
 * individual brics.
 */
std::int64_t
TrackTouched::countIO(const task_t& task)
{
  return 9 * task.zsize;
}

/**
 * \brief
 * Get a list of everything that needs to be generated.
 *
 * \details
 * This is a public interface. The method updates the state.
 * The bricks will prematurely be marked as done. This is more
 * efficient, but might not be acceptable to the caller.
 */
TrackTouched::tasklist_t
TrackTouched::getWorkAndClear(const LevelAndForce& lodmode, const indexlist_t& roi)
{
  //showDirty(0, false);
  //TrackTouched::tasklist_t result = listUpToLevel(maxlevel, force, roi);
  //logger_(0, "getWorkAndClear: " + std::to_string(result.size()));
  //return result;
  // Could lock a mutex here, and hold it only while cloning.
  return listUpToLevel(lodmode, roi);
}

/**
 * \brief
 * Get a list of everything that needs to be generated.
 *
 * \details
 * This is a public interface.
 */
TrackTouched::tasklist_t
TrackTouched::getWork(const LevelAndForce& lodmode,const indexlist_t& roi) const
{
  //showDirty(0, false);
  const std::int32_t nlods = (std::int32_t)this->lodsizes_.size();
  if (lodmode.level == 0 || lodmode.force == 0 || (lodmode.force != 3 && isAllClean()) || nlods == 1)
    return tasklist_t{};

#if 0
  if (false && lodmode.level == 1) {
    // TODO-WIP-BrickedAPI: Not sure I want to enable (and test!) this tweak.
    // With just a single layer to process, major shortcuts are possble
    // and the clone isn't needed.
    indexlist_t next = listOneLevel(/*lod=*/1, lodmode.force, roi);
    tasklist_t result;
    for (const index3_t& pos : next)
      result.push_back(task_t{pos, /*lod=*/1, this->lodsizes_[/*lod=*/1][2]});
    return result;
  }
#endif

  // Use a clone because listUpToLevel() clears the flags as a side effect.
  std::shared_ptr<TrackTouched> clone;
  {
    // Could lock a mutex here, and hold it only while cloning.
    clone = std::make_shared<TrackTouched>(*this, /*dummy=*/true);
  }
  //TrackTouched::tasklist_t result = clone->listUpToLevel(maxlevel,force,roi);
  //logger_(0, "getWorkAndClear: " + std::to_string(result.size()));
  //return result;
  return clone->listUpToLevel(lodmode, roi);
}

/**
 * \brief For debugging, show dirty bricks.
 *
 * \details
 * The ascii art display shows brick-columns (combines all bricks in Z)
 * while the numbers printed after each LOD are per brick.
 */
void
TrackTouched::showDirty(bool lowlevel, int loglevel) const
{
  if (!logger_(loglevel))
    return;

#ifndef _WIN32
  const bool color = ::isatty(1);
#else
  const bool color = false;
#endif

  static auto flagname =
    [](BrickPhysFlag flag) -> std::string
    {
      switch (flag) {
      case BrickPhysFlag::Clean: return "Clean";
      case BrickPhysFlag::PartlyWritten: return "PartlyWritten";
      case BrickPhysFlag::FullyWritten: return "FullyWritten";
      case BrickPhysFlag::WasWritten: return "WasWritten";
      default: return std::to_string(int(flag));
      }
    };
  if (isAllClean()) {
    logger_(loglevel, "All bricks are clean");
  }
  else if (std::is_sorted(modified_bricks_.begin(), modified_bricks_.end()) &&
      modified_bricks_.front() == modified_bricks_.back()) {
    logger_(loglevel, std::stringstream().flush()
            << "All " << modified_bricks_.size() << " bricks are "
            << flagname(modified_bricks_.front()));
  }
  else {
    logger_(loglevel, lowlevel ?
            "Dumping dirty bricks":
            "Dumping bricks ready to generate (looks just 1 level down)");
    for (std::int32_t lod = (std::int32_t)lodsizes_.size() - 1; lod >= 0; --lod) {
      logger_(loglevel, "LOD " + std::to_string(lod));
      //int num_clean{0}, num_partly{0}, num_fully{0};
      for (std::int64_t jj=0; jj < lodsizes_[lod][1]; ++jj) {
        std::stringstream line;
        for (std::int64_t ii = 0; ii < lodsizes_[lod][0]; ++ii) {
          std::int64_t index0 = getBrickLookupIndex(index3_t{ii, jj, 0}, lod);
          BrickPhysFlag dirty = getPhysByLinearIndex(index0, false);
          BrickVirtFlag state = getVirtByBrickPosition(index3_t{ii ,jj,0 }, lod, /*col=*/true);
          const char *pflag = "?";
          const char *vflag = "?";
          if (index0 < 0 && state == BrickVirtFlag::NotNeeded)
            pflag = ">"; // Outside
          else if (index0 < 0)
            pflag = "?"; // Error, outside should be clean.
          else {
            switch (dirty) {
            case BrickPhysFlag::Clean:         pflag = "."; break;
            case BrickPhysFlag::PartlyWritten: pflag = "o"; break;
            case BrickPhysFlag::FullyWritten:  pflag = "*"; break;
            case BrickPhysFlag::WasWritten:    pflag = "-"; break;
            default:                           pflag = "?"; break;
            }
          }
          switch (state) {
          case BrickVirtFlag::NotNeeded:   vflag = "."; break;
          case BrickVirtFlag::PartlyReady: vflag = "o"; break;
          case BrickVirtFlag::FullyReady:  vflag = "*"; break;
          default:                         vflag = "?"; break;
          }
          if (!color) {
            line << (lowlevel ? pflag : vflag);
          }
          else {
            switch (state) {
            case BrickVirtFlag::NotNeeded:   vflag = "\033[47;1m"; break;
            case BrickVirtFlag::PartlyReady: vflag = "\033[43;1m"; break;
            case BrickVirtFlag::FullyReady:  vflag = "\033[42;1m"; break;
            default:                         vflag = ""; break;
            }
            line << vflag << pflag << "\033[0m";
          }
        }
        logger_(loglevel, line.str());    // After each inline
      }
      // After each LOD
      //logger_(loglevel, std::stringstream().flush()
      //        << "Total " << num_clean << " clean, "
      //        << num_partly << " partly, and "
      //        << num_fully << " written at lod " << lod << ".");
    }
  }
}

/**
 * \brief
 * Convert from brick number to linear index.
 * The brick number is relative to the specified lod.
 *
 * \details
 * The mapping has been chosen to be the same as the main lookup
 * table, so this call can be forwarded to class LookupTable even
 * though the two aren't related.
 */
std::int64_t
TrackTouched::getBrickLookupIndex(index3_t pos, std::int64_t lod) const
{
  return LookupTable::getBrickLookupIndex(
    pos[0], pos[1], pos[2], lod,
    this->lodsizes_, this->brickoffsets_,
    /*throw_on_error=*/false);
}

/**
 * \brief Get all linear indices for the parents (i.e. lod-1) of the
 * specified bricks. Input indices are block numbers relative to the
 * given lod, NOT lod-1.
 *
 * \details
 * Brick outside range, whether caused by bad input parameters or by
 * trying to get 8 parents close to the survey edge, are quietly ignored.
 *
 * To help code re-use, the function has been moved to class LookupTable
 * even though the two classes aren't related.
 */
std::vector<std::int64_t>
TrackTouched::getParentBrickIndices(
     std::int64_t i0, std::int64_t j0, std::int64_t k0,
     std::int64_t ni, std::int64_t nj, std::int64_t nk,
     std::int64_t lod) const
{
  return LookupTable::getParentBrickIndices
    (i0, j0, k0, ni, nj, nk, lod,
     this->lodsizes_, this->brickoffsets_,
     /*throw_on_error=*/false);
}

/**
 * Duplicated from bulk.cpp, see there for details.
 */
bool
TrackTouched::logger_(int priority, const std::string& message) const
{
  return loggerfn_(priority, message);
}

/**
 * Duplicated from bulk.cpp, see there for details.
 */
bool
TrackTouched::logger_(int priority, const std::ios& ss) const
{
  auto sstream = dynamic_cast<const std::stringstream*>(&ss);
  return logger_(priority, sstream ? sstream->str() : std::string());
}

} // namespace

// TODO-WIP-BrickedAPI: The buffers used by genlod need to be kept
// somewhere to avoid lots of malloc / free. Not sure whether they
// belong in ZgyInternalBulk, GenLodB, or in TrackTouched.
//
// TODO-WIP-BrickedAPI: Make the class more self contained
// by duplicating and moving code from LookupTable.
//
// TODO-WIP-BrickedAPI: Nice to have: keep track of the number
// of times each brick has been written. Warn if more than one.
// Might not be needed if validating state changes.
//
// TODO-WIP-BrickedAPI: As a consistency check, after processing the
// entire file in finalize, check that all physical flags are Clean or
// WasWritten before (maybe) resetting the entire tracker.
