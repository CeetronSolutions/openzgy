// Copyright 2017-2023, Schlumberger
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

#include "readwriterelay.h"

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

namespace OpenZGY { namespace Tools {
#if 0
}}
#endif

/**
 * Create a large virtual survey by logically mirroring the source
 * survey a specified number of times.
 *
 * Intercept all calls to the public API that needs to return a
 * different result compared to the input.
 *
 * This is meant for a tool that artificially inflates a ZGY file by
 * mirroring it. Exclusively for testing.
 *
 * The idea is that the mirroring feature (and possibly later the
 * upsampling feature) can be added to the existing zgycopyc tool with
 * just a few lines of code. So it won't bloat that tool even more.
 *
 * There are both benefits and drawbacks to this strategy compared to
 * having a separate tool. Like the old ZGY implementation does.
 *
 *  - Most of the existing bells and whistles in zgycopyc will be
 *    available when creating a mirrored cube. Such as multi-threading,
 *    and reading and/or writing directly from sdms. And using shortcuts
 *    when encountering dead traces. This is obviously a good thing.
 *
 *  - Writing directly to sdms will still slower than it ought to be.
 *    Because low resolution compute reads back lod0 from the output.
 *    This isn't related to the choice of using a virtual wrapper; it is
 *    a consequence of using OpenZGY. Which is needed for upload in any
 *    case. A combination of smarter iteration pattern in zgycopy and
 *    the ongoing OpenZGY performance work will hopefully fix this.
 *
 *    ZgyReaderMirror will automatically crop the input cube to the
 *    brick size. This ought to have been the copy chunk size;
 *    the current implementation will only work with that alignment.
 *    So, you may need to use 64,64,0 as the chunk size if things don't
 *    line up properly.
 *
 *    ZgyReaderMirror can be used by itself. For copying entire files,
 *    it is better to combine it with ZgyWriterMirror. This prevents
 *    the same bricks to be read multiple times.
 *
 * If the upsampling feature should also be ported, this can be done
 * in a similar virtual wrapper. Upsampling (typically vertically) can
 * be chained with mirroring (typically horizontally) if desired.
 * Upsampling should be first in the chain. For this use it becomes
 * even more important that the copy loop only reads the real input
 * file once. I.e. remember to make use of ZgyWriterMirror.
 */
class ZgyReaderMirror : public ZgyReaderRelay
{
public:
  ZgyReaderMirror(
       const std::shared_ptr<IZgyReader>& relay,
       const std::array<int,3>& mirrors);
  virtual ~ZgyReaderMirror();

  // From IZgyReader:
  void read(
       const size3i_t& start, const size3i_t& size,
       float* data, int lod = 0) const override;
  void read(
       const size3i_t& start, const size3i_t& size,
       std::int16_t* data, int lod = 0) const override;
  void read(
       const size3i_t& start, const size3i_t& size,
       std::int8_t* data, int lod = 0) const override;
  std::pair<bool,double> readconst(
       const size3i_t& start, const size3i_t& size,
       int lod = 0, bool as_float = true) const override;

  // From IZgyMeta:
  size3i_t         size()            const override { return inflated_size_; }
  const corners_t  corners()         const override { return worldcorners_; }
  const corners_t  indexcorners()    const override { return indexcorners_; }
  const corners_t  annotcorners()    const override { return annotcorners_; }
  std::vector<size3i_t> brickcount() const override { return brickcount_; }
  int32_t          nlods()           const override { return nlods_; }
  std::string      verid()           const override { return verid_; }
  SampleStatistics statistics()      const override { return statistics_; }
  SampleHistogram  histogram()       const override { return histogram_; }
  std::shared_ptr<const FileStatistics> filestats() const override{return filestats_; }

private:
  std::array<int,3> getReplicant(const size3i_t& start, const size3i_t& size) const;
  std::array<std::int64_t,3> getStart(const size3i_t& start, const size3i_t& size, const std::array<int,3>& replicant) const;

private:
  size3i_t               inflated_size_;
  size3i_t               orig_cropped_size_;
  corners_t              worldcorners_;
  corners_t              indexcorners_;
  corners_t              annotcorners_;
  std::vector<size3i_t>  brickcount_;
  int32_t                nlods_;
  std::string            verid_;
  SampleStatistics       statistics_;
  SampleHistogram        histogram_;
  std::shared_ptr<const FileStatistics> filestats_;
};

class ZgyWriterMirror : public ZgyWriterRelay
{
public:
  typedef ZgyWriterRelay base;
  ZgyWriterMirror(
       const std::shared_ptr<IZgyWriter>& relay,
       const std::array<int,3>& mirrors);
  virtual ~ZgyWriterMirror();
  static std::shared_ptr<IZgyWriter>open(
       const std::string& filename,
       const IOContext* iocontext,
       const std::array<int,3>& mirrors);

private:
  std::array<int,3>      mirrors_;
  size3i_t               orig_cropped_size_;
  static void missing(const std::string& name);

public:
  size3i_t size() const override { return orig_cropped_size_; }
  // Implement for completeness, but nobody needs these yet.
  const corners_t  corners()         const override { missing("corners");      return base::corners(); }
  const corners_t  indexcorners()    const override { missing("indexcorners"); return base::indexcorners(); }
  const corners_t  annotcorners()    const override { missing("annotcorners"); return base::annotcorners(); }
  std::vector<size3i_t> brickcount() const override { missing("brickcount");   return base::brickcount(); }
  int32_t          nlods()           const override { missing("nlods");        return base::nlods(); }

  void write(const size3i_t& start, const size3i_t& size, const float* data)              override;
  void write(const size3i_t& start, const size3i_t& size, const std::int16_t *data)       override;
  void write(const size3i_t& start, const size3i_t& size, const std::int8_t* data)        override;
  void writeconst(const size3i_t& start, const size3i_t& size, const float* data)         override;
  void writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data) override;
  void writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data)   override;
};

}} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
