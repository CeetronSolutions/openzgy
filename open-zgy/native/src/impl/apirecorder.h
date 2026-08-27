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

#pragma once

#include "../api.h"
#include "../tools/readwriterelay.h"

#include <mutex>
#include <string>
#include <sstream>
#include <ostream>
#include <memory>
#include <utility>
#include <atomic>
#include <cstdint>
#include <array>

#ifdef _MSC_VER
#pragma warning(push)
// The warning is due to me using pure virtual classes as "real" interfaces,
// having a similar inheritance scheme for interfaces and implementation.
// I believe this is safe as long as the interfaces are truly pure.
#pragma warning(disable:4250) // inherits via dominance
#endif

namespace OpenZGY { namespace Telemetry {
#if 0
}}
#endif

class ApiRecorder
{
public:
  typedef std::array<std::int64_t,3> size3i_t;
  explicit ApiRecorder(const std::string& name);
  ~ApiRecorder();
  void logger(const std::string& str) const;
  void logger(const std::ios& ss) const;
  void logio(
       const std::string prefix,
       const size3i_t& start,
       const size3i_t& size,
       int lod) const;
  static bool enabled();

private:
  static std::atomic<int> next_id_;
  static std::mutex global_mutex_;
  static std::weak_ptr<std::ostream> global_os_;

private:
  int id_;
  std::string name_;
  std::shared_ptr<std::ostream> os_;

private:
  static std::shared_ptr<std::ostream> getVerboseFileFromEnv(int id);
  static double timestamp();
};

class ApiReaderRecorder : public OpenZGY::Tools::ZgyReaderRelay
{
public:
    typedef std::array<std::int64_t,3> size3i_t;

private:
  std::shared_ptr<ApiRecorder> recorder_;
  void logger(const std::string& str) const { recorder_->logger(str); }
  void logger(const std::ios& ss) const { recorder_->logger(ss); }
  void logio(
       const std::string prefix,
       const size3i_t& start,
       const size3i_t& size,
       int lod) const
  {
    recorder_->logio(prefix, start, size, lod);
  }

  ApiReaderRecorder(const std::shared_ptr<IZgyReader>& relay, const std::string& name);

public:
  virtual ~ApiReaderRecorder();

  void read(const size3i_t& start, const size3i_t& size, float*        data, int lod) const override;
  void read(const size3i_t& start, const size3i_t& size, std::int16_t* data, int lod) const override;
  void read(const size3i_t& start, const size3i_t& size, std::int8_t*  data, int lod) const override;
  std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, int lod, bool as_float) const override;
  void close() override;

public:
  static std::shared_ptr<IZgyReader> inject(
       std::shared_ptr<IZgyReader> in,
       const std::string& name);
};

class ApiWriterRecorder : public OpenZGY::Tools::ZgyWriterRelay
{
public:
  //typedef std::array<std::int64_t,3> size3i_t;
  typedef OpenZGY::Tools::ZgyWriterRelay::size3i_t size3i_t;

private:
  std::shared_ptr<ApiRecorder> recorder_;
  void logger(const std::string& str) const { recorder_->logger(str); }
  void logger(const std::ios& ss) const { recorder_->logger(ss); }
  void logio(
       const std::string prefix,
       const size3i_t& start,
       const size3i_t& size,
       int lod) const
  {
    recorder_->logio(prefix, start, size, lod);
  }

  ApiWriterRecorder(const std::shared_ptr<IZgyWriter>& relay, const std::string& name);

public:
  virtual ~ApiWriterRecorder();

  void read(const size3i_t& start, const size3i_t& size, float*        data) const override;
  void read(const size3i_t& start, const size3i_t& size, std::int16_t* data) const override;
  void read(const size3i_t& start, const size3i_t& size, std::int8_t*  data) const override;
  std::pair<bool,double> readconst(const size3i_t& start, const size3i_t& size, bool as_float) const override;

  void write(const size3i_t& start, const size3i_t& size, const float* data)              override;
  void write(const size3i_t& start, const size3i_t& size, const std::int16_t *data)       override;
  void write(const size3i_t& start, const size3i_t& size, const std::int8_t* data)        override;
  void writeconst(const size3i_t& start, const size3i_t& size, const float* data)         override;
  void writeconst(const size3i_t& start, const size3i_t& size, const std::int16_t * data) override;
  void writeconst(const size3i_t& start, const size3i_t& size, const std::int8_t* data)   override;
  void close_incomplete() override;
  void close() override;
  void finalize(
       const std::vector<DecimationType>& decimation,
       const std::function<bool(std::int64_t,std::int64_t)>& progress,
       FinalizeAction action,
       bool force) override;

public:
  static std::shared_ptr<IZgyWriter> inject(
       std::shared_ptr<IZgyWriter> in,
       const std::string& name);
};

}} // namespace

#ifdef _MSC_VER
#pragma warning(pop)
#endif
