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

using System.Runtime.InteropServices;

namespace OpenZGY.Managed
{
  public class SampleStatistics : ZgyBase
  {
    internal SampleStatistics(ZgyHandle handle) : base(handle)
    {
    }

    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_statistics_cnt(ZgyHandle handle, ref long value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_statistics_sum(ZgyHandle handle, ref double value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_statistics_ssq(ZgyHandle handle, ref double value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_statistics_min(ZgyHandle handle, ref double value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_statistics_max(ZgyHandle handle, ref double value);

    public long cnt
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_statistics_cnt(Handle, ref value));
        return value;
      }
    }

    public double sum
    {
      get
      {
        double value = 0;
        Tools.freeHandle(oz_statistics_sum(Handle, ref value));
        return value;
      }
    }

    public double ssq
    {
      get
      {
        double value = 0;
        Tools.freeHandle(oz_statistics_ssq(Handle, ref value));
        return value;
      }
    }

    public double min
    {
      get
      {
        double value = 0;
        Tools.freeHandle(oz_statistics_min(Handle, ref value));
        return value;
      }
    }

    public double max
    {
      get
      {
        double value = 0;
        Tools.freeHandle(oz_statistics_max(Handle, ref value));
        return value;
      }
    }
  }
}
