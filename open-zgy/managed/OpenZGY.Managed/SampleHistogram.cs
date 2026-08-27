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
  public class SampleHistogram : ZgyBase
  {
    internal SampleHistogram(ZgyHandle handle) : base(handle)
    {
    }

    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_histogram_samplecount(ZgyHandle handle, ref long value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_histogram_minvalue(ZgyHandle handle, ref double value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_histogram_maxvalue(ZgyHandle handle, ref double value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_histogram_bins(ZgyHandle handle, [Out] long[] data, ref int actual_size_return, int allocated_size);

    public long samplecount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_histogram_samplecount(Handle, ref value));
        return value;
      }
    }

    public double minvalue
    {
      get
      {
        double value = 0;
        Tools.freeHandle(oz_histogram_minvalue(Handle, ref value));
        return value;
      }
    }

    public double maxvalue
    {
      get
      {
        double value = 0;
        Tools.freeHandle(oz_histogram_maxvalue(Handle, ref value));
        return value;
      }
    }

    public long[] bins()
    {
      int actual_size = 0;
      int allocated_size = 0;
      long[] bins = null;
      oz_histogram_bins(Handle, bins, ref actual_size, allocated_size);
      bins = new long[actual_size];
      allocated_size = actual_size;
      oz_histogram_bins(Handle, bins, ref actual_size, allocated_size);
      return bins;
    }
  }
}
