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
  public class ZgyMeta : ZgyBase
  {
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_size(ZgyHandle handle, ref long ni_return, ref long nj_return, ref long nk_return);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_datatype(ZgyHandle handle, ref int enum_as_int);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_datarange(ZgyHandle handle, ref float lo_return, ref float hi_return);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_raw_datarange(ZgyHandle handle, ref float lo_return, ref float hi_return);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_zunitdim(ZgyHandle handle, ref int enum_as_int);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_hunitdim(ZgyHandle handle, ref int enum_as_int);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_zunitname(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_hunitname(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_zunitfactor(ZgyHandle handle, ref double value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_hunitfactor(ZgyHandle handle, ref double value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_zstart(ZgyHandle handle, ref float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_zinc(ZgyHandle handle, ref float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_annotstart(ZgyHandle handle, ref float il_return, ref float xl_return);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_annotinc(ZgyHandle handle, ref float il_return, ref float xl_return);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_corners(ZgyHandle handle, ref double X00, ref double Y00, ref double XN0, ref double YN0, ref double X0M, ref double Y0M, ref double XNM, ref double YNM);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_indexcorners(ZgyHandle handle, ref double X00, ref double Y00, ref double XN0, ref double YN0, ref double X0M, ref double Y0M, ref double XNM, ref double YNM);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_annotcorners(ZgyHandle handle, ref double X00, ref double Y00, ref double XN0, ref double YN0, ref double X0M, ref double Y0M, ref double XNM, ref double YNM);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_bricksize(ZgyHandle handle, ref long ni_return, ref long nj_return, ref long nk_return);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_nlods(ZgyHandle handle, ref int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_verid(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_statistics(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_histogram(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_filestats(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_meta_toString(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_tools_indexToAnnot(ZgyHandle handle, ref double x_return, ref double y_return, double x, double y);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_tools_indexToWorld(ZgyHandle handle, ref double x_return, ref double y_return, double x, double y);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_tools_annotToIndex(ZgyHandle handle, ref double x_return, ref double y_return, double x, double y);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_tools_annotToWorld(ZgyHandle handle, ref double x_return, ref double y_return, double x, double y);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_tools_worldToIndex(ZgyHandle handle, ref double x_return, ref double y_return, double x, double y);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_tools_worldToAnnot(ZgyHandle handle, ref double x_return, ref double y_return, double x, double y);

    protected ZgyMeta(ZgyHandle handle) : base(handle)
    {
    }

    protected void validateIOArgs(long[] start, long[] size, object data)
    {
      validateIOArgs(start, size, data, false);
    }

    protected void validateIOArgs(long[] start, long[] size, object data, bool scalar)
    {
      if (start == null)
        throw new System.ArgumentNullException("start", "start should be a long[3] but was null");
      else if (start.Length != 3)
        throw new System.ArgumentException(string.Format("start should be a long[3] but was long[{0}]", start.Length), "start");
      if (size == null)
        throw new System.ArgumentNullException("size", "size should be a long[3] but was null");
      else if (size.Length != 3)
        throw new System.ArgumentException(string.Format("size should be a long[3] but was long[{0}]", start.Length), "size");
      long expect_len = scalar ? 1 : size[0] * size[1] * size[2];
      if (data == null)
        throw new System.NullReferenceException(string.Format("read() expected array[{0}], got null", expect_len));
      long actual_len;
      if (data is float[])
        actual_len = (data as float[]).LongLength;
      else if (data is short[])
        actual_len = (data as short[]).LongLength;
      else if (data is byte[])
        actual_len = (data as byte[]).LongLength;
      else
        throw new System.ArgumentException("data should be float[], short[], or byte[], not " + data.GetType().Name);
      if (expect_len != actual_len)
        throw new System.ArgumentException(string.Format("expected array[{0}], got array[{1}]", expect_len, actual_len));
    }

    public long[/*3*/] size
    {
      get
      {
        {
          var ret = new long[3];
          Tools.freeHandle(oz_meta_size(Handle, ref ret[0], ref ret[1], ref ret[2]));
          return ret;
        }
      }
    }

    public SampleDataType datatype
    {
      get
      {
        int ret = default(int);
        Tools.freeHandle(oz_meta_datatype(Handle, ref ret));
        return (SampleDataType)ret;
      }
    }

    public float[/*2*/] datarange
    {
      get
      {
        var ret = new float[2];
        Tools.freeHandle(oz_meta_datarange(Handle, ref ret[0], ref ret[1]));
        return ret;
      }
    }

    public float[/*2*/] raw_datarange
    {
      get
      {
        var ret = new float[2];
        Tools.freeHandle(oz_meta_raw_datarange(Handle, ref ret[0], ref ret[1]));
        return ret;
      }
    }

    public UnitDimension zunitdim
    {
      get
      {
        int ret = default;
        Tools.freeHandle(oz_meta_zunitdim(Handle, ref ret));
        return (UnitDimension)ret;
      }
    }

    public UnitDimension hunitdim
    {
      get
      {
        int ret = default;
        Tools.freeHandle(oz_meta_hunitdim(Handle, ref ret));
        return (UnitDimension)ret;
      }
    }

    public string zunitname
    {
      get
      {
        return Tools.resultGetString(Tools.checkHandle(oz_meta_zunitname(Handle)));
      }
    }

    public string hunitname
    {
      get
      {
        return Tools.resultGetString(Tools.checkHandle(oz_meta_hunitname(Handle)));
      }
    }

    public double zunitfactor
    {
      get
      {
        double ret = default;
        Tools.freeHandle(oz_meta_zunitfactor(Handle, ref ret));
        return ret;
      }
    }

    public double hunitfactor
    {
      get
      {
        double ret = default;
        Tools.freeHandle(oz_meta_hunitfactor(Handle, ref ret));
        return ret;
      }
    }

    public float zstart
    {
      get
      {
        float ret = default;
        Tools.freeHandle(oz_meta_zstart(Handle, ref ret));
        return ret;
      }
    }

    public float zinc
    {
      get
      {
        float ret = default;
        Tools.freeHandle(oz_meta_zinc(Handle, ref ret));
        return ret;
      }
    }

    public float[/*2*/] annotstart
    {
      get
      {
        var ret = new float[2];
        Tools.freeHandle(oz_meta_annotstart(Handle, ref ret[0], ref ret[1]));
        return ret;
      }
    }

    public float[/*2*/] annotinc
    {
      get
      {
        var ret = new float[2];
        Tools.freeHandle(oz_meta_annotinc(Handle, ref ret[0], ref ret[1]));
        return ret;
      }
    }

    public double[/*4*/][/*2*/] corners
    {
      get
      {
        double[][] ret = new double[4][];
        ret[0] = new double[2];
        ret[1] = new double[2];
        ret[2] = new double[2];
        ret[3] = new double[2];
        Tools.freeHandle(oz_meta_corners(Handle,
          ref ret[0][0], ref ret[0][1],
          ref ret[1][0], ref ret[1][1],
          ref ret[2][0], ref ret[2][1],
          ref ret[3][0], ref ret[3][1]));
        return ret;
      }
    }

    public double[/*4*/][/*2*/] indexcorners
    {
      get
      {
        double[][] ret = new double[4][];
        ret[0] = new double[2];
        ret[1] = new double[2];
        ret[2] = new double[2];
        ret[3] = new double[2];
        Tools.freeHandle(oz_meta_indexcorners(Handle,
          ref ret[0][0], ref ret[0][1],
          ref ret[1][0], ref ret[1][1],
          ref ret[2][0], ref ret[2][1],
          ref ret[3][0], ref ret[3][1]));
        return ret;
      }
    }

    public double[/*4*/][/*2*/] annotcorners
    {
      get
      {
        double[][] ret = new double[4][];
        ret[0] = new double[2];
        ret[1] = new double[2];
        ret[2] = new double[2];
        ret[3] = new double[2];
        Tools.freeHandle(oz_meta_annotcorners(Handle,
          ref ret[0][0], ref ret[0][1],
          ref ret[1][0], ref ret[1][1],
          ref ret[2][0], ref ret[2][1],
          ref ret[3][0], ref ret[3][1]));
        return ret;
      }
    }

    public long[/*3*/] bricksize
    {
      get
      {
        var ret = new long[3];
        Tools.freeHandle(oz_meta_bricksize(Handle, ref ret[0], ref ret[1], ref ret[2]));
        return ret;
      }
    }

    public int nlods
    {
      get
      {
        int ret = default;
        Tools.freeHandle(oz_meta_nlods(Handle, ref ret));
        return ret;
      }
    }

    public string verid
    {
      get
      {
        return Tools.resultGetString(oz_meta_verid(Handle));
      }
    }

    /*
        Replace ^\(.*\)\((.*\)$
        public TYPE \1
        {
          get
          {
            TYPE ret = default;
            Tools.freeHandle(\1(Handle, ref ret)); // \2
            return ret;
          }
        }
    */

    public string toString()
    {
      return Tools.resultGetString(oz_meta_toString(Handle));
    }

    public SampleStatistics statistics()
    {
      return new SampleStatistics(Tools.checkHandle(oz_meta_statistics(Handle)));
    }

    public SampleHistogram histogram()
    {
      return new SampleHistogram(Tools.checkHandle(oz_meta_histogram(Handle)));
    }

    public FileStatistics filestats()
    {
      return new FileStatistics(Tools.checkHandle(oz_meta_filestats(Handle)));
    }

    public double[/*2*/] indexToAnnot(double[/*2*/] input)
    {
      var result = new double[2];
      Tools.freeHandle(oz_tools_indexToAnnot(Handle, ref result[0], ref result[1], input[0], input[1]));
      return result;
    }

    public double[/*2*/] indexToWorld(double[/*2*/] input)
    {
      var result = new double[2];
      Tools.freeHandle(oz_tools_indexToWorld(Handle, ref result[0], ref result[1], input[0], input[1]));
      return result;
    }

    public double[/*2*/] annotToIndex(double[/*2*/] input)
    {
      var result = new double[2];
      Tools.freeHandle(oz_tools_annotToIndex(Handle, ref result[0], ref result[1], input[0], input[1]));
      return result;
    }

    public double[/*2*/] annotToWorld(double[/*2*/] input)
    {
      var result = new double[2];
      Tools.freeHandle(oz_tools_annotToWorld(Handle, ref result[0], ref result[1], input[0], input[1]));
      return result;
    }

    public double[/*2*/] worldToIndex(double[/*2*/] input)
    {
      var result = new double[2];
      Tools.freeHandle(oz_tools_worldToIndex(Handle, ref result[0], ref result[1], input[0], input[1]));
      return result;
    }

    public double[/*2*/] worldToAnnot(double[/*2*/] input)
    {
      var result = new double[2];
      Tools.freeHandle(oz_tools_worldToAnnot(Handle, ref result[0], ref result[1], input[0], input[1]));
      return result;
    }
  }
}
