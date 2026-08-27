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
  // TODO, finalizer? Or use SafeHandle which is more robust.

  [System.Security.SuppressUnmanagedCodeSecurity()]
  public class ZgyWriterArgs : ZgyBase
  {
    // Only for ToString(), we don't otherwise need it.
    internal string FileName { get; set; } = null;

    public ZgyWriterArgs() : base(Tools.checkHandle(oz_writerargs_create()))
    {
    }

    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_create();
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_filename(ZgyHandle handle, string name);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_iocontext(ZgyHandle handle, ZgyHandle iocontext);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_zfp_compressor(ZgyHandle handle, float sqnr);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_zfp_lodcompressor(ZgyHandle handle, float sqnr);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_size(ZgyHandle handle, long ni, long nj, long nk);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_bricksize(ZgyHandle handle, long ni, long nj, long nk);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_datatype(ZgyHandle handle, int enum_as_int);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_datarange(ZgyHandle handle, float lo, float hi);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_zunit(ZgyHandle handle, int enum_as_int, string name, double factor);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_hunit(ZgyHandle handle, int enum_as_int, string name, double factor);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_ilstart(ZgyHandle handle, float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_ilinc(ZgyHandle handle, float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_xlstart(ZgyHandle handle, float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_xlinc(ZgyHandle handle, float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_zstart(ZgyHandle handle, float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_zinc(ZgyHandle handle, float value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_corners(ZgyHandle handle, double X00, double Y00, double XN0, double YN0, double X0M, double Y0M, double XNM, double YNM);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_metafrom(ZgyHandle handle, ZgyHandle other);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_merge(ZgyHandle handle, ZgyHandle other);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writerargs_toString(ZgyHandle handle);

    public ZgyWriterArgs filename(string name)
    {
      Tools.freeHandle(oz_writerargs_filename(Handle, name));
      FileName = name;
      return this;
    }

    public ZgyWriterArgs iocontext(IOContext iocontext)
    {
      Tools.freeHandle(oz_writerargs_iocontext(Handle, ZgyBase.GetHandle(iocontext)));
      Logger = iocontext == null ? null : iocontext.Logger;
      TokenCB = iocontext == null ? null : iocontext.TokenCB;
      return this;
    }

    public ZgyWriterArgs zfp_compressor(float sqnr)
    {
      Tools.freeHandle(oz_writerargs_zfp_compressor(Handle, sqnr));
      return this;
    }

    public ZgyWriterArgs zfp_lodcompressor(float sqnr)
    {
      Tools.freeHandle(oz_writerargs_zfp_lodcompressor(Handle, sqnr));
      return this;
    }

    public ZgyWriterArgs size(long ni, long nj, long nk)
    {
      Tools.freeHandle(oz_writerargs_size(Handle, ni, nj, nk));
      return this;
    }

    public ZgyWriterArgs bricksize(long ni, long nj, long nk)
    {
      Tools.freeHandle(oz_writerargs_bricksize(Handle, ni, nj, nk));
      return this;
    }

    public ZgyWriterArgs datatype(SampleDataType dt)
    {
      Tools.freeHandle(oz_writerargs_datatype(Handle, (int)dt));
      return this;
    }

    public ZgyWriterArgs datarange(float lo, float hi)
    {
      Tools.freeHandle(oz_writerargs_datarange(Handle, lo, hi));
      return this;
    }

    public ZgyWriterArgs zunit(UnitDimension dim, string name, double factor)
    {
      Tools.freeHandle(oz_writerargs_zunit(Handle, (int)dim, name, factor));
      return this;
    }

    public ZgyWriterArgs hunit(UnitDimension dim, string name, double factor)
    {
      Tools.freeHandle(oz_writerargs_hunit(Handle, (int)dim, name, factor)); // (ZgyHandle handle, ENUM dimension, const char* name, double factor)
      return this;
    }

    public ZgyWriterArgs ilstart(float value)
    {
      Tools.freeHandle(oz_writerargs_ilstart(Handle, value));
      return this;
    }

    public ZgyWriterArgs ilinc(float value)
    {
      Tools.freeHandle(oz_writerargs_ilinc(Handle, value));
      return this;
    }

    public ZgyWriterArgs xlstart(float value)
    {
      Tools.freeHandle(oz_writerargs_xlstart(Handle, value));
      return this;
    }

    public ZgyWriterArgs xlinc(float value)
    {
      Tools.freeHandle(oz_writerargs_xlinc(Handle, value));
      return this;
    }

    public ZgyWriterArgs zstart(float value)
    {
      Tools.freeHandle(oz_writerargs_zstart(Handle, value));
      return this;
    }

    public ZgyWriterArgs zinc(float value)
    {
      Tools.freeHandle(oz_writerargs_zinc(Handle, value));
      return this;
    }

    public ZgyWriterArgs corners(double[/*4*/][/*2*/] corners)
    {
      Tools.freeHandle(oz_writerargs_corners(
        Handle,
        corners[0][0], corners[0][1],
          corners[1][0], corners[1][1],
          corners[2][0], corners[2][1],
          corners[3][0], corners[3][1]));
      return this;
    }

    public ZgyWriterArgs metafrom(ZgyReader other)
    {
      Tools.freeHandle(oz_writerargs_metafrom(Handle, other.Handle));
      return this;
    }

    public ZgyWriterArgs merge(ZgyWriterArgs other)
    {
      Tools.freeHandle(oz_writerargs_merge(Handle, other.Handle));
      return this;
    }

    public string toString()
    {
      return Tools.resultGetString(oz_writerargs_toString(Handle));
    }
  }
}
