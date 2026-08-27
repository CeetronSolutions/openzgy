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
  [System.Security.SuppressUnmanagedCodeSecurity()]
  public class ZgyReader : ZgyMeta
  {
    private ZgyReader(ZgyHandle handle): base(handle)
    {
    }

    private string filename_ = default; // For ToString() only.
    public override string ToString()
    {
      if (Handle == null || Handle.IsInvalid)
        return "ZgyReader closed instance";
      else if (string.IsNullOrWhiteSpace(filename_))
        return "ZgyReader instance";
      else
        return string.Format("ZgyReader(\"{0}\")", filename_);
    }

    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_reader_open(string name, ZgyHandle iocontext);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_reader_close(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_reader_read_float(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [Out] float[] data, int lod);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_reader_read_int16(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [Out] short[] data, int lod);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_reader_read_int8(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [Out] byte[] data, int lod);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_reader_readconst(ZgyHandle handle, ref int is_const_intflag, ref double value, long i0, long j0, long k0, long ni, long nj, long nk, int lod, int as_float_intflag);

    public static ZgyReader open(string name, IOContext iocontext)
    {
      // Do this for methods that don't need an existing handle.
      Tools.checkLibraryVersion();
      var ret = new ZgyReader(Tools.checkHandle(oz_reader_open(name, ZgyBase.GetHandle(iocontext)))); 
      ret.filename_ = name; // For ToString() only.
      ret.Logger = iocontext == null ? null : iocontext.Logger; // For object lifetime only.
      ret.TokenCB = iocontext == null ? null : iocontext.TokenCB; // Ditto.
      return ret;
    }

    public void close()
    {
      // The application should always use "using" to ensure proper cleanup.
      // Technically the close() could be skipped (since Dispose calls it)
      // but for consistency with the C++ API it should probably be used.
      Dispose();
    }

    public void read(long[] start, long[] size, [Out] float[] data, int lod)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_reader_read_float(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data, lod));
    }

    public void read(long[] start, long[] size, [Out] short[] data, int lod)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_reader_read_int16(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data, lod));
    }

    public void read(long[] start, long[] size, [Out] byte[] data, int lod)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_reader_read_int8(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data, lod));
    }

    public bool readconst(ref double value, long[] start, long[] size, int lod, bool as_float)
    {
      int is_const_intflag = 0;
      Tools.freeHandle(oz_reader_readconst(Handle, ref is_const_intflag, ref value, start[0], start[1], start[2], size[0], size[1], size[2], lod, as_float ? 1 : 0));
      return is_const_intflag != 0;
    }

    /// <summary>
    /// Release all resources, even if something throws.
    /// Safe to call more than once.
    /// If the application forgot to close before disposing this will
    /// silently be done now. This matches the behavior in C++ where
    /// a reader or writer going out of scope will silently close.
    /// </summary>
    public override void Dispose()
    {
      if (Handle != null && !Handle.IsInvalid)
      {
        try
        {
          // What gets freed here is the result of close(),
          // not the writer instance itself.
          Tools.freeHandle(oz_reader_close(Handle));
        }
        finally
        {
          // Dispose the handle now, freeing the reader instance.
          base.Dispose();
        }
      }
    }
  }
}
