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
  public class ZgyWriter : ZgyMeta
  {
    private string filename_ = default; // For ToString() only.
    public override string ToString()
    {
      if (Handle == null || Handle.IsInvalid)
        return "ZgyWriter closed instance";
      else if (string.IsNullOrWhiteSpace(filename_))
        return "ZgyWriter instance";
      else
        return string.Format("ZgyWriter(\"{0}\")", filename_);
    }

    private ZgyWriter(ZgyHandle handle) : base(handle)
    {
    }

    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_open(ZgyHandle args);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_reopen(ZgyHandle args);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_finalize(ZgyHandle handle, [In] int[] decimation, int num_decimation, InteropProgressCallback progress, int action, int force_intflag);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_close(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_close_incomplete(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_read_float(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [Out] float[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_read_int16(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [Out] short[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_read_int8(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [Out] byte[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_readconst(ZgyHandle handle, ref int is_const_intflag, ref double value, long i0, long j0, long k0, long ni, long nj, long nk, int as_float_intflag);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_write_float(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [In] float[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_write_int16(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [In] short[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_write_int8(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [In] byte[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_writeconst_float(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [In] float[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_writeconst_int16(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [In] short[] data);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_writer_writeconst_int8(ZgyHandle handle, long i0, long j0, long k0, long ni, long nj, long nk, [In] byte[] data);

    public static ZgyWriter open(ZgyWriterArgs args)
    {
      // Do this for methods that don't need an existing handle.
      Tools.checkLibraryVersion();
      var ret = new ZgyWriter(Tools.checkHandle(oz_writer_open(args.Handle)));
      ret.filename_ = args.FileName; // For ToString() only.
      ret.Logger = args.Logger; // For object lifetime only.
      ret.TokenCB = args.TokenCB; // Ditto.
      return ret;
    }

    public static ZgyWriter reopen(ZgyWriterArgs args)
    {
      // Do this for methods that don't need an existing handle.
      Tools.checkLibraryVersion();
      var ret = new ZgyWriter(Tools.checkHandle(oz_writer_reopen(args.Handle)));
      ret.filename_ = args.FileName; // For ToString() only.
      ret.Logger = args.Logger; // For object lifetime only.
      ret.TokenCB = args.TokenCB; // Ditto.
      return ret;
    }

    // read() and readconst() are identical to those in ZgyWriter except the LOD arguent is missing.
    public void read(long[] start, long[] size, [Out] float[] data)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_writer_read_float(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void read(long[] start, long[] size, [Out] short[] data)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_writer_read_int16(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void read(long[] start, long[] size, [Out] byte[] data)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_writer_read_int8(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public bool readconst(ref double value, long[] start, long[] size, bool as_float)
    {
      int is_const_intflag = 0;
      Tools.freeHandle(oz_writer_readconst(Handle, ref is_const_intflag, ref value, start[0], start[1], start[2], size[0], size[1], size[2], as_float ? 1 : 0));
      return is_const_intflag != 0;
    }

    public void write(long[] start, long[] size, [In] float[] data)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_writer_write_float(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void write(long[] start, long[] size, [In] short[] data)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_writer_write_int16(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void write(long[] start, long[] size, [In] byte[] data)
    {
      validateIOArgs(start, size, data);
      Tools.freeHandle(oz_writer_write_int8(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void writeconst(long[] start, long[] size, float[/*1*/] data)
    {
      validateIOArgs(start, size, data, true);
      Tools.freeHandle(oz_writer_writeconst_float(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void writeconst(long[] start, long[] size, short[/*1*/] data)
    {
      validateIOArgs(start, size, data, true);
      Tools.freeHandle(oz_writer_writeconst_int16(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void writeconst(long[] start, long[] size, byte[/*1*/] data)
    {
      validateIOArgs(start, size, data, true);
      Tools.freeHandle(oz_writer_writeconst_int8(Handle, start[0], start[1], start[2], size[0], size[1], size[2], data));
    }

    public void finalize(DecimationType[] decimation, ProgressDelegate progress, FinalizeAction action, bool force)
    {
      int[] int_decimation = null;
      int int_decimation_size = 0;
      if (decimation != null && decimation.Length > 0)
      {
        int_decimation = new int[decimation.Length];
        int_decimation_size = decimation.Length;
        for (int ii = 0; ii < decimation.Length; ++ii)
          int_decimation[ii] = (int)decimation[ii];
      }
      InteropProgressCallback progressf = Callbacks.forwardProgress(progress);
      Tools.freeHandle(oz_writer_finalize(Handle, int_decimation, int_decimation_size, progressf, (int)action, force ? 1 : 0));
    }

    public void close()
    {
      if (Handle != null && !Handle.IsInvalid)
      {
        try
        {
          // What gets freed here is the result of close(),
          // not the writer instance itself.
          Tools.freeHandle(oz_writer_close(Handle));
        }
        finally
        {
          // Dispose the handle now, freeing the writer instance.
          // No need to wait for the user to dispose. Because
          // the instance cannot be used for anything anyway.
          // Only dispose the base class to avoid infinite recursion.
          base.Dispose();
        }
      }
    }

    public void close_incomplete()
    {
      // See close() for explanations.
      if (Handle != null && !Handle.IsInvalid)
      {
        try
        {
          Tools.freeHandle(oz_writer_close_incomplete(Handle));
        }
        finally
        {
          base.Dispose();
        }
      }
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
      close(); // will also call base.Dispose()
    }
  }
}
