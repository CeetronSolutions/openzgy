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
  public class FileStatistics : ZgyBase
  {
    internal FileStatistics(ZgyHandle handle) : base(handle)
    {
    }

    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_fileVersion(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_fileSize(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_headerSize(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_dataStart(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_alphaNormalCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_alphaNormalSizePerEntry(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_alphaCompressedCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_alphaCcompressedSize(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_alphaMissingCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_alphaConstantCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_brickNormalCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_brickNormalSizePerEntry(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_brickCompressedCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_brickCompressedSize(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_brickMissingCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_brickConstantCount(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_usedSize(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_usedIfUncompressed(ZgyHandle handle, ref long result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_compressionFactor(ZgyHandle handle, ref double result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_isCompressed(ZgyHandle handle, ref int result);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_filestats_toString(ZgyHandle handle);

    public long fileVersion
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_fileVersion(Handle, ref value));
        return value;
      }
    }

    public long fileSize
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_fileSize(Handle, ref value));
        return value;
      }
    }

    public long headerSize
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_headerSize(Handle, ref value));
        return value;
      }
    }

    public long dataStart
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_dataStart(Handle, ref value));
        return value;
      }
    }

    public long alphaNormalCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_alphaNormalCount(Handle, ref value));
        return value;
      }
    }

    public long alphaNormalSizePerEntry
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_alphaNormalSizePerEntry(Handle, ref value));
        return value;
      }
    }

    public long alphaCompressedCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_alphaCompressedCount(Handle, ref value));
        return value;
      }
    }

    public long alphaCcompressedSize
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_alphaCcompressedSize(Handle, ref value));
        return value;
      }
    }

    public long alphaMissingCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_alphaMissingCount(Handle, ref value));
        return value;
      }
    }

    public long alphaConstantCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_alphaConstantCount(Handle, ref value));
        return value;
      }
    }

    public long brickNormalCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_brickNormalCount(Handle, ref value));
        return value;
      }
    }

    public long brickNormalSizePerEntry
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_brickNormalSizePerEntry(Handle, ref value));
        return value;
      }
    }

    public long brickCompressedCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_brickCompressedCount(Handle, ref value));
        return value;
      }
    }

    public long brickCompressedSize
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_brickCompressedSize(Handle, ref value));
        return value;
      }
    }

    public long brickMissingCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_brickMissingCount(Handle, ref value));
        return value;
      }
    }

    public long brickConstantCount
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_brickConstantCount(Handle, ref value));
        return value;
      }
    }

    public long usedSize
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_usedSize(Handle, ref value));
        return value;
      }
    }

    public long usedIfUncompressed
    {
      get
      {
        long value = 0;
        Tools.freeHandle(oz_filestats_usedIfUncompressed(Handle, ref value));
        return value;
      }
    }

    public double compressionFactor
    {
      get
      {
        double value = 0;
        Tools.freeHandle(oz_filestats_compressionFactor(Handle, ref value));
        return value;
      }
    }

    public bool isCompressed
    {
      get
      {
        int value_intflag = 0;
        Tools.freeHandle(oz_filestats_isCompressed(Handle, ref value_intflag));
        return value_intflag != 0;
      }
    }

    public string toString()
    {
      return Tools.resultGetString(oz_filestats_toString(Handle));
    }
  }
}
