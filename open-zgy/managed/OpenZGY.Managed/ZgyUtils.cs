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
  public class ZgyUtils : ZgyBase
  {
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_utils_utils(string prefix, ZgyHandle iocontext);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_utils_close(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_utils_deletefile(ZgyHandle handle, string filename, int missing_ok_intflag);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_utils_alturl(ZgyHandle handle, string filename);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_utils_idtoken(ZgyHandle handle);

    private ZgyUtils(ZgyHandle handle) : base(handle)
    {
    }

    public override string ToString()
    {
      return "ZgyUtils instance";
    }

    public static ZgyUtils utils(string prefix, IOContext iocontext)
    {
      // Do this for methods that don't need an existing handle.
      Tools.checkLibraryVersion();
      var ret = new ZgyUtils(Tools.checkHandle(oz_utils_utils(prefix, ZgyBase.GetHandle(iocontext))));
      ret.Logger = iocontext == null ? null : iocontext.Logger; // For object lifetime only.
      ret.TokenCB = iocontext == null ? null : iocontext.TokenCB; // Ditto.
      return ret;
    }

    public void deleteFile(string filename, bool missing_ok)
    {
      Tools.freeHandle(oz_utils_deletefile(Handle, filename, missing_ok ? 1 : 0));
    }

    public string alturl(string filename)
    {
      return Tools.resultGetString(Tools.checkHandle(oz_utils_alturl(Handle, filename)));
    }

    public string idtoken(string filename)
    {
      return Tools.resultGetString(Tools.checkHandle(oz_utils_idtoken(Handle)));
    }

    public override void Dispose()
    {
      if (Handle != null && !Handle.IsInvalid)
      {
        try
        {
          // What gets freed here is the result of close(),
          // not the writer instance itself.
          Tools.freeHandle(oz_utils_close(Handle));
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
