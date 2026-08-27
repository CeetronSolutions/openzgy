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

  public class IOContext : ZgyBase
  {
    protected IOContext(ZgyHandle handle) : base(handle)
    {
    }
  }

  [System.Security.SuppressUnmanagedCodeSecurity()]
  public class SeismicStoreIOContext : IOContext
  {
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_create();
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_sdurl(ZgyHandle handle, string value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_sdapikey(ZgyHandle handle, string value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_sdtoken(ZgyHandle handle, string value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_sdtokencb(ZgyHandle handle, InteropTokenCallbackV3 callback);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_maxsize(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_maxhole(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_aligned(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_segsize(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_segsplit(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_iothreads(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_cputhreads(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_writethreads(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_legaltag(ZgyHandle handle, string value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_writeid(ZgyHandle handle, string value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_seismicmeta(ZgyHandle handle, string value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_setRoAfterWrite(ZgyHandle handle, int value_intflag);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_forceRoBeforeRead(ZgyHandle handle, int value_intflag);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_forceRwBeforeWrite(ZgyHandle handle, int value_intflag);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_retryCount(ZgyHandle handle, int value);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_logger(ZgyHandle handle, InteropLoggerCallback callback);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_toString(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_credentialsFrom(ZgyHandle handle, ZgyHandle utils_handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_ssiocontext_clone(ZgyHandle handle);

    public SeismicStoreIOContext() : base(Tools.checkHandle(oz_ssiocontext_create()))
    {
    }

    private SeismicStoreIOContext(ZgyHandle handle) : base(handle)
    {
    }

    public SeismicStoreIOContext sdurl(string value)
    {
      Tools.freeHandle(oz_ssiocontext_sdurl(Handle, value));
      return this;
    }

    public SeismicStoreIOContext sdapikey(string value)
    {
      Tools.freeHandle(oz_ssiocontext_sdapikey(Handle, value));
      return this;
    }

    public SeismicStoreIOContext sdtoken(string value)
    {
      Tools.freeHandle(oz_ssiocontext_sdtoken(Handle, value));
      return this;
    }

    public SeismicStoreIOContext maxsize(int value)
    {
      Tools.freeHandle(oz_ssiocontext_maxsize(Handle, value));
      return this;
    }

    public SeismicStoreIOContext maxhole(int value)
    {
      Tools.freeHandle(oz_ssiocontext_maxhole(Handle, value));
      return this;
    }

    public SeismicStoreIOContext aligned(int value)
    {
      Tools.freeHandle(oz_ssiocontext_aligned(Handle, value));
      return this;
    }

    public SeismicStoreIOContext segsize(int value)
    {
      Tools.freeHandle(oz_ssiocontext_segsize(Handle, value));
      return this;
    }

    public SeismicStoreIOContext segsplit(int value)
    {
      Tools.freeHandle(oz_ssiocontext_segsplit(Handle, value));
      return this;
    }

    public SeismicStoreIOContext iothreads(int value)
    {
      Tools.freeHandle(oz_ssiocontext_iothreads(Handle, value));
      return this;
    }

    public SeismicStoreIOContext cputhreads(int value)
    {
      Tools.freeHandle(oz_ssiocontext_cputhreads(Handle, value));
      return this;
    }

    public SeismicStoreIOContext writethreads(int value)
    {
      Tools.freeHandle(oz_ssiocontext_writethreads(Handle, value));
      return this;
    }

    public SeismicStoreIOContext legaltag(string value)
    {
      Tools.freeHandle(oz_ssiocontext_legaltag(Handle, value));
      return this;
    }

    public SeismicStoreIOContext writeid(string value)
    {
      Tools.freeHandle(oz_ssiocontext_writeid(Handle, value));
      return this;
    }

    public SeismicStoreIOContext seismicmeta(string value)
    {
      Tools.freeHandle(oz_ssiocontext_seismicmeta(Handle, value));
      return this;
    }

    public SeismicStoreIOContext setRoAfterWrite(bool value)
    {
      Tools.freeHandle(oz_ssiocontext_setRoAfterWrite(Handle, value ? 1 : 0));
      return this;
    }

    public SeismicStoreIOContext forceRoBeforeRead(bool value)
    {
      Tools.freeHandle(oz_ssiocontext_forceRoBeforeRead(Handle, value ? 1 : 0));
      return this;
    }

    public SeismicStoreIOContext forceRwBeforeWrite(bool value)
    {
      Tools.freeHandle(oz_ssiocontext_forceRwBeforeWrite(Handle, value ? 1 : 0));
      return this;
    }

    public SeismicStoreIOContext retryCount(int value)
    {
      Tools.freeHandle(oz_ssiocontext_retryCount(Handle, value));
      return this;
    }

    public SeismicStoreIOContext clone()
    {
      return new SeismicStoreIOContext(Tools.checkHandle(oz_ssiocontext_clone(Handle)));
    }

    public SeismicStoreIOContext logger(LoggerDelegate fn)
    {
      Logger = fn;
      Tools.freeHandle(oz_ssiocontext_logger(Handle, Callbacks.forwardLogger(fn)));
      return this;
    }

    public SeismicStoreIOContext sdtokencb(TokenDelegate fn)
    {
      TokenCB = fn;
      Tools.freeHandle(oz_ssiocontext_sdtokencb(Handle, Callbacks.forwardTokenV3(fn)));
      return this;
    }

    public string toString()
    {
      return Tools.resultGetString(oz_ssiocontext_toString(Handle));
    }

    public SeismicStoreIOContext credentialsFrom(ZgyUtils utils)
    {
      Tools.freeHandle(oz_ssiocontext_credentialsFrom(Handle, utils.Handle));
      return this;
    }
  }
}
