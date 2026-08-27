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

using System.Runtime.ConstrainedExecution;
using System.Runtime.InteropServices;
using System.Security.Permissions;

// See https://docs.microsoft.com/en-us/dotnet/standard/native-interop/best-practices
// See https://docs.microsoft.com/en-us/dotnet/api/system.runtime.interopservices.safehandle?view=net-6.0

namespace OpenZGY.Managed
{
  [SecurityPermission(SecurityAction.InheritanceDemand, UnmanagedCode = true)]
  public sealed class ZgyHandle : System.Runtime.InteropServices.SafeHandle
  {
    // On Windows: "OpenZGY.dll", case insensitive.
    // On Linux: "libopenzgy.so", case sensitive.
    // The following should work in both cases.
    // If not, consider SetDllImportResolver()
    internal const string OpenZGYLibrary = "openzgy";

    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern void oz_freeHandle(System.IntPtr handle);

    //private static int debug_lastid_ = 0;
    //private readonly int debug_id_;

    [ReliabilityContract(Consistency.WillNotCorruptState, Cer.MayFail)]
    internal ZgyHandle() : base(System.IntPtr.Zero, true)
    {
      //debug_id_ = ++debug_lastid_;
      //System.Console.WriteLine("Create   ZgyHandle {0}", debug_id_);
    }

    public override bool IsInvalid
    {
      [System.Security.SecurityCritical]
      get { return this.handle == System.IntPtr.Zero; }
    }

    protected override bool ReleaseHandle()
    {
      if (handle != System.IntPtr.Zero)
      {
        //System.Console.WriteLine("Release  ZgyHandle {0} {1:x}", debug_id_, handle.ToInt64());
        oz_freeHandle(this.handle);
        return true;
      }
      else
      {
        //System.Console.WriteLine("Release  ZgyHandle {0} (null)", debug_id_);
        return true; // TODO: Not sure... but probably cannot happen.
      }
    }

    protected override void Dispose(bool disposing)
    {
      //System.Console.WriteLine("{0} ZgyHandle {1} {2:x}", disposing ? "Dispose " : "Finalize", debug_id_, handle.ToInt64());
      base.Dispose(disposing);
    }
  }
}

/*
The library was changed to use ZgyHandle:SafeHandle instead of IntPtr.
The switch was not straight forward so it does have some risk.
In case it went wrong, here is the checklist that was used.

In general, any handle expecting to be freed by oz_freeHandle will change.

What needs to change:

* Remove all aliases ZgyHandle = System.IntPtr

* Implement ZgyHandle deriving from SafeHandle

* ZgyHandle calls oz_freeHandle and is the only place where this is done.

* Existing calls to oz_freeHandle(h) change to h.Dispose()
  This include freeHandle (and freeHandleNoThrow, but that goes away)
  Try to ensure the handle goes out of scope. If associated with a ZgyBase
  its stored handle should be set to null.

* A few checks for handle != IntPtr.Zero or ZgyHandle.Zero now test != null.

* All classes deriving from ZgyBase are affected:
  - Main: ZgyMeta, ZgyReader, ZgyWriter, ZgyUtils
  - Returned from ZgyMeta: FileStatistics, SampleHistogram, SampleStatistics
  - Passed from application: ZgyWriterArgs, IOContext

* Other affected ZgyHandles NOT deriving from ZgyBase:
  - Ephemeral: SUCCESS, ERROR, CLEANUP, STRING.

EXCEPTIONS that keep using System.IntPtr

* oz_malloc still returns IntPtr because it is not freed by managed code.

* oz_resultGet{String,ExceptionType,ExceptionMessage} still returns IntPtr
  because they get freed indirectly by the ZgyHandle that owns them.

* InteropTokenCallbackV3 delegate still returns an IntPtr, because this
  is returning a string allocated by oz_malloc() and ownership later
  transfered to the C++ layer. No oz_freeHandle() expected here.

Benefits

* Ensures that all unmanaged resources get freed. This is especially
  important for the lightweight classes such as SampleStatistics
  because it is unreasonable to ask consumers to dispose those after use.

Risks

- SafeHandle is more difficult to use compared to a raw IntPtr. Not as
  difficult as trying to write finalizers by hand that guarantee
  cleanup. But still tricky enough that the change introduces some
  risk.

- The reader and writer classes should really be explicitly disposed,
  preferably with the using() construct. If the finalizer is set up to
  close these automatically then this might hide bugs in the
  application. Because the behavior depends on when or if the
  finalizers decide to run. So, disposing or finalizing just the
  SafeHandle, i.e. not doing this from the ZgyBase Dispose(), should
  do the very minimal amount of work. Tricky because destructing a
  reader or writer on the C++ end will in fact close it.

- Closing a reader or writer from a finalizer might cause SDAPI code
  to be called after connections have been torn down.

Mitigation

- close(), whether called explicitly or from ZgyBase.Dispose(),
  causes the smart pointer on the C++ side to become empty.
  (or maybe I need a bool). (or some check for isOpen).

- ZgyHandle.Dispose called without close/ZgyBase.Dispose() means
  the C++ side can see this is not a controlled shutdown.
  
- ZgyHandleT gets a new member closed_ set when reader/writer/utils close.
  Or maybe just reset the pointer, making free a no-op

- ZgyUtils is missing a close(); is the destructor sufficient? Probably
  BUT I need a close in capi so I can flag a controlled shutdown.

- ZgyHandleBase::freeHandle must ignore an empty pointer.

- ZgyHandleBase::freeHandle print a warning if READER or WRITER and
  the pointer is not empty. And then don't do anything that causes
  network access. This means read and write locks will not be released.
  
- Currently the only way to abandon a reader or writer, preventing
  closing SDAPI handle and releasing locks etc. is to leak the handle.
  And maybe this is good enough. An alternative would be a new
  close_abandon() method.
*/