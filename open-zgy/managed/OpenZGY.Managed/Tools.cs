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
  internal static class Tools
  {
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern ZgyHandle oz_checkLibraryVersion(uint version);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern int oz_resultIsSuccess(ZgyHandle handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern System.IntPtr oz_resultGetExceptionType(ZgyHandle Handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern System.IntPtr oz_resultGetExceptionMessage(ZgyHandle Handle);
    [DllImport(ZgyHandle.OpenZGYLibrary)] private static extern System.IntPtr oz_resultGetString(ZgyHandle Handle);

    internal static void checkLibraryVersion()
    {
      freeHandle(oz_checkLibraryVersion(0x000101));
    }

    /// <summary>
    /// Convert error messages back to exceptions.
    /// The handle will be freed if an exception is thrown.
    /// Otherwise the handle remains valid and is returned
    /// as a convenience.
    /// </summary>
    internal static ZgyHandle checkHandle(ZgyHandle handle)
    {
      if (handle == null || handle.IsInvalid)
      {
        // An external function that should return a long lived handle
        // should never return null. For short lived handles null means
        // success, so test for that before calling checkHandle().
        throw new Errors.ZgyInternalError("OpenZGY.Managed: unexpected invalid handle");
      }
      else if (oz_resultIsSuccess(handle) == 0)
      {
        System.IntPtr typename_ptr = oz_resultGetExceptionType(handle);
        System.IntPtr message_ptr = oz_resultGetExceptionMessage(handle);
        string typename = Marshal.PtrToStringAnsi(typename_ptr);
        string message = Marshal.PtrToStringAnsi(message_ptr);
        handle.Dispose(); // Also frees the unmanaged string pointers.
        // The assumption is that handle now goes out of scope.
        throw Errors.ZgyError.MakeException(typename, message);
      }
      return handle;
    }

    /**
     * Convert a ZgyHandle, returned from a function that in C++
     * returns a string, back to an actual string in C#.
     * Does not throw, but crashes on garbage pointer.
     * oz_resultGetString() returns "" if handle is wrong,
     * typically an error handle.
     * Free the handle unconditionally. The assumption is
     * that the caller will not retain a reference to it.
     */
    internal static string resultGetString(ZgyHandle handle)
    {
      if (handle == null || handle.IsInvalid)
        throw new Errors.ZgyInternalError("resultGetString: Invalid handle");
      checkHandle(handle); // Disposes and throws on error.
      System.IntPtr result_ptr = oz_resultGetString(handle);
      string result = Marshal.PtrToStringAnsi(result_ptr);
      handle.Dispose(); // Free result_ptr and the handle itself.
      return result;
    }

    /// <summary>
    /// Convert error messages back to exceptions.
    /// Free the handle unconditionally.
    /// </summary>
    /// <details>
    /// If the handle has already been checked once for errors then it is
    /// safe to assume it will not throw. Note that if checkHandle() throws
    /// then the call to Dispose() is skipped. This is correct because on
    /// error it is the responsibility of checkHandle() to free it.
    /// </details>
    internal static void freeHandle(ZgyHandle handle)
    {
      if (handle != null)
        checkHandle(handle).Dispose();
    }
  }
}
