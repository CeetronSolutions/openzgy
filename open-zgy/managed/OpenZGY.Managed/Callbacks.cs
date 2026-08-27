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
  public delegate bool ProgressDelegate(long done, long total);
  public delegate bool LoggerDelegate(int level, string message);
  public delegate string TokenDelegate();

  // Marshal bool as int in the C API.
  // loggercb: bool(*progress)(int, std::string)
  // progresscb: bool(*progress)(int)
  // tokencb: const char*(*tokencb)() or std::size_t(* tokencb)(char*, std::size_t)

  [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  public delegate int InteropLoggerCallback(int level, string message);

  [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  public delegate int InteropProgressCallback(long done, long total);

  [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  public delegate string InteropTokenCallbackV1();

  [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  public delegate long InteropTokenCallbackV2(long buffersize, string message);

  [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  public delegate System.IntPtr InteropTokenCallbackV3();

  public class Callbacks
  {
    [DllImport(ZgyHandle.OpenZGYLibrary)] internal static extern System.IntPtr oz_malloc(long size);
    [DllImport(ZgyHandle.OpenZGYLibrary)] internal static extern ZgyHandle oz_test_example_v1(InteropLoggerCallback logger, InteropProgressCallback progress, InteropTokenCallbackV1 tokencb);
    [DllImport(ZgyHandle.OpenZGYLibrary)] internal static extern ZgyHandle oz_test_example_v2(InteropLoggerCallback logger, InteropProgressCallback progress, InteropTokenCallbackV2 tokencb);
    [DllImport(ZgyHandle.OpenZGYLibrary)] internal static extern ZgyHandle oz_test_example_v3(InteropLoggerCallback logger, InteropProgressCallback progress, InteropTokenCallbackV3 tokencb);

    static public LoggerDelegate EmptyCallback()
    {
      return (int level, string message) => false;
    }

    static public LoggerDelegate StandardCallback(int verbose, string prefix)
    {
      return (int level, string message) =>
      {
        if (level >= verbose && !string.IsNullOrEmpty(message))
        {
          string[] split = message.Split(new char[] { '\r', '\n' }, System.StringSplitOptions.RemoveEmptyEntries);
          var ss = new System.Text.StringBuilder();
          foreach (var line in split)
            ss.AppendLine(prefix + line);
          if (!message.EndsWith(System.Environment.NewLine))
            message += System.Environment.NewLine;
          System.Console.Write(ss.ToString());
          System.Console.Out.Flush();
        }
        return level >= verbose;
      };
    }

    /**
     * The V1 token callback is directly usable while V1 and V2 require
     * manual marshalling and forwarding to the application's callback type.
     */
    public static void TestExampleV1(LoggerDelegate logger, ProgressDelegate progress, TokenDelegate tokencb)
    {
      Tools.freeHandle(oz_test_example_v1(forwardLogger(logger), forwardProgress(progress), forwardTokenV1(tokencb)));
    }

    /* forwardTokenV2 not written yet
    public static void TestExampleV2(LoggerDelegate logger, ProgressDelegate progress, TokenDelegate tokencb)
    {
      Tools.freeHandle(oz_test_example_v2(forwardLogger(logger), forwardProgress(progress), forwardTokenV2(tokencb)));
    }
    */

    public static void TestExampleV3(LoggerDelegate logger, ProgressDelegate progress, TokenDelegate tokencb)
    {
      Tools.freeHandle(oz_test_example_v3(forwardLogger(logger), forwardProgress(progress), forwardTokenV3(tokencb)));
    }

    static internal InteropLoggerCallback forwardLogger(LoggerDelegate app_callback)
    {
      if (app_callback == null)
        return null;
      return (int level, string message) =>
      {
        try
        {
          //System.Console.WriteLine("C#:   ForwardLogger");
          return app_callback(level, message) ? 1 : 0;
        }
        catch
        {
          System.Console.WriteLine("Exception in C# logger callback");
          return 0; // Should not / did not print anything.
        }
      };
    }

    static internal InteropProgressCallback forwardProgress(ProgressDelegate app_callback)
    {
      if (app_callback == null)
        return null;
      return (long done, long total) =>
      {
        try
        {
          //System.Console.WriteLine("C#:   ForwardProgess");
          return app_callback(done, total) ? 1 : 0;
        }
        catch
        {
          System.Console.WriteLine("Exception in C# progress callback");
          return 1; // Continue execution, this is not fatal.
        }
      };
    }

    static internal InteropTokenCallbackV1 forwardTokenV1(TokenDelegate app_callback)
    {
      if (app_callback == null)
        return null;
      return () =>
      {
        try
        {
          //System.Console.WriteLine("C#:   ForwardTokenV1");
          return app_callback();
        }
        catch {
          System.Console.WriteLine("Exception in C# token callback V1");
          return ""; // Empty token.
        }
      };
    }

    /**
     * Strategy #3 for dealing with returned strings.
     * Ask the C++ layer to allocate memory for the returned string,
     * then copy the result into that unmanaged buffer.
     * Return the buffer without marshaling (i.e. declare as IntPtr)
     * and it will (a) be usable as-is cast to a char* and
     * (b) the C++ code now knows how to free it.
     * */
    static internal InteropTokenCallbackV3 forwardTokenV3(TokenDelegate app_callback)
    {
      if (app_callback == null)
        return null;
      return () =>
      {
        //System.Console.WriteLine("C#:   ForwardTokenV3");
        System.IntPtr hglobal_buffer = System.IntPtr.Zero;
        System.IntPtr malloc_buffer = System.IntPtr.Zero;
        try
        {
          string token = app_callback(); // might throw
          hglobal_buffer = Marshal.StringToHGlobalAnsi(token);
          if (hglobal_buffer != System.IntPtr.Zero)
          {
            malloc_buffer = oz_malloc(token.Length + 1);
            if (malloc_buffer != System.IntPtr.Zero)
            {
              // Both buffers are unmanaged memory, but it is not clear
              // whether it is safe to deallocate hglobal_buffer with free()
              // on the C++ side.
              for (int ii = 0; ii < token.Length; ++ii)
                Marshal.WriteByte(malloc_buffer, ii, Marshal.ReadByte(hglobal_buffer, ii));
              Marshal.WriteByte(malloc_buffer, token.Length, 0); // zero terminate
              return malloc_buffer; // Transfers ownership to the C++ side.
            }
          }
          System.Console.WriteLine("Allocation failure in C# token callback V3");
          return System.IntPtr.Zero;
        }
        catch
        {
          System.Console.WriteLine("Exception in C# token callback V3");
          return System.IntPtr.Zero;
        }
        finally
        {
          if (hglobal_buffer != System.IntPtr.Zero)
            Marshal.FreeHGlobal(hglobal_buffer);
          // Cannot get here with malloc_buffer empty.
        }
      };
    }
  }
}
