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

namespace OpenZGY.Managed
{
  public class ZgyBase : System.IDisposable
  {
    private static ZgyHandle empty_ = new ZgyHandle();
    internal protected ZgyHandle Handle { get; private set; } = null;

    // These are only used to manage object lifetimes. The C++ code
    // keeps references to the unmanaged std::function by copying it.
    // IOContext and ZgyWriterArgs are ephemeral, so the C# delegates
    // need to be saved in IOContext and ZgyWriterArgs.iocontext(),
    // and copied to Zgy{Reader,Writer,Utils} when the
    // ZgyWriterArgs or IOContext is used to create one of those.
    internal LoggerDelegate Logger { get; set; } = null;
    internal TokenDelegate TokenCB { get; set; } = null;

    protected ZgyBase(ZgyHandle handle)
    {
      Handle = handle;
    }

    public virtual void Dispose()
    {
      if (Handle != null)
      {
        var victim = Handle;
        Handle = null;
        victim.Dispose();
      }
    }

    /// <summary>
    /// Get the associated SafeHandle, converting nulls to an empty handle.
    /// </summary>
    internal static ZgyHandle GetHandle(ZgyBase instance)
    {
      return instance != null && instance.Handle != null ? instance.Handle : empty_;
    }
  }
}
