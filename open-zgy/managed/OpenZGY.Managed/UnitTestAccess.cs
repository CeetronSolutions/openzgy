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
  /// <summary>
  /// Back door to allow unit tests to access selected internal methods.
  /// Only harmless stuff here. If there is a real risk of exposing something
  /// that application might actually use then consider ExternalsVisibkeTo.
  /// </summary>
  public static class UniTestAccess
  {
    /// <summary>
    /// Check lifetime of callback.
    /// </summary>
    /// <remarks>
    /// Make sure the logger that was specified in IOContext, possibly
    /// indirectly by setting an iocontext in a ZgyWriterArgs,
    /// has been copied all the way up to the Zgy{Reader.Writer,Utils}.
    /// It is important to this because a failure here might cause the
    /// delegate to be garbage collected early. Causing infrequent crashes.
    /// </remarks>
    public static bool LoggerIsSetTo(ZgyBase instance, LoggerDelegate logger)
    {
      return object.ReferenceEquals(instance.Logger, logger);
    }

    /// <summary>
    /// Check lifetime of callback.
    /// </summary>
    public static bool TokenCBIsSetTo(ZgyBase instance, TokenDelegate tokenCB)
    {
      return object.ReferenceEquals(instance.TokenCB, tokenCB);
    }
  }
}
