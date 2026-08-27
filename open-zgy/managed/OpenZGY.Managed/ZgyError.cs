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

namespace OpenZGY.Managed.Errors
{
  public class ZgyError : System.ApplicationException
  {
    protected internal ZgyError(string message) : base(message)
    {
    }

    static public System.Exception MakeException(string type, string message)
    {
      switch (type)
      {
        case "ZgyNeedOldLibrary": return new ZgyNeedOldLibrary(message);
        case "ZgyUpdateRules": return new ZgyUpdateRules(message);
        case "ZgyFormatError": return new ZgyFormatError(message);
        case "ZgyCorruptedFile": return new ZgyCorruptedFile(message);
        case "ZgyUserError": return new ZgyUserError(message);
        case "ZgyInternalError": return new ZgyInternalError(message);
        case "ZgyEndOfFile": return new ZgyEndOfFile(message);
        case "ZgySegmentIsClosed": return new ZgySegmentIsClosed(message);
        case "ZgyAborted": return new ZgyAborted(message);
        case "ZgyMissingFeature": return new ZgyMissingFeature(message);
        case "ZgyIoError": return new ZgyIoError(message);
        case "ZgyNotReadOnlyError": return new ZgyNotReadOnlyError(message);
        case "ZgyError": return new ZgyError(message);
        default: return new System.ApplicationException(message);
      }
    }
  }

  public class ZgyFormatError : ZgyError
  {
    protected internal ZgyFormatError(string message) : base(message)
    {
    }
  }

  public class ZgyNeedOldLibrary : ZgyFormatError
  {
    protected internal ZgyNeedOldLibrary(string message) : base(message)
    {
    }
  }

  public class ZgyUpdateRules : ZgyFormatError
  {
    protected internal ZgyUpdateRules(string message) : base(message)
    {
    }
  }

  public class ZgyCorruptedFile : ZgyError
  {
    protected internal ZgyCorruptedFile(string message) : base(message)
    {
    }
  }

  public class ZgyUserError : ZgyError
  {
    protected internal ZgyUserError(string message) : base(message)
    {
    }
  }

  public class ZgyInternalError : ZgyError
  {
    protected internal ZgyInternalError(string message) : base(message)
    {
    }
  }

  public class ZgyEndOfFile : ZgyError
  {
    protected internal ZgyEndOfFile(string message) : base(message)
    {
    }
  }

  public class ZgySegmentIsClosed : ZgyError
  {
    protected internal ZgySegmentIsClosed(string message) : base(message)
    {
    }
  }

  public class ZgyAborted : ZgyError
  {
    protected internal ZgyAborted(string message) : base(message)
    {
    }
  }

  public class ZgyMissingFeature : ZgyError
  {
    protected internal ZgyMissingFeature(string message) : base(message)
    {
    }
  }

  public class ZgyIoError : ZgyError
  {
    protected internal ZgyIoError(string message) : base(message)
    {
    }
  }

  public class ZgyNotReadOnlyError : ZgyError
  {
    protected internal ZgyNotReadOnlyError(string message) : base(message)
    {
    }
  }
}
