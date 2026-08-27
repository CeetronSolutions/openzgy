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

using System;
using System.Collections.Generic;
using System.Text;

namespace OpenZGY.Managed
{
  public class ProgressWithDots
  {
    private int dots_printed_;
    private int length_;
    private object lock_ = new object();

    /**
     * \brief Simple progress bar.
     */
    public ProgressWithDots(int length)
    {
      dots_printed_ = 0;
      length_ = length;
    }

    public ProgressDelegate get()
    {
      return (long gone, long total) =>
      {
        if (length_ < 1)
          return true;
        lock (lock_)
        {
          if (dots_printed_ == 0)
          {
            var ss = new System.Text.StringBuilder();
            ss.Append("[");
            for (int ii = 0; ii < length_; ++ii)
              ss.Append(" ");
            ss.Append("]\r[");
            System.Console.Write(ss.ToString());
          }
          long needed = (total <= 0) ? 1 : 1 + ((gone * (length_ - 1)) / total);
          if (needed > dots_printed_)
          {
            while (needed > dots_printed_)
            {
              System.Console.Write(".");
              dots_printed_ += 1;
            }
            System.Console.Out.Flush();
          }
          if (gone == total)
            System.Console.WriteLine();
          System.Console.Out.Flush();
          return true;
        }
      };
    }
  }
}
