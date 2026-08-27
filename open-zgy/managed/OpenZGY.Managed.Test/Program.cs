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

namespace OpenZGY.Managed.Test
{
  internal class UnitTests
  {
    internal static void AssertEqual(object expect, object actual)
    {
      string message = null;
      if (expect == null && actual == null)
      {
      }
      else if (expect == null)
      {
        message = "Expected null, got something else.";
      }
      else if (actual == null)
      {
        message = "Expected something, got null.";
      }
      else if (object.ReferenceEquals(expect, actual))
      {
      }
      else if (expect.Equals(actual) && actual.Equals(expect))
      {
      }
      else
      {
        string exp = "(" + expect.GetType().Name + ")" + expect.ToString();
        string act = "(" + actual.GetType().Name + ")" + actual.ToString();
        message = $"Expected {exp} but got {act}";
      }
      if (!string.IsNullOrWhiteSpace(message))
      {
        System.Console.WriteLine(message);
        throw new System.ApplicationException("Unit test failure.");
      }
      else
      {
        //System.Console.WriteLine("Ok: {0}", actual);
      }
    }

    internal static void AssertEqualFloat(double expect, double actual, double epsilon)
    {
      if (System.Math.Abs(expect - actual) > epsilon)
      {
        System.Console.WriteLine($"Expected {expect} but got {actual}. Tolerance is {epsilon}, difference{expect-actual}.");
        throw new System.ApplicationException("Unit test failure.");
      }
    }

    internal static void AssertTrue(bool cond, string format, params object[] args)
    {
      if (!cond)
        throw new System.ApplicationException("Unit test failure. " + string.Format(format, args));
    }

    private static string cloud_temp_prefix_ = null;
    /// <summary>
    /// Prefix used for temporary files on the cloud.
    /// Within a single test run it returns the same value every time.
    /// </summary>
    internal static string CloudTempPrefix
    {
      get
      {
        if (cloud_temp_prefix_ == null)
        {
          // From native/src/test_utils.cpp; should work the same way because
          // there are tools that recognize and purge any old garbage.
          //ss << "/tmp-"
          //   << std::hex << std::setw(8) << std::setfill('0')
          //   << std::uint32_t(time(nullptr))
          //   << "-";
          string project = System.Environment.GetEnvironmentVariable("OPENZGY_SDTESTSINK");
          if (string.IsNullOrWhiteSpace(project))
            project = "sd://opendes/slb-testsink/a";
          long random = System.DateTimeOffset.UtcNow.ToUnixTimeSeconds();
          cloud_temp_prefix_ = $"{project}/tmp-{random:x8}-";
          //System.Console.WriteLine($"Cloud prefix is {cloud_temp_prefix_}");
        }
        return cloud_temp_prefix_;
      }
    }

    /**
     * Check almost all of SeismicStoreIOContext by creating an instancs
     * (setting one property at a time) and then dumping the entire struct
     * (formatted on the C++ side, returned as a multi-line string).
     *
     * Not covered here: logger, sdtokencb, retryCount.
     * The last one is just missing from toString().
     * logger and sdtokencb will be tested elsewhere.
     *
     * Covered:
     * aligned clone cputhreads create forceRoBeforeRead forceRwBeforeWrite
     * iothreads legaltag maxhole maxsize sdapikey sdtoken sdurl segsize
     * segsplit seismicmeta setRoAfterWrite toString writeid writethreads
     */
    internal static void TestSDContext()
    {
      // ARRANGE
      string expect1 =
"SD context:\n" +
"  sdurl:    \"sd-url\"\n" +
"  sdapikey: \"sd-apikey\"\n" +
"  sdtoken:  \"sd-token\"\n" +
"  maxsize:  3 MB\n" +
"  maxhole:  4 MB\n" +
"  aligned:  512 MB\n" +
"  segsize:  448 MB\n" +
"  segsplit: 7\n" +
"  threads:  4 I/O, 5 CPU, 1 Write\n" +
"  legaltag: \"my-legaltag\"\n" +
"  writeid:  \"my-writeid\"\n" +
"  seismeta: \"my-seismicmeta\"\n" +
"  romode:   ro_after_writero_before_readrw_before_write\n"
+ "";
      string expect2 = expect1.Replace("sd-token", "updated-token");

      // ACT
      var context = new SeismicStoreIOContext()
        .sdurl("sd-url")
        .sdapikey("sd-apikey")
        .sdtoken("sd-token")
        //.sdtokencb(callback)
        .maxsize(3)
        .maxhole(4)
        .aligned(512)
        .segsize(64)
        .segsplit(7)
        .iothreads(4)
        .cputhreads(5)
        .writethreads(1)
        .legaltag("my-legaltag")
        .writeid("my-writeid")
        .seismicmeta("my-seismicmeta")
        .setRoAfterWrite(true)
        .forceRoBeforeRead(true)
        .forceRwBeforeWrite(true)
        .retryCount(3)
        .logger(Callbacks.EmptyCallback())
      ;
      var context2 = context.clone().sdtoken("updated-token");

      //System.Console.WriteLine(context.toString());
      //System.Console.WriteLine(context2.toString());

      // ASSERT
      try
      {
        string actual1 = context.toString();
        string actual2 = context2.toString();
        AssertEqual(expect1, actual1);
        AssertEqual(expect2, actual2);
      }
      finally
      {
        // Disposing these is not critical, unlike ZgyReader, ZgyUtils,
        // and in particular ZgyWriter. But it is a good habit where
        // unmanaged handles are involved.
        if (context2 != null)
          context2.Dispose();
        if (context != null)
          context.Dispose();
      }
    }

    /**
     * Same strategy as TestContext.
     *
     * Not covered here: iocontext, metafrom.
     *
     * Covered:
     * bricksize corners create datarange datatype filename hunit
     * ilinc ilstart merge size toString xlinc xlstart zfp_compressor
     * zfp_lodcompressor zinc zstart zunit
     */
    internal static void TestZgyWriterArgs()
    {
      var args = new ZgyWriterArgs()
  .filename("my pretty file")
  .size(111, 222, 333)
  .bricksize(32, 32, 128)
  .datatype(SampleDataType.int16)
  .datarange(-1000, +1200)
  .zunit(UnitDimension.time, "ms", 0.0001)
  .hunit(UnitDimension.length, "ft", 0.3048)
  .ilstart(101).ilinc(2)
  .xlstart(500).xlinc(5)
  .zstart(-100).zinc(4)
  .zfp_compressor(30)
  .zfp_lodcompressor(50)
  .corners(new double[][] {
    new double[] { 91, 92 },
    new double[] { 93, 94 },
    new double[] { 95, 96 },
    new double[] { 97, 98 }
  });
      var args2 = new ZgyWriterArgs().filename("two").merge(args);
      // Next is a no-op because merge isn't supoosed to merge defaults.
      args.merge(new ZgyWriterArgs());

      string expect0 = @"ZgyWriterArgs
  filename:    """"
  iocontext:   (null)
  compressor:  (null)
  lodcompress: (null)
  size:        (0,0,0)
  bricksize:   (64,64,64)
  datatype:    1003
  datarange:   0 to -1
  zunit:       2000 """" 1
  hunit:       2000 """" 1
  ilstart/inc: 0 / 0
  xlstart/inc: 0 / 0
  zstart/inc:  0 / 0
  corner0:     0, 0
  corner1:     0, 0
  corner2:     0, 0
  corner3:     0, 0
".Replace("\r\n", "\n");

      string expect1 = @"ZgyWriterArgs
  filename:    ""my pretty file""
  iocontext:   (null)
  compressor:  *
  lodcompress: *
  *size:        (111,222,333)
  *bricksize:   (32,32,128)
  *datatype:    1002
  *datarange:   -1000 to 1200
  *zunit:       2001 ""ms"" 0.0001
  *hunit:       2002 ""ft"" 0.3048
  *ilstart/inc: 101 / 2
  *xlstart/inc: 500 / 5
  *zstart/inc:  -100 / 4
  *corner0:     91, 92
  *corner1:     93, 94
  *corner2:     95, 96
  *corner3:     97, 98
".Replace("\r\n", "\n");

      string expect2 = @"ZgyWriterArgs
  filename:    ""two""
  iocontext:   (null)
  compressor:  (null)
  lodcompress: (null)
  *size:        (111,222,333)
  *bricksize:   (32,32,128)
  *datatype:    1002
  *datarange:   -1000 to 1200
  *zunit:       2001 ""ms"" 0.0001
  *hunit:       2002 ""ft"" 0.3048
  *ilstart/inc: 101 / 2
  *xlstart/inc: 500 / 5
  *zstart/inc:  -100 / 4
  *corner0:     91, 92
  *corner1:     93, 94
  *corner2:     95, 96
  *corner3:     97, 98
".Replace("\r\n", "\n");

      AssertEqual(expect0, new ZgyWriterArgs().toString());
      AssertEqual(expect1, args.toString());
      AssertEqual(expect2, args2.toString());
    }

    /**
     * Check that if a logger- or token callback is provided in an IOContext,
     * which in turn is provided directly in a ZgyReader or ZgyUtils, or
     * via a ZgyWriterArgs instance passed to a ZgyWriter, then those
     * delegates will be referencend from the ZgyReader/Writer/Utils on the
     * mamaged side. This is to keep the managed delegates stay in scope
     * as long as the unmanaged ones might be retained at the C++ level.
     *
     * The progress callback does not have the same issues because its
     * scope is only a single function.
     *
     * The test ought to work both in a seismic store context and on-prem,
     * by using the kludge of passing a SeismicStoreIOContext also to an
     * on-prem writer.
     */
    internal static void TestCallbackLifetime(string filename)
    {
      var token = System.Environment.GetEnvironmentVariable("OPENZGY_TOKEN");
      var logger = new LoggerDelegate((int level, string message) => false);
      var tokencb = new TokenDelegate(() => token);
      var context = new SeismicStoreIOContext().logger(logger).sdtokencb(tokencb);
      var args = new ZgyWriterArgs().iocontext(context).filename(filename).size(4, 4, 4);

      if (!UniTestAccess.LoggerIsSetTo(context, logger))
        throw new System.ApplicationException("IOContext did not capture logger delegate");
      if (!UniTestAccess.TokenCBIsSetTo(context, tokencb))
        throw new System.ApplicationException("IOContext did not capture tokencb delegate");

      if (!UniTestAccess.LoggerIsSetTo(args, logger))
        throw new System.ApplicationException("ZgyWriterArgs did not capture logger delegate");
      if (!UniTestAccess.TokenCBIsSetTo(args, tokencb))
        throw new System.ApplicationException("ZgyWriterArgs did not capture tokencb delegate");

      // For local files the IOContext is not actually used, but used "enough" that we can test.
      using (ZgyWriter writer = ZgyWriter.open(args))
      {
        if (!UniTestAccess.LoggerIsSetTo(writer, logger))
          throw new System.ApplicationException("ZgyWriter did not capture logger delegate");
        if (!UniTestAccess.TokenCBIsSetTo(writer, tokencb))
          throw new System.ApplicationException("ZgyWriter did not capture tokencb delegate");
        writer.close_incomplete();
      }

      using (ZgyReader reader = ZgyReader.open(filename, context))
      {
        if (!UniTestAccess.LoggerIsSetTo(reader, logger))
          throw new System.ApplicationException("ZgyReader did not capture logger delegate");
        if (!UniTestAccess.TokenCBIsSetTo(reader, tokencb))
          throw new System.ApplicationException("ZgyReader did not capture tokencb delegate");
      }

      using (var utils = ZgyUtils.utils("/", context))
      {
        if (!UniTestAccess.LoggerIsSetTo(utils, logger))
          throw new System.ApplicationException("ZgyUtils did not capture logger delegate");
        if (!UniTestAccess.TokenCBIsSetTo(utils, tokencb))
          throw new System.ApplicationException("ZgyUtils did not capture tokencb delegate");
        utils.deleteFile(filename, false);
      }
    }

    internal static void TestCreateFile(string filename)
    {
      var log = new System.Text.StringBuilder();
      long done = -1;
      long total = -1;
      int ping = 0;
      var logger = new LoggerDelegate((int level, string message) =>
      {
        if (!string.IsNullOrWhiteSpace(message))
          log.AppendFormat("{0}: {1}{2}", level, message,
            message.EndsWith("\n") ? "" : System.Environment.NewLine);
        return true;
      });
      var progress = new ProgressDelegate((long new_done, long new_total) =>
      {
        done = new_done;
        total = new_total;
        ++ping;
        return true;
      });
      using (var iocontext = new SeismicStoreIOContext().logger(logger))
      using (var args = new ZgyWriterArgs()
          .filename(filename)
          .size(111, 222, 333)
          .bricksize(32, 32, 128)
          .datatype(SampleDataType.int16)
          .datarange(-1000, +1200)
          .zunit(UnitDimension.time, "ms", 0.0001)
          .hunit(UnitDimension.length, "ft", 0.3048)
          .ilstart(101).ilinc(2)
          .xlstart(500).xlinc(5)
          .zstart(-100).zinc(4)
          .iocontext(iocontext)
          .corners(new double[][] {
           new double[] { 10000, 20000 },
           new double[] { 11000, 20000 },
           new double[] { 10000, 22000 },
           new double[] { 11000, 22000 }
          }))
      using (ZgyWriter writer = ZgyWriter.open(args))
      {
        writer.writeconst(new long[3] { 0, 0, 0 }, writer.size, new float[] { 42.0f });
        writer.write(new long[3] { 0, 0, 0 }, new long[3] { 1, 1, 5 }, new float[5] { 47, 51, 41, 00, 98 });
        writer.write(new long[3] { 0, 1, 0 }, new long[3] { 1, 1, 5 }, new short[5] { -500, 0, 1000, 4000, 32760 });
        writer.finalize(new DecimationType[] { DecimationType.Average }, progress, FinalizeAction.BuildDefault, false);
        writer.close();
      }
      //System.Console.WriteLine($"Written {done}/{total}, stats reported {ping} times");
      //System.Console.WriteLine(log.ToString());
      AssertTrue(log.ToString().Contains("open for write:"), "Logger doesn't seem to work.");
      AssertTrue(done > 0 && done == total && ping > 0, $"Written {done}/{total}, stats reported {ping} times");
    }

    internal static void TestReadFile(string filename)
    {
      using (ZgyReader reader = ZgyReader.open(filename, null))
      {
        var start = new long[3] { 0, 0, 0 };
        var size = new long[3] { 1, 2, 5 };
        var samplecount = size[0] * size[1] * size[2];
        var bulk = new float[samplecount];
        double const_value = -999.25;
        bool is_const = reader.readconst(ref const_value, start, size, 0, true);
        AssertTrue(!is_const, "The data at the start of the cube should not be const.");
        var next = new long[3] { 0, 0, 128 };
        bool is_next_const = reader.readconst(ref const_value, next, size, 0, true);
        AssertTrue(is_next_const, "The data in the next brick SHOULD be const.");
        AssertTrue(System.Math.Abs(const_value - 42) < 0.01, $"Const value should have been 42, was {const_value}");
        reader.read(start, size, bulk, 0);
        var expect = new float[] { 47, 51, 41, 00, 98, 83.23f, 100.02f, 133.59f, 234.30f, 1199.77f };
        for (int ii = 0; ii < samplecount; ++ii)
          AssertTrue(System.Math.Abs(bulk[ii] - expect[ii]) <= 0.05, $"bulk data (float) mismatch at {ii} expected {expect[ii]} actual {bulk[ii]}");
        var ibulk = new short[samplecount];
        reader.read(start, size, ibulk, 0);
        var iexpect = new short[] { -1579, -1460, -1758, -2979, -60, -500, 0, 1000, 4000, 32760 };
        for (int ii = 0; ii < samplecount; ++ii)
          AssertTrue(ibulk[ii] == iexpect[ii], $"bulk data (short) mismatch at {ii} expected {iexpect[ii]} actual {ibulk[ii]}");
      }
    }

    /**
     * Read from a file open for write.
     *
     * The ZgyWriter has read() and readconst() methods but they are not
     * the same as those in the reader. Also they lack the lod parameter.
     */
    internal static void TestReadWritableFile(string filename)
    {
      using (var args = new ZgyWriterArgs().filename(filename))
      {
        using (ZgyWriter writer = ZgyWriter.reopen(args))
        {
          //System.Console.WriteLine(writer.toString());
          //System.Console.WriteLine(writer.filestats().toString());
          var start = new long[3] { 0, 0, 0 };
          var size = new long[3] { 1, 2, 5 };
          var samplecount = size[0] * size[1] * size[2];
          var bulk = new float[samplecount];
          double const_value = -999.25;
          bool is_const = writer.readconst(ref const_value, start, size, /*0,*/ true);
          AssertTrue(!is_const, $"The data at the start of the cube should not be const {const_value}.");
          var next = new long[3] { 0, 0, 128 };
          bool is_next_const = writer.readconst(ref const_value, next, size, /*0,*/ true);
          AssertTrue(is_next_const, "The data in the next brick should be const.");
          AssertTrue(System.Math.Abs(const_value - 42) < 0.01, $"Const value should have been 42, was {const_value}");
          writer.read(start, size, bulk/*, 0*/);
          var expect = new float[] { 47, 51, 41, 00, 98, 83.23f, 100.02f, 133.59f, 234.30f, 1199.77f };
          for (int ii = 0; ii < samplecount; ++ii)
            AssertTrue(System.Math.Abs(bulk[ii] - expect[ii]) <= 0.05, $"bulk data (float) mismatch at {ii} expected {expect[ii]} actual {bulk[ii]}");
          var ibulk = new short[samplecount];
          writer.read(start, size, ibulk/*, 0*/);
          var iexpect = new short[] { -1579, -1460, -1758, -2979, -60, -500, 0, 1000, 4000, 32760 };
          for (int ii = 0; ii < samplecount; ++ii)
            AssertTrue(ibulk[ii] == iexpect[ii], $"bulk data (short) mismatch at {ii} expected {iexpect[ii]} actual {ibulk[ii]}");
          AssertEqual(1, writer.nlods);
        }
      }
    }

    internal static void CopyFile(string input, string output, bool use_token_callback)
    {
      using (SeismicStoreIOContext context = new SeismicStoreIOContext()
        .sdurl(System.Environment.GetEnvironmentVariable("OPENZGY_SDURL"))
        .sdapikey(System.Environment.GetEnvironmentVariable("OPENZGY_SDAPIKEY"))
        .iothreads(4).cputhreads(8).retryCount(4)
        /*.logger(Callbacks.StandardCallback(2, "logger callback: "))*/)
      {
        if (use_token_callback)
          context.sdtokencb((() => System.Environment.GetEnvironmentVariable("OPENZGY_TOKEN")));
        else
          context.sdtoken(System.Environment.GetEnvironmentVariable("OPENZGY_TOKEN"));
        using (ZgyReader reader = OpenZGY.Managed.ZgyReader.open(input, context))
        {
          using (ZgyWriterArgs args = new ZgyWriterArgs().metafrom(reader).filename(output).iocontext(context))
          {
            using (ZgyWriter writer = ZgyWriter.open(args))
            {
              long samples = reader.size[0] * reader.size[1] * reader.size[2];
              float[] data = new float[samples];
              var origin = new long[3] { 0, 0, 0 };
              reader.read(origin, reader.size, data, 0);
              writer.write(origin, reader.size, data);
              writer.finalize(new DecimationType[1] { DecimationType.Average }, null, FinalizeAction.BuildDefault, false);
              writer.close();
            }
          }
        }
      }
    }

    private static bool FileExists(string filename)
    {
      try
      {
        using (var context = new SeismicStoreIOContext())
        {
          using (ZgyReader reader = OpenZGY.Managed.ZgyReader.open(filename, context))
          {
          }
          return true;
        }
      }
      catch (OpenZGY.Managed.Errors.ZgyIoError)
      {
        return false;
      }
      catch (OpenZGY.Managed.Errors.ZgyInternalError)
      {
        // The cloud accessor might throw this exception.
        // It probably shouldn't.
        return false;
      }
    }

    internal static void TestDeleteFile(string filename)
    {
      AssertTrue(FileExists(filename), $"{filename} did not exist");
      using (var context = new SeismicStoreIOContext())
      {
        using (var utils = ZgyUtils.utils(filename, context))
        {
          utils.deleteFile(filename, false);
        }
      }
      AssertTrue(!FileExists(filename), $"{filename} did not get deleted");
    }

    /**
     * Testing utils_alturl.
     * Only the API is being tested here, so run altUrl on a
     * regular file where it will always return its input.
     */
    internal static void TestAltUrl(string filename)
    {
      string alternate;
      using (var utils = ZgyUtils.utils("/", null))
      {
        alternate = utils.alturl(filename);
        //System.Console.WriteLine($"Alturl: {alternate}");
      }
      AssertEqual(filename, alternate);
    }

    /**
     * This test should work also if the token is expired.
     * It might not work for a bogus string, see
     * CallbackAuthProvider::getServiceAuthTokenImpl() and
     * utils:getAuthTokenExpiration()
     */
    internal static void TestIdToken(string filename)
    {
      var token = System.Environment.GetEnvironmentVariable("OPENZGY_TOKEN");
      //var token = "FakeTOKEN";
      var context = new SeismicStoreIOContext().sdtoken(token);
      string token2;
      using (var utils = ZgyUtils.utils("sd://", context))
      {
        token2 = utils.idtoken(filename);
        //System.Console.WriteLine($"IdToken: {token}");
      }
      AssertEqual(token, token2);
    }

    internal static void TestCredentialsFrom(string filename)
    {
      using (var context = new SeismicStoreIOContext()
        .sdurl(System.Environment.GetEnvironmentVariable("OPENZGY_SDURL"))
        .sdapikey(System.Environment.GetEnvironmentVariable("OPENZGY_SDAPIKEY"))
        .sdtoken(System.Environment.GetEnvironmentVariable("OPENZGY_TOKEN")))
      {
        using (var utils = ZgyUtils.utils(filename, context))
        {
          using (var shared = new SeismicStoreIOContext().credentialsFrom(utils))
          {
            using (ZgyReader reader = OpenZGY.Managed.ZgyReader.open(filename, shared))
            {
              AssertEqual((long)111, reader.size[0]);
              AssertEqual((long)222, reader.size[1]);
              AssertEqual((long)333, reader.size[2]);
            }
          }
        }
      }
    }

    private static double Norm(double x, double y)
    {
      return System.Math.Sqrt(x * x + y * y);
    }

    private static double Delta(double[] a, double[] b)
    {
      return Norm(a[0] - b[0], a[1] - b[1]);
    }

    /**
     * Test ZgyTools (parent of ZgyReader and ZgyWriter).
     * Covers all 6 converson methods. The general transform
     * is not available in the managed API.
     */
    internal static void TestCoordinateConversion(string filename)
    {
      using (ZgyReader reader = OpenZGY.Managed.ZgyReader.open(filename, null))
      {
        var ipoint = new double[2] { 0, 0 };
        var apoint = reader.indexToAnnot(ipoint);
        var wpoint = reader.annotToWorld(apoint);
        var apoint2 = reader.worldToAnnot(wpoint);
        var ipoint2 = reader.annotToIndex(apoint2);
        var wpoint2 = reader.indexToWorld(ipoint2);
        var ipoint3 = reader.worldToIndex(wpoint2);
        AssertTrue(Delta(ipoint, ipoint2) <= 0.0001, "First coordinamne conversion round trip failed.");
        AssertTrue(Delta(ipoint, ipoint3) <= 0.0001, "First coordinamne conversion round trip failed.");
        //System.Console.WriteLine($"Round trip: ({ipoint3[0]}, {ipoint3[1]}) annot ({apoint[0]}, {apoint[1]}) world ({wpoint[0]:f2}, {wpoint[1]:f2})");
      }
    }

    /// <summary>
    /// Output most of a file's metadata.
    /// Formatting done in the C++ code.
    /// </summary>
    internal static string CXXDumpMeta(ZgyReader r)
    {
      return r.toString();
    }

    /// <summary>
    /// Output most of a file's metadata.
    /// Tries to make the output as close as possible to that from
    /// ZgyReader.toString() so the two can later be diffed.
    /// </summary>
    internal static string ManagedDumpMeta(ZgyReader r)
    {
      //System.Console.WriteLine("\nC++ dump{0}\n", r.toString());
      using (SampleStatistics stat = r.statistics())
      using (SampleHistogram hist = r.histogram())
      using (FileStatistics filestats = r.filestats())
      {
        var icorners = r.indexcorners;
        var acorners = r.annotcorners;
        var wcorners = r.corners;
        var nl = System.Environment.NewLine;
        var ss = new System.Text.StringBuilder();
        ss.AppendFormat("File format and version        = SampleDataType::{1} ZGY version {2}{0}", nl, r.datatype, filestats.fileVersion);
        ss.AppendFormat("Current data Version           = {1}{0}", nl, r.verid);
        ss.AppendFormat("Size I,J,K                     = ({1}, {2}, {3}){0}", nl, r.size[0], r.size[1], r.size[2]);
        ss.AppendFormat("Brick size I,J,K               = ({1}, {2}, {3}){0}", nl, r.bricksize[0], r.bricksize[1], r.bricksize[2]);
        ss.AppendFormat("Number of bricks I,J,K         = ({1}, {2}, {3}){0}", nl,
          (r.size[0] + r.bricksize[0] - 1) / r.bricksize[0],
          (r.size[1] + r.bricksize[1] - 1) / r.bricksize[1],
          (r.size[2] + r.bricksize[2] - 1) / r.bricksize[2]);
        ss.AppendFormat("Number of LODs                 = {1}{0}", nl, r.nlods);
        ss.AppendFormat("Coding range min/max           = {1:g6} {2:g6} (raw: {3:g6} {4:g6}) {5}{0}", nl,
            r.datarange[0], r.datarange[1],
            r.raw_datarange[0], r.raw_datarange[1],
            r.size[0] * r.size[1] * r.size[2]);
        ss.AppendFormat("Statistical min/max/count/avg  = {1:g6} {2:g6} {3} {4:g6} {5:g6}{0}", nl, stat.min, stat.max, stat.cnt, stat.sum / stat.cnt, (System.Math.Sqrt(stat.ssq / stat.cnt)));
        ss.AppendFormat("Histogram range min/max/count  = {1:g6} {2:g6} {3} bincount {4}{0}", nl, hist.minvalue, hist.maxvalue, hist.samplecount, hist.bins().Length);
        ss.AppendFormat("Inline start/increment/count   = {1:g6} {2:g6} {3}{0}", nl, r.annotstart[0], r.annotinc[0], r.size[0]);
        ss.AppendFormat("Xline  start/increment/count   = {1:g6} {2:g6} {3}{0}", nl, r.annotstart[1], r.annotinc[1], r.size[1]);
        ss.AppendFormat("Sample start/increment/count   = {1:g6} {2:g6} {3}{0}", nl, r.zstart, r.zinc, r.size[2]);
        ss.AppendFormat("Horizontal dim/factor/name     = UnitDimension::{1} {2} '{3}'{0}", nl, r.hunitdim, r.hunitfactor, r.hunitname);
        ss.AppendFormat("Vertical dim/factor/name       = UnitDimension::{1} {2} '{3}'{0}", nl, r.zunitdim, r.zunitfactor, r.zunitname);
        ss.AppendFormat("Ordered Corner Points Legend   = {1}{0}", nl, "[  <i>,   <j>] { <inline>,   <xline>} (  <easting>,  <northing>)");
        for (int ix = 0; ix < 4; ++ix)
        {
          ss.AppendFormat("Ordered Corner Point {7}         = [{0,5}, {1,5}] {{{2,9}, {3,9}}} ({4,11:f2}, {5,11:f2}){6}",
            icorners[ix][0], icorners[ix][1],
            acorners[ix][0], acorners[ix][1],
            wcorners[ix][0], wcorners[ix][1],
            System.Environment.NewLine, ix);
        }
        return ss.ToString().Replace(System.Environment.NewLine, "\n");
      }
    }

    /**
     * This test gives good coverage of metadata access
     * but it is extremely fragile.
     */
    internal static void CompareDumpMeta(ZgyReader reader)
    {
        string cxx = CXXDumpMeta(reader);
        string local = ManagedDumpMeta(reader);
        if (cxx != local)
        {
          System.Console.WriteLine($"\nC++ dump:\n{cxx}");
          System.Console.WriteLine($"\nManaged dump:\n{local}");
          int len = System.Math.Min(cxx.Length, local.Length);
          for (int ii = 0; ii < len; ++ii)
          {
            if (cxx[ii] != local[ii])
            {
              System.Console.WriteLine($"Strings differ at position {ii}");
            }
          }
          AssertTrue(cxx.Length == local.Length, "Mismatch between dump from managed and C++ code. Length {0}/{1}", cxx.Length, local.Length);
          AssertTrue(cxx == local, "Mismatch between dump from managed and C++ code.");
        }
      }

    /**
     * As CompareDumpMeta but assert the contents of one property at a time.
     * Given that the coverage is the same this makes CompareDumpMeta redundant.
     */
    internal static void TestVerifyMeta(string filename)
    {
      using (var r = ZgyReader.open(filename, null))
      using (SampleStatistics stat = r.statistics())
      using (SampleHistogram hist = r.histogram())
      using (FileStatistics filestats = r.filestats())
      {
        var icorners = r.indexcorners;
        var acorners = r.annotcorners;
        var wcorners = r.corners;
        AssertEqual(SampleDataType.int16, r.datatype);
        AssertTrue(!string.IsNullOrWhiteSpace(r.verid), "Empty verid"); // UUID cganges each time.
        AssertEqual((long)111, r.size[0]);
        AssertEqual(222L, r.size[1]);
        AssertEqual(333L, r.size[2]);
        AssertEqual(32L, r.bricksize[0]);
        AssertEqual(32L, r.bricksize[1]);
        AssertEqual(128L, r.bricksize[2]);
        AssertEqual(4, r.nlods);
        AssertEqualFloat(-1000.0, r.datarange[0], 0.01);
        AssertEqualFloat(+1200.0, r.datarange[1], 0.01);
        AssertEqualFloat(-1000.0, r.raw_datarange[0], 0.01);
        AssertEqualFloat(+1200.0, r.raw_datarange[1], 0.01);
        AssertEqualFloat(101, r.annotstart[0], 0.01);
        AssertEqualFloat(2, r.annotinc[0], 0.01);
        AssertEqualFloat(500, r.annotstart[1], 0.01);
        AssertEqualFloat(5, r.annotinc[1], 0.01);
        AssertEqualFloat(-100, r.zstart, 0.01);
        AssertEqualFloat(4, r.zinc, 0.001);
        AssertEqual(UnitDimension.length, r.hunitdim);
        AssertEqualFloat(0.3048, r.hunitfactor, 0.0001);
        AssertEqual("ft", r.hunitname);
        AssertEqual(UnitDimension.time, r.zunitdim);
        AssertEqualFloat(0.0001, r.zunitfactor, 0.0000001);
        AssertEqual("ms", r.zunitname);
        // Statistics()
        AssertEqualFloat(0.0122072, stat.min, 0.0001);
        AssertEqualFloat(1199.77, stat.max, 0.05);
        AssertEqual(8205786L, stat.cnt);
        AssertEqualFloat(42.0083, stat.sum / stat.cnt, 0.0002);
        AssertEqualFloat(42.0103, (System.Math.Sqrt(stat.ssq / stat.cnt)), 0.0002);
        // Histogram
        AssertEqualFloat(-1000.0, hist.minvalue, 0.01);
        AssertEqualFloat(+1200.0, hist.maxvalue, 0.01);
        AssertEqual(r.size[0] * r.size[1] * r.size[2], hist.samplecount);
        AssertTrue(hist.bins() != null, "Histogram bins are missing");
        AssertEqual(256, hist.bins().Length);
        // Lattice
        double[][] actual_index = r.indexcorners;
        double[][] actual_annot = r.annotcorners;
        double[][] actual_world = r.corners;
        double[][] expect_index = new double[4][] {
           new double[2] { 0, 0 },
           new double[2] { r.size[0] - 1, 0 },
           new double[2] { 0, r.size[1] - 1 },
           new double[2] { r.size[0] - 1, r.size[1] - 1 }
        };
        double[][] expect_annot = new double[4][] {
           new double[2] { 101, 500 },
           new double[2] { 321, 500 },
           new double[2] { 101, 1605 },
           new double[2] { 321, 1605 }
        };
        double[][] expect_world = new double[4][] {
           new double[2] { 10000, 20000 },
           new double[2] { 11000, 20000 },
           new double[2] { 10000, 22000 },
           new double[2] { 11000, 22000 }
        };
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 2; ++jj)
            AssertEqualFloat(expect_index[ii][jj], actual_index[ii][jj], 0.0001);
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 2; ++jj)
            AssertEqualFloat(expect_annot[ii][jj], actual_annot[ii][jj], 0.01);
        for (int ii = 0; ii < 4; ++ii)
          for (int jj = 0; jj < 2; ++jj)
            AssertEqualFloat(expect_world[ii][jj], actual_world[ii][jj], 0.01);
        // Filestats
        long samples_per_brick = r.bricksize[0] * r.bricksize[1] * r.bricksize[2];
        //System.Console.WriteLine(filestats.toString());
        AssertEqual(0L, filestats.alphaCcompressedSize);
        AssertEqual(0L, filestats.alphaCompressedCount);
        AssertEqual(0L, filestats.alphaConstantCount);
        AssertEqual(39L, filestats.alphaMissingCount);
        AssertEqual(0L, filestats.alphaNormalCount);
        AssertEqual(4L, filestats.brickNormalCount);
        AssertEqual(2 * samples_per_brick, filestats.brickNormalSizePerEntry);
        AssertEqualFloat(1.0, filestats.compressionFactor, 0.01);
        AssertEqual(2 * samples_per_brick, filestats.dataStart);
        AssertEqual(10 * samples_per_brick, filestats.fileSize);
        AssertEqual(3L, filestats.fileVersion);
        AssertEqual(3555L, filestats.headerSize);
        AssertTrue(!filestats.isCompressed, "should not be compressed");
        AssertEqual(1052131L, filestats.usedIfUncompressed);
        AssertEqual(1052131L, filestats.usedSize);
      }
    }

    internal static void TestCompareDumpMeta(string filename)
    {
      using (var reader = ZgyReader.open(filename, null))
      {
        CompareDumpMeta(reader);
      }
    }

    /**
     * This represents a function in the application that needs to
     * call an OpenZGY method expecting one or more callbacks.
     *
     * NOTE! If the C++ layer keeps a reference to the callback
     * (which it doesn't do here, but might for logger and token)
     * then the caller is responsible for keeping the delegate
     * on the C# side alive.
     *
     * NOTE! Callbacks need to be thread safe.
     */
    internal static void TestExampleCallbacks(int version)
    {
      int loggerCB_calls = 0;
      LoggerDelegate loggerCB = (int level, string message) =>
      {
        ++loggerCB_calls;
        //System.Console.WriteLine("C#:   Logger({0}, \"{1}\")", level, message);
        return true;
      };

      int progressCB_calls = 0;
      ProgressDelegate progressCB = (long gone, long total) =>
      {
        ++progressCB_calls;
        //System.Console.WriteLine("C#:   Progress {0}/{1}", gone, total);
        return true;
      };

      int tokenCB_calls = 0;
      TokenDelegate tokenCB = () =>
      {
        ++tokenCB_calls;
        //System.Console.WriteLine("C#:   Requested token");
        return "Hello, here is your token!";
      };

      switch (version)
      {
        case 1: Callbacks.TestExampleV1(loggerCB, progressCB, tokenCB); break;
        case 3: Callbacks.TestExampleV3(loggerCB, progressCB, tokenCB); break;
        default:
          System.Console.WriteLine("Uknown version {0} of TestExampleCallbacks.", version);
          break;
      }

      AssertTrue(loggerCB_calls > 0, "Logger callback was not invoked.");
      AssertTrue(progressCB_calls > 0, "Progress callback was not invoked.");
      AssertTrue(tokenCB_calls > 0, "Token callback was not invoked.");

      // Possibly redundant in this case but it doesn't hurt.
      // In the real code the managed delegates need to be
      // referenced by the ZgyBase class so they don't go out
      // of scope while the C++ API still has a reference.
      System.GC.KeepAlive(loggerCB);
      System.GC.KeepAlive(progressCB);
      System.GC.KeepAlive(tokenCB);
    }
  }

  class Program
  {
    /// <summary>
    /// Dump the first few bytes on the file.
    /// </summary>
    public static void DumpBulk(ZgyReader reader)
    {
      var start = new long[3] { 0, 0, 6 };
      var size = new long[3] { 2, 2, 10 };
      var bulk = new float[size[0] * size[1] * size[2]];
      double const_value = -999.25;
      bool is_const = reader.readconst(ref const_value, start, size, 0, true);
      reader.read(start, size, bulk, 0);
      var sb = new System.Text.StringBuilder();
      if (is_const)
        sb.AppendLine($"Region is constant value {const_value:f1}.");
      else
        sb.AppendLine("Region is not known to be constant.");
      for (int ii = 0; ii < 40; ii += 10)
      {
        sb.Append("   ");
        for (int jj = 0; jj < 10; ++jj)
          sb.AppendFormat("{0,8:f1}", bulk[ii + jj]);
        sb.Append(System.Environment.NewLine);
      }
      System.Console.Write(sb.ToString());
    }

    private static void Hello(string message)
    {
      //System.Console.WriteLine("......................................................................");
      System.Console.WriteLine(">>> " + message);
      //System.Console.WriteLine("......................................................................");
    }

    static void Main(string[] args)
    {
      //System.Console.WriteLine("Starting up, please attach a debugger and press return");
      //System.Console.ReadLine();
      //System.Console.WriteLine("Thank you.");
      bool have_sd = true;
      foreach (var it in args)
        if (it == "--nosd")
          have_sd = false;
      string guid = System.Guid.NewGuid().ToString();
      string local_file = System.IO.Path.GetTempPath() + "OpenZGY.Managed.Test." + guid + ".zgy";
      string cloud_file = have_sd ? UnitTests.CloudTempPrefix + "OpenZGY.Managed.Test.zgy" : null;
      try
      {
        Hello("TestSDContext");
        UnitTests.TestSDContext();
        Hello("TestZgyWriterArgs");
        UnitTests.TestZgyWriterArgs();
        Hello("TestCallbackLifetime");
        UnitTests.TestCallbackLifetime(local_file);
        Hello("TestCreateFile");
        UnitTests.TestCreateFile(local_file);
        Hello("TestReadFile");
        UnitTests.TestReadFile(local_file);
        Hello("TestReadWritableFile");
        UnitTests.TestReadWritableFile(local_file);
        Hello("TestAltUrl");
        UnitTests.TestAltUrl(local_file);
        Hello("TestCoordinateConversion");
        UnitTests.TestCoordinateConversion(local_file);
        Hello("TestCompareDumpMeta");
        UnitTests.TestCompareDumpMeta(local_file);
        Hello("TestVerifyMeta on local file");
        UnitTests.TestVerifyMeta(local_file);
        if (have_sd)
        {
          Hello("Test Upload to cloud");
          UnitTests.CopyFile(local_file, cloud_file, true);
          Hello("Test Download same file");
          UnitTests.CopyFile(local_file, cloud_file, false);
          Hello("TestIdToken");
          UnitTests.TestIdToken(local_file);
          Hello("TestCredentialsFrom");
          UnitTests.TestCredentialsFrom(cloud_file);
        }
        Hello("TestVerifyMeta after round trip");
        UnitTests.TestVerifyMeta(local_file);
        Hello("TestReadFile after round trip");
        UnitTests.TestReadFile(local_file);
        Hello("TestDeleteFile (local)");
        UnitTests.TestDeleteFile(local_file);
        if (have_sd)
        {
          Hello("TestDeleteFile (cloud)");
          UnitTests.TestDeleteFile(cloud_file);
        }
        Hello("TestExampleCallbacks");
        UnitTests.TestExampleCallbacks(3);
        Hello("All tests completed.");
      }
      catch (System.Exception ex)
      {
        System.Console.WriteLine("ERROR! {0}: {1}", ex.GetType().Name, ex.Message);
        try
        {
          System.IO.File.Delete(local_file);
          Hello("Cleaned up the local temp file.");
          // If cloud_file is orphaned it will hopefully be deleted
          // after a few days from one of the tests in the C++ code.
        }
        catch
        {
        }
        System.Environment.Exit(1);
      }
    }
  }
}
