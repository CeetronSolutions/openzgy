This is a very verbose documentation of:
- ZgyWriterArgsV2::guidsfrom() in OpenZGY
- checkCanKeepGuids() in both copy tools
- ZgyWriterArgs::checkSameMeta (never actually checked in)

These functions should ONLY be used if the target file has, or will
have, the exact same contents. One possible use case is when a ZGY
file is downloaded from or uploaded to the cloud.

These small differences are acceptable:
 - Bricks re-ordered on file.
 - Small changes to low resolution data, outside application's control.
 - Changes in compression or algorithm for low resolution data.

These differences are not acceptable:
 - Metadata such as size, geometry, units, annotation, etc.
 - Brick size, even though this is normally transparent.
 - Different storage type (int8/int16/float), even if widening.
 - Data range, if storage is int8 or int16.
 - Output to lossy compression, as that always changes data.
 - Input from lossy or lossless compression.
 - Anything else leaving target not bit-for-bit equal in lod0.

Enforcing these rules can be done inside Petrel, inside the copying
tools, and/or or inside OpenZGY. By "copying tools" I mean zgycopy,
zgytool, and any new code with similar functionality. The lower level
a check is done at, the better it protects against accidental use.
Checks in Petrel won't help if a user runs one of the copy tools from
the command line. Checks in the copying tools won't protect against
somebody writing a new copying tool or making an ill advised change to
the existing ones. Checks in OpenZGY would appear to be best. But most
of the rules are not possible to enforce in the OpenZGY API. Or it
would be quite complicated to do so.

When one of the copying tools requests to keep the guids of the input,
it is promising that the bulk data it will eventually send to the
writer will be exactly the same as the bulk data it reads from the
input file. This cannot be verified in the api. We need to trust the
copying tools to enforce the bulk related rules.

Inside the copying tools there should be checks that the following
options should not be used to change how the bulk is handled. Not all
these options are supported in both ZgyCopy and ZgyTool. Even fewer
options are currently used from Petrel. So if we had not been
concerned about abuse from the commend line tools, the list of options
to check would have been much shorter.

To summarize, at which levels should the rules be enforced?

- Protect against Petrel asking ZgyTool for --keepguids when it shouldn't.
- Also protect against users of command line tools with more options.
- Debatable: Protect against apps bypassing the tools and using OpenZGY.

Checking inside OpenZGY is arguably best. But only a small subset of
the rules are actually testable inside OpenZGY. And of those, some are
rather tricky to check at that level.

So, the copying tools can, and should, check their parameters as much
as possible. This protects from abuse by end users and by Petrel.
Petrel does not need to check anything if it only calls the ZgyTool
library and is ok with getting an exception when --keepguids is not
appropriate.

Here are the options broken down by where they might be checked.

- No OpenZGY check possible, bad things happen if not checked:
  - --noise, --brickcount, --input=(compressed file)
- No OpenZGY check but probably not critical:
  - --operation!=zgy2zgy, --input=(mocked), --output=(mocked).
- No OpenZGY check in the unlikely case where the size ends up unchanged:
  - --lod, --cropstart, --cropsize, --upsample, --mirror
- OpenZGY Would get a false positive on lossless output compression:
  - --sqnr*
- These can be checked both in the tools and in OpenZGY:
  - --osamplesize, --obricksize, --update*, --size
- Not yet implemented in tools, might be enforced in OpenZGY just in case:
  - --datarange, --geometry, --annotation, --units

**Do not try to enforce anything inside OpenZGY.**
That is my suggestion. Checking there adds little value because so few
options can be verified at that level. And it is tricky because the
input file is not available to OpenZGY when the output file is being
created. So the call to copyguids() would need to save the relevant
properties somewhere that can be checked in create. Or save a
reference to the reader. Yuck.

Footnote regarding --update: The current checks forbids that option
completely. Because it is only used for testing, and would IMHO be
pointless in combination with --keepguids. If I am wrong, it is
possible to open both the input and output for read and then check
those parameters that update would take from the existing output file
instead of copying from the input.

Footnote regarding --sqnr: The copying tools know whether they have
requested lossy or lossless compression. The OpenZGY API does not.
Because it only sees a compression callback. So, a check inside
OpenZGY would also need to reject writing a losslessly compressed
output file.

To check for input compression, the input file needs to be opened and
the brick lookup table scanned. Distinguishing between lossy and
lossless compression in the input file is not feasible. If that
becomes possible, the ban on lossless compression in the input file
might be removed.

Implementation note: If using metafrom() and not overriding those
settings later, that takes care of the first 4 bullets in the list
close to the start of this document. As of the current implementation.
Conversely, the current implementation of metafrom() can provide a
starting point for a function that checks those 4 bullets.

Note that any update to a file, even if just changing trivial
metadata, will cause the guids to be changed unless the file is opened
or reopened using guidfrom(). So, using guidfrom() on something that
isn't the last write access will not have any effect.

The API will only allow copying the guids from another file. It is not
possible to set arbitrary guids. It is also not supposed to be
possible to tell reopen() to keep its existing guids (i.e. they can
only be copied from somewhere else). A sufficiently devious programmer
might mock the input in order to set an arbitrary id. But the
copy-only design is to prevent accidental abuse of the option. There
would be nothing accidental about mocking the input file.

Implementation notes for checkCanKeepGuids in both copy tools:

Both tools use .metafrom() at some point. This means that the code
doesn't need to worry about what the default value type and range is
when not specified in the argument package. The code also doesn't need
to worry about anything (annotation, geometry, etc.) where zgycopy and
zgytool doesn't set these in any other way than copying. If we later
need to allow --update and --keepguids together, the tests get more
complex.

The test for input- and output compression, and changes to data type
and brick size, are quite relevant for Petrel. Most of the other
command line options will only be seen when invoking zgycopy or
zgytool via the command line.

As explained above, there are currently no checks for valid keepguid
in OpenZGY itself. I started writing such checks but only a few of the
possible checks have been implemented. Out of the few that are
possible there at all. In case checks turns out to be needed, I
include the unfinished code here.

```c++
/**
 * Check that the metadata is similar enough to allow keeping the old
 * guids when copying. Throw if they don't match. The tests are
 * necessary but not sufficient for keeping the guids.
 *
 * The checks compare floating point numbers for exact equal. This is
 * probably not an issue, as the numbers are probably just a plain
 * copy in the first place.
 *
 * The current implementation happens to check exactly the same fields
 * as metafrom() will set. If the application code uses args.metafrom()
 * and does not reset any of those parameters later, then the tests
 * will succeed. And the test here would hardly be needed.
 *
 * Typical usage will be to create a temporary ZgyWriterArgs and call
 * tempargs.metafrom(input). Then create the ZgyWriterArgs that will
 * actually be used to create the output. Call checkSameMeta()
 * immediately before creating the file.
 *
 * It is easier to do the check in spplication code or zgycopy/zgytool
 * because the input file is available there. Otherwise, guidsfrom()
 * might store a copy of the input metadata. Which is a bit messy.
 *
 * This doesn't allow for updating a file instead of creating it. In
 * that case, maybe use metafrom() on the output file as well? And, if
 * the real create args are manipulated later, make sure the bogus
 * writer-arguments get the same update. Calling merge() might help,
 * but be very careful with how that works.
 */
void
ZgyWriterArgs::checkSameMeta(const ZgyWriterArgs& a, const ZgyWriterArgs& b)
{
  const char* msg{nullptr};
  if      (a._size        != b._size)        msg = "size differs";
  else if (a._bricksize   != b._bricksize)   msg = "bricksize differs";
  else if (a._datatype    != b._datatype)    msg = "datatype differs";
  else if (a._datarange   != b._datarange)   msg = "datarange differs";
  else if (a._zunitdim    != b._zunitdim)    msg = "units differ";
  else if (a._hunitdim    != b._hunitdim)    msg = "units differ";
  else if (a._zunitname   != b._zunitname)   msg = "units differ";
  else if (a._hunitname   != b._hunitname)   msg = "units differ";
  else if (a._zunitfactor != b._zunitfactor) msg = "units differ";
  else if (a._hunitfactor != b._hunitfactor) msg = "units differ";
  else if (a._annotstart  != b._annotstart)  msg = "annotation differ";
  else if (a._annotinc    != b._annotinc)    msg = "annotation differ";
  else if (a._zstart      != b._zstart)      msg = "annotation differ";
  else if (a._zinc        != b._zinc)        msg = "annotation differ";
  else if (a._corners     != b._corners)     msg = "geometry differs";
  else msg = nullptr;
  if (msg != nullptr)
    throw OpenZGY::Errors::ZgyUserError("Cannot copy guids because " + std::string(msg));
}
```
