# Universally Unique Identifiers

uuids and guids (the names are usually used interchangeably) are
essentially just 128 random bits but are notoriously difficult to deal
with. I need a function to generate a new uuid that can be stored on a
ZGY file as 16 consecutive bytes. I also need a function to convert
those 16 consecutive bytes to a human readable string. That can't be
too difficult, right?

See [RFC 4122](https://tools.ietf.org/html/rfc4122)
and [Wikipedia](https://en.wikipedia.org/wiki/Universally_unique_identifier).

## Non standard representation in ZGY files

The uuids in a ZGY file on disk are stored not big-endian as
[RFC 4122](https://tools.ietf.org/html/rfc4122)
requires. Instead they are stored piecewise little endian.
This affects how the raw bytes are converted to a canonical string.

In more detail: uuids can be stored as:
1. A canonical 36-character string xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx.
2. An array uint8[16] corresponding 1:1 to the canonical string.
3. A struct {uint32,uint16[2],uint8[8]} of size 16 bytes.

Since the representation in (3) contains multi-byte integers the
layout will differ between little-endian and big-endian machines.
[RFC 4122](https://tools.ietf.org/html/rfc4122)
says the canonical representation is big-endian so (2) and (3) will be
byte-by-byte identical only on a big-endian architecture. ZGY stores
uuids as little-endian instead.

This cannot easily be changed because that would break code that reads
an uuid from a ZGY file and converts it to a string. Existing ZGY
files would then appear to have suddenly changed their uuids.

Note that code that tries to "convert" between (2) and (3) by casting
is wrong in any case. It will only work on one (little-endian or
big-endian) architecture.

To avoid confusion in ZGY an access library should prevent application
code from seeing the raw bytes. Only the canonical string
representation should be accessible. This isn't quite the case for the
old ZGY-Public accessor but the issue is manageable,

* ZGY-Public hides the uuids completely from the application code

* OpenZGY/C++ will hide the raw bytes in the API and only expose the
  stringized version that is computed in native/src/impl/guid.cpp.

* OpenZGY/Python will try to do the same, or possibly just hide the
  uuids completely. In Python the data hiding is implemented just as a
  naming convention so it isn't quite as well protected as in
  OpenZGY/C++.

* ZGY-Internal, used only by Petrel via Salmon, exposes raw bytes but
  they are quickly converted to strings by
  Salmon/Shared/PublicTypes/UniqueId.cpp.
  To my knowledge this is the only place in Salmon that converts.
  This one is more difficult to verify since the raw bytes are in fact
  allowed to escape from the ZGY-Internal API. Note that the code in
  UniqueId.cpp is technically incorrect and will give wrong results if
  compiled on a big-endian machine. This is not a problem for Petrel.
  But please don't copy/paste that code anywhere.

# API for dealing with uuids

On Linux there are two libraries, libuuid.so and libossp-uuid.so, that
can convert from raw 16 bytes to a string and on Windows there is at
least one. These algorithms **DO NOT MATCH**. This is where the
problems start. I can roll my own library that try to mimic either of
those algorithms but there is a risk that I end up implementing yet
another incompatible version.

I believe cause of the problem is this: 

An uuid stored in memory is traditionally laid out as
uint32,unit16[2],uint8[8]. On disk it just looks like 16 consecutive
bytes. Files on disk clearly need to be interoperable between
little-endian and big-endian machines. It is less obvious that the
long and the two shorts need to be accessible as such, but this has to
do with uuids based on a timestamp instead of random bits. The
standard says to use big-endian for interop. So, on the x86
architecture there is some non-obvious byte swapping going on to
convert between a uuid struct and a 16-byte array that can be written
to disk.

If the APIs you are working on expect uuids to be passed as uint8[16]
buffers that are already in big endian order then no byte swapping is
needed.

If you are on a little endian machine and are using an api that
expects or returns a struct uuid and you simply cast that to an
uint8[16] then you are doing something very wrong. This is what Petrel
(or actually Salmon) is doing today in UniqueId.cpp.

So the root cause is that the api on windows expect a struct gui while
those on linux expect an uint8[16].

To complicate the issue, the
[wikipedia page](https://en.wikipedia.org/wiki/Universally_unique_identifier)
very clearly explains that the older type of guids created by
microsoft, tagged as variant 2, have the opposite byte ordering.
Testing suggests that this is incorrect.
[RFC 4122](https://tools.ietf.org/html/rfc4122)
is less clear on that subject. I choose to trust my testing and assume
wikipedia is wrong or at least what they say doesn't apply to that
particular API that I am using.

## UUID types

The simple way of generating a new uuid is to use a good random number
generator to create ~128 bits of random data. The assumption being
that the random number generator being used has enough entropy to
virtually guarantee there will never be a collision. In spite of the
birthday paradox.
Unfortunately, sufficient entropy **CANNOT BE GUARANTEED**. The
quality depends on the implementation of the random generator used.
Which might be OS dependent. It can even depend on hardware. The
alternative is to use a much more complex algorithm using time,
ethernet MAC address, and other information. This has its own set of
problems. Such as identifying the machine on which the uuid was
created (privacy issue) and what to do on a machine where there is no
ethernet address or (for virtual machines) where the ethernet MAC
address is not unique.

## OpenZGY choices

Both the generate and the format algorithms should be implemented
locally in the OpenZGY code. Because interoperability between Linux
and Windows is the most important criteria.

## Generate

When generating a new uuid I trust that the random number generator my
code uses has sufficient entropy. It helps that the user typically
doesn't deal with very many ZGY files. They might be numbered in the
thousands but hardly in the millions. RISK: Bad number generators and
the birthday paradox. MITIGATION: Be as specific as possible when
choosing the random generator to use. MITIGATION: Share the random
number generator between threads.

I will generate uuids tagged as variant 1 (DCE) and version 4
(random). This also appears to be the default on windows when calling
UuidCreate().

This also mitigates the uncertainty about byteswapping of variant 2 uuids.

## Format

For backwards compatibility try to make the storage to string
conversion match what is done in Salmon today. This means treating the
uuid stored on file as piecewise little-endian in violation of RFC
4122. The reason is that only the ZGY-Internal API exposes the uuids
from the ZGY file, and only Petrel uses ZGY-Internal. So only Petrel /
Salmon will have done an explicit raw to string conversion of uuids.
And Petrel, last time I checked, only runs on Windows. ZGY-Public
which is built on top if ZGY-Internal hides the uuids on read. RISK: I
might not manage to replicate the Windows algorithm correctly. RISK:
There might be a mismatch between OpenZGY/Python and OpenZGY/C++.

OpenZGY should only expose string uuids in its api. This makes the
conversion from bytes stored on file to string an implementation
detail. Just like the choice of using mostly little-endian storage for
integers. RISK: This in itself won't help with files created by the
old ZGY-Public or ZGY-Internal libraries. Hence the previous bullet.

Mitigation: Identify the place or places in Salmon where the
conversion is done. It looks like this is just one place: class
UniqueId in Salmon/Shared/PublicTypes/UniqueId.cpp is a candidate.

## Deja vu

It is not surprising that this issue has shown up before. See the mail
thread with subject gv34oban_218install.bin from 2009. The solution
being discussed was to fix it on the seismic server side. RISK: It is
possible that a kludge was put into Salmon or Petrel instead. The last
entry in that mail thread indicates that the plan was to fix it inside
the seismic server. But, that mail was sent while the Trondheim office
was being shut down and several of the developers were leaving the
company.

Mitigation: If we are confident that
Salmon/Shared/PublicTypes/UniqueId.cpp is the only place where
conversion is done then whatever happened with that issue is
irrelevant.

## TL;DR

All of these steps taken together should end up with a very minor residual risk.

## APIs and packaging on Linux

The following is only relevant if we do not use a private
implementation for guids. So this is about understanding the old ZGY
accessor. It is probably a bad idea to change anything there. Except
possibly moving away from libossp-uuid in CentOS or at least linking
it statically.

### APIs on Linux

On Linux systems there are two different uuid libraries to choose
from. Note that the naming convention is really really bad making it
easy to confuse the two.

The situation on Debian and Ubuntu is as follows:
- libuuid1, uuid-dev, uuid-runtime:
  - installs libuuid.so, /usr/bin/uuidgen, <uuid/uuid.h>, etc.
- libossp-uuid16, libossp-uuid-dev, uuid
  - installs libossp-uuid.so.16.*, /usr/bin/uuid, <ossp/uuid.h>, etc.
- Note that the "uuid" command line tool is actually an "ossp-uid".
- Note that the "uuid-dev" is NOT the dev package (headers) for "uuid".
- Note that in CentOS the OSSP headers go in <uuid.h> not <ossp/uuid.h>

The situation on CentOS and Fedora is as follows:

- libuuid, libuuid-devel, util-linux:
  - installs libuuid.so, /usr/bin/uuidgen, <uuid/uuid.h>, etc.
- uuid, uuid-devel
  - installs libossp-uuid.so.16.*, /usr/bin/uuid, <uuid.h>, etc.
  - On CentOS 8 you need to enable powertools to get uuid-devel.

It looks like libuuid1 is always installed by default on all distros.

In the Salmon build, which uuid library is used depends partly on which
packages have been installed in the build environment and partly on how
the tests in two Makefiles work. Thanks to Docker the behavior is now
deterministic but it takes some time to understand what is going on.

On Debian and Ubuntu: Salmon and ZGY chooses to use libuuid1.
Installing the OSSP version in the build environment will not change
this because the package ends up in a different place than where the
Makefile looks for it.

On CentOS and Fedora: Salmon and ZGY chooses to use the ossp version.
Removing that package will automatically switch to use libuuid1.
Nowadays when all official builds use Docker the ossp libraries will
always be used for CentOS and Fedora. Prior to that it was a bit
arbitrary which library was used since it depended on the build
machine. Ouch.

All platforms link dynamically with the chosen library, libuuid.so.1
or libossp-uuid.so.16

With 20/20 hindsight Salmon ought to have used libuuid1
unconditionally. The reason it doesn't is probably the confusing
naming convention. Making it less than obvious that libuuid1 would
always be available at runtime. And the only tricky bit is to figure
out how to install the needed headers.

The situation today is that users of ZGY-Public on CentOS / RedHat
will need to install the "uuid" package. If we really need to fix this
then either switch to using libuuid1 for all distros or switch to
linking libossp-uuid statically. Probably not worth the trouble.

### APIs on Windows.

There is a method UuidToStringA() in <rpc.h> and rpcrt4.lib. Beware
that the Windows API expects struct UUID as input while the Linux apis
deal with 16-byte arrays. There is also class UUID in .NET but that is
not relevant for OpenZGY.
