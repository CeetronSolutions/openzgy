# ZGY APIs other than C++

## C# API

The design goal for the C# API is to make it as similar as possible to
the C++ API and as complete as possible.

Names have the same capitalization as in C++ even where that violates
C# style guidelines. Some functions may be changed to properties where
this seems intuitive. Function signatures and argument types may only
be changed if there is a very good reason to.

The dump() debug function in several classes has been renamed toString()
and will return a multi-line string instead of writing to an ostream.

The C# API is implemented as a wrapper around the C API using P/Invoke,
which in turn will is a wrapper around the C++ API. No C++/CLI needed.
This makes the code simpler to link into applications and also makes
is compatible with Mono. The trade off is that the wrapper project
becomes fairly large and there is less checking of function signatures
and types.

Implementation using SWIG was briefly considered but rejected.
Admittedly not much research was done before dropping that idea.

## C API

This is designed as a chatty interface with some 150 extern "C"
functions. One for each member function in C++ plus a few more.

The C API is expected to be a tool to write wrappers in higher
level languages. Not to be called directly from applications.

## Python API

As with the C# API it is possible to write a Python API on top of the
new C API. This could be implemented entirely in Python using the
ctypes module. It is not clear whether this would be a good idea.
Given that there already is a Python wrapper implemented as a C++
Python extension. A new Python wrapper would definitely not be high
priority.

## Second C++ API

It is possible but probably pointless to create a new C++ API
as a wrapper around the new C API. The only use for this would
be to run the entire regular test suite with that second C++
API and thereby get a really good test of the new C (but not C#)
API. Yes this still sounds like mostly useless work because
the wrapper itself shouldn't require that much testing.
Just forget that I mentioned it ;-)

## C implementation

In general, one extern "C" function is written for each free
function and each member function. Including overloads. The estimate as of
the start of this project is to write about 150 functions in each of
the C and the C++ projects. Most of the functions will be practically
one-liners.

### Version control

There are two approaches being used. The first one is that any user of
the API is expected to call oz_checkVersion(OZ_C_API_CURRENT_VERSION)
before any other function. If this returns an exception handle then
the entire API is unavailable due to a version mismatch.

The second approach is that any existing oz_XXX function might be
deprecated and possibly replaced with oz_XXX_v2, oz_XXX_v3, et cetera.
The deprecated function must be kept but is allowed to return error
when called. This kind of change is backwards compatible.

### Calling conventions

Nearly all functions expect an opaque handle as the first argument
and returns another opaque handle.

The client is always required to free the returned handle when done
with it. Some handles are long lived (e.g. an open file handle) while
others are ephemeral and only used to marshal complicated function
results.

Internally there are some 12 different types of handle. The client is
supposed to keep track of what each handle is for. So there is no
specific function to get the type of a handle. With the exception of
the exception handle, see below.

Passing a wrong handle type to a function will return an error (an
exception handle) which the client can recognize. So this is not a big
deal. Passing a freed handle or other garbage triggers undefined
behavior. The wrapper tries to detect this and abort the application
but that cannot be guaranteed.

For debugging it might be possible to reliably detect use of a freed
handle (by keeping them in an internal pool) and a leaked ephemeral
handle (if not already freed when the associated long lived handle is
freed). This feature is not implemented and might not ever be.

All methods that should return a particular handle can return an
exception handle instead. This handle contains information about an
exception caught by the wrapper and converted to a return status. This
is the one case where the client is allowed to check the type of a
handle. Calling oz_resultIsSuccess(handle) is equivalent to checking
handle->type != EXCEPTION and can be invoked on any valid handle.
As a reminder, oz_resultIsSuccess() will, as any other method, crash
if the handle is invalid.

Functions that return a handle as the function result (almost all of
them) will return a scalar result via a single \[Out] argument. Short
fixed size arrays such as float\[2] or long long\[3] will be returned
via multiple \[Out] arguments. A string result will be returned as a
string handle in the function result, from which the actual string may
be retrieved. Fixed size bulk data for read and write buffers and for
the histogram are allocated and pinned by the client and passed to the
library which will copy-in and copy-out as needed.

It is possible for the wrapper to allocate data itself and pass it as
a pointer to the client. If so, that data will be released when the
handle is closed by the client. The client should not free any buffers
from an \[Out] parameter even when it knows that the library allocated
them. Most likely the implementation won't need to use this feature.

The STATISTICS, HISTOGRAM, and FILESTATS handles represent POD structs.
Initially these also use a chatty interface were members are retrieved
one at a time. If this proves to be too expensive then an alternative
access method might be provided that copies out the entire struct in
one go. This would probably be more hassle to maintain.

### Ephemeral handle types

| Type    | Used for |
| ------- | -------- |
| SUCCESS | Represents a successful result. Holds no other information. Might actually use a nullptr instead, to be decided later. |
| CLEANUP | From the clients point of view, the same as SUCCESS. Internally it allows freeing any bulk data allocated by the library. |
| ERROR   | Represents a caught exception. May be returned when another type is expected. Details about the exception can be retrieved. |
| STRING  | Wraps a returned string. Only needed a couple of places. The string itself can be retrieved from the handle, in a way that simplifies memory management. |

In C#, wrapping functions that only returns a pass/fail status i.e.
SUCCESS, CLEANUP, or ERROR will be almost exclusively boiler plate code.
With a suitable helper function they might even be one-liners.

```c++
TYPE wrap_XXX(ZgyHandle handle, ...)
{
  TYPE result;
  ZgyHandle ret = oz_XXX(handle, [Out]result, ...);
  if (!oz_isOk(ret)) {
    // marshal and copy the error message
    string msg = oz_getErrorString(ret);
    oz_freeHandle(ret);
    // Depending on how error handling works in this language
    throw msg;
  }
  else {
    oz_freeHandle(ret);
    return result;
  }
}
```

### Long lived handle types

| Handle type | Represents            | # of members   | C linkage prefix |
| ----------- | --------------------- | -------------- | ---------------- |
|READER+WRITER| IZgyMeta,IZgyTools    | 24 getters     | oz\_meta\_       |
|READER       | IZgyReader            | 7              | oz\_reader\_     |
|WRITER       | IZgyWriter            | 16             | oz\_writer\_     |
|UTILS        | IZgyUtils             | 4              | oz\_utils\_      |
|STATISTICS   | SampleStatistics      | 5 getters      | oz\_statistics\_ |
|HISTOGRAM    | SampleHistogram       | 4 getters      | oz\_histogram\_  |
|FILESTATS    | FileStatistics        | 20 getters     | oz\_filestats\_  |
|IOCONTEXT    | SeismicStoreIOContext | 24 setters     | oz\_ssiocontext\_  |
|WRITERARGS   | ZgyWriterArgs         | 23 setters     | oz\_writerargs\_ |

### Return handle types

Almost all C functions take a handle as input and returns another handle.
The following table shows the handle types returned. The expected type of
the input handle can be derived from the prefix as shown in the table above.

| C++ member name      | Returns    | Input handle  |
| -------------------- | ---------- | ------------- |
| (any not listed)     | SUCCESS    | (see prefix)  |
| ZgyMeta::brickcount  | (error)    | READER,WRITER |
| ZgyMeta::zunitname   | STRING     | READER,WRITER |
| ZgyMeta::hunitname   | STRING     | READER,WRITER |
| ZgyMeta::verid       | STRING     | READER,WRITER |
| ZgyMeta::statistics  | STATISTICS | READER        |
| ZgyMeta::histogram   | HISTOGRAM  | READER        |
| ZgyMeta::filestats   | FILESTATS  | READER,WRITER |
| ZgyReader::open      | READER     | (none)        |
| ZgyWriter::create    | WRITER     | (none)        |
| ZgyWriter::clone     | WRITER     | (none)        |
| ZgyWriter::reopen    | WRITER     | (none)        |
| ZgyUtils::utils      | UTILS      | (none)        |
| ZgyUtils::alturl     | STRING     | UTILS         |
| IOContext::new       | IOCONTEXT  | (none)        |
| IOContext::dbg_trace | (error)    | IOCONTEXT     |

### Main functions not wrapping a specific C++ method and often not returning a handle

```c++
ZgyHandle oz_checkLibraryVersion();
bool      oz_resultIsSuccess(ZgyHandle);
char*     oz_resultGetExceptionMessage(ZgyHandle);
char*     oz_resultGetExceptionType(ZgyHandle);
char*     oz_resultGetString(ZgyHandle);
void      oz_freeHandle(ZgyHandle);
```

### Estimated size of the C and C# API

See the tables above.

22 main methods, 53 getters, 47 setters, and maybe 8 special purpose
functions that are the only ones not returning a handle. Grand total
**130**.

In addition there are **4** enums that will be treated as integers in
the C API and cast to a real enum in C#, with the enum tags manually
kept in sync.

```c++
enum SampleDataType;
enum UnitDimension;
enum DecimationType;
enum FinalizeAction.;
```

There are **12** different ZGY exceptions that can be thrown as well
as some SDAPI exceptions that can be passed thru. The simplest
approach is to ignore the specific exception type and report just the
message. It would not be too much extra work to distinguish between
the 12 known types and perhaps a couple of the SDAPI ones.

There are **3** callback function types that need to be handled.

```c++
typedef std::function<std::string()> tokencb_t;
typedef std::function<bool(int, const std::string&)> logger_t;
typedef std::function<bool(std::int64_t,std::int64_t)> progress_t;
```

So, including about 20 things that are not functions the checklist for
what needs to be written contains 150 bullet points for the C wrapper.
And presumably about as many in the C# wrapper. So, total of 300 TODO
items. Fortunately almost all of them trivial.

## Performance concerns

The C API is a chatty interface and assumes that the cost of interop
will be negligible. The C# API and other potential language bindings
are not expected to cache any data for the same reason. If this
assumption turns out to be wrong then consider the following.

ZgyMeta getters, SampleStatistics, SampleHistogram, and FileStatistics
can be cached on the C# side. At least for files opened for read. For
writable files there will need to be logic to invalidate the cache
after writes. There may be increased maintenance cost.

The above classes in the C API can be changed to pass entire POD
structs instead of making one call per field. There is a significant
risk of getting the structs out of sync. And more issues with
compatibility between versions. Also, while C# ought to be able to
marshal these complex types, other languages might not.

## Checklist after a new function has been added to the C++ API.

Implement an extern "C" function in native/src/capi/capi.cpp named
oz_CLASS_NAME. Scalar results are returned via ref or \[Out]
parameters. Enums and booleans are treated as int in the C API.

Declare the function as external in managed/OpenZGY.Managed/CLASS.cs.
Be careful to make the function signature match. The compilers will
*not* help with this.

Implement a C# CLASS.NAME() function or CLASS.NAME getter in
managed/OpenZGY.Managed/CLASS.cs.

#### Functions returning scalars

The simplest are functions returning void or scalars, including short
fixed size arrays where the return is done via multiple \[Out]
parameters.

```c++
// In the C++ project with C linkage
return protect([&]() -> void {ZgyCLASSHandle::get(handle)->NAME());});

// In the C# class
[DllImport("OpenZGY.dll")] private static extern ZgyHandle
oz_CLASS_NAME(ZgyHandle handle, ref T return_value, ..., T arg_value, ...)

// In the C# class inside the getter or function CLASS.NAME
TYPE ret = default; // Not if NAME() in C++ returns void.
Tools.freeHandle(oz_CLASS_NAME(Handle, ref ret, ...));
return ret;
```

The C++ lambda always returns void because the actual function return
from the C API is a handle that only reports pass (no details) or fail
(with message). Returning any real results via \[Out] parameters adds
a little more code.

#### Functions reading or writing bulk data

As above, with a pointer to the buffer owned by the C# code and passed
as \[In] or \[Out]. NAME() in C++ typically returns void so the C#
wrapper also does that.

#### Functions returning strings

```c++
// In the C++ project with C linkage
return protectReturn([&]() {
  return new ZgyStringHandle
    (ZgyCLASSHandle::get(handle)->NAME());});

// In the C# class
[DllImport("OpenZGY.dll")] private static extern ZgyHandle
oz_CLASS_NAME(ZgyHandle handle, T arg_value, ...)

// In the C# class inside the getter or function CLASS.NAME
return Tools.resultGetString(Tools.checkHandle(oz_CLASS_NAME(Handle)));
```

#### Functions returning a class instance

```c++
// In the C++ project with C linkage
return protectReturn([&]() {
  return new ZgyRESULTCLASSHandle
    (ZgyCLASSHandle::get(handle)->NAME());});

// In the C# class
[DllImport("OpenZGY.dll")] private static extern ZgyHandle
oz_CLASS_NAME(ZgyHandle handle, T arg_value, ...)

// In the C# class inside the getter or function CLASS.NAME
return new CLASS(Tools.checkHandle(oz_CLASS_NAME(Handle)));
```

#### Functions expecting a callback parameter

Definitely not trivial, especially if the callback itself returns
variable length data such as a string. See the source code for
examples with the logger callback, progress callback, and the token
refresh callback.

