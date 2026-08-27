
## C# API

This is the C# API for OpenZGY.

The design goal for the C# API is to make it as similar as possible to
the C++ API and as complete as possible.

Names have the same capitalization as in C++ even where that violates C#
style guidelines. Some functions may be changed to properties where
this seems intuitive. Function signatures and argument types may only
be changed if there is a very good reason to.

The dump() debug function in several classes has been renamed toString()
and will return a multi-line string instead of writing to an ostream.

Note that callback functions implemented by application code
(logger, progress reporter, and token refresh) must be thread safe.

The C# API is implemented as a wrapper around the C API using P/Invoke,
which in turn will is a wrapper around the C++ API. No C++/CLI needed.
This makes the code simpler to link into applications and also makes
is compatible with Mono. The trade off is that the wrapper project
becomes fairly large and there is less checking of function signatures
and types.

Function documentation has in most cases not been copied from C++.
Developers need to refer to the C++ documentation instead.
Copying it would clearly be a benefit when using intellisense.
But this needs to be weighed against the near certainty of the
text becoming out of date after a while.

For more details see [README.md](../native/src/capi/README.md)
in the C API implementation.

### Methods that do not need to be wrapped

| Method | Reason |
| ------ | ------ |
| FileStatistics::segmentSizes() | Unlikely to be used |
| IZgyMeta::brickcount() | Unlikely to be used |
| IZgyTools::transform() | Convenience, rather inconvenient with interop. |
| IZgyTools::transform1() | Convenience, rather inconvenient with interop. |
| ZgyWriterArgs::compressor(compressor_t) | Would be too expensive to use. |
| ZgyWriterArgs::lodcompressor(compressor_t);
| ZgyWriterArgs::compressor(name, args[]); | More convenient to use zfp_compressor |
| ZgyWriterArgs::lodcompressor(name, args); | More convenient to use zfp_lodcompressor |
| SeismicStoreIOContext::sdtoken(value, type); | Deprecated overload |
| SeismicStoreIOContext::sdtokencb(value, type); | Deprecated overload |
| SeismicStoreIOContext::debug_trace(); | For debugging, pointless here. |
| SeismicStoreIOContext::segsizebytes() | For unit test only. |
