# Migration guide from ZGY-Public / ZGY-Cloud to OpenZGY

## ZGY-Public (C++)

The API is completely different. Since the core API (i.e. without the
cloud extensions) is very small I don't think this will be a problem.

## ZGY-Public (Python)

The new API is very similar to the old one. Some stringly typed enums
(zunitdim, hunitdim, datatype, sourcetype) have been changed to become
real Python enums.

## ZGY-Cloud (C++ and Python)

The old ZGY accessor had a separate ZGY-Cloud module with its own
"back door" api for cloud specific operations such as setting the
credentials. The ZGY-Cloud module (both the C++ and Python API)
have not been published but are used in Petrel and in the ML project.

The following list shows the replacements. Note that the
basic create/read/write api is unchanged.

Unless otherwise noted, zgy_someFunction() in the old ZGY-Cloud module
corresponds to zgycloud.someFunction() in the old Python wrapper.

```
configure(backend, name, sdkey, sdurl)
setToken(backend, token, type, global)
setTokenCb(backend, callback, type) # C++ only
setSegmentSize(size_in_mb)
  => Information to be specified in the iocontext.

deleteFile(name, token)
  => pure python: openzgy.api.ZgyUtils.delete()
     or sdglue.SdUtil(*cred).delete(name)
  => wrapper: openzgycpp.ZgyUtils.delete
  => native: ZgyUtils.deleteFile (TODO-Low test?)

enableRealCache(on, size_in_mb)
setCacheHooks(hooks) # C++ only
getCacheHooks() # C++ only
  => Not supported by OpenZGY.

getTestToken(backend)
  => pass "FILE:carbon.slbapp.com" as the IOContext token.

lastError()
resetError()
  => N/A because all errors are reported synchronously.

listObjects(name, recursive, token)
listBuckets(project) # C++ only
listProject() # C++ only
  => Not supported, and never was officially.

setLogger(callback, level)
getCacheStats(...) # C++ only
  => Debugging and logging is somewhat different.

getLegalTag(name, token) # C++: zgy_{get,set}Tag
setLegalTag(name, token, legaltag)
getMetaData(name, token)
setMetaData(name, token, metadata)
getSeismicMeta(name, token) # C++: zgy_{get,set}SeisMeta
setSeismicMeta(name, token, seismicmeta)
  => TODO-Medium investigate further.
     legaltag, seismicmeta, writeid settable in iocontext
     but cannot(?) be modified later.
     metadata might be missing

getContext(backend, name, global)
setContext(backend, name, value, global)
  => TODO-Medium investigate further.
  => This is a general back door mechanism, used for ...
  => and replaced by ...

zgy_checkAccess(name, token) # C++ only
  => TODO-Low might need to support this.
```
