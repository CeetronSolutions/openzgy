// Copyright 2017-2021, Schlumberger
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

#ifdef _WIN32
#pragma warning(disable: 4251) // Needs dll-interface
#endif

#define PY_SSIZE_T_CLEAN
#undef NDEBUG

// Major kludge ahead. Azure DevOps doesn't provide debug versions of
// python3.dll. There are several equally bad choices.
//
//   (1) Don't provide debug builds of the wrapper at all.
//   (2) Use build_ext --debug and the trick below.
//       Wrapper and OpenZGY are debug, Python is release.
//   (3) Turn off build_ext --debug. The wrapper is compiled in release
//       mode, linking agaist release Python and debug OpenZGY.
//       Probably the riskiest approach since the wrapper makes no
//       attempt to be safe with mixing debug and release.
//       Does python3.dll do that? Wishful thinking.
//   (4) Build on a self hosted server or in docker and install the
//       complete Python including the debug library.

#if defined(_WIN32)
#pragma push_macro("_DEBUG")
#undef _DEBUG
#endif
#include <Python.h>
#include <structmember.h>
#if defined(_WIN32)
#pragma pop_macro("_DEBUG")
#endif

#include <api.h>
#include <writerargs.h>
#include <iocontext.h>
#include <exception.h>
//#include <impl/logger.h> // inlined below, BAD SMELL!

#include <atomic>
#include <memory>
#include <omp.h>
#include <iostream>
#include <assert.h>
#include <functional>
#include <limits>

#ifdef _WIN32
#include <process.h> // getpid()
#endif

/**
 * KLUDGE! Duplicated from ../native/src/impl/logger.h.
 * When installing an sdist using a pre-installed OpenZGY
 * this internal header might not be available.
 */
namespace InternalZGY {
class OPENZGY_API LoggerBase
{
public:
  typedef std::function<bool(int, const std::string&)> LoggerFn;
  static int getVerboseFromEnv(const char *envname);
  static bool logger(const LoggerFn& logger, int priority, const std::string& str = std::string());
  static bool logger(const LoggerFn& logger, int priority, const std::ios& ss);
  static LoggerFn emptyCallback();
  static LoggerFn standardCallback(int level, const std::string& prefix, const std::string& suffix);
};
}

/*
 * MODULE_NAME vs. FULL_MODULE_NAME:
 *
 * The Pickle module will save __module__ + '.' + __name__ of types, and
 * expects __module__ to be a full name i.e. "openzgycpp.wrapper" and not
 * just "wrapper". This module name is set several places:
 *
 *    static PyTypeObject {tp_name=...} // E.g. for ZgyReaderClassType.
 *    module arg to procedural enum constructor
 *    PyErr_NewException(name, ...)
 *
 * As far as I can see, __name__ should always be unqualified. When
 * statically constructing a PyTypeObject with tp_name "a.b.c" the
 * runtime assigns "a.b" to __module__ and just "c" to __name__. So
 * I assume the same applies to other places that set __name__.
 *
 * This means that PyModuleDef.m_name, set whan statically constructing
 * the module itself, should use "wrapper" and not "openzgycpp.wrapper".
 * The fully qualified module name will be __package__ + '.' +__name__.
 * (Why not __module here?)
 */
#define MODULE_NAME "wrapper"
#define FULL_MODULE_NAME "openzgycpp." MODULE_NAME

using namespace OpenZGY;
using namespace OpenZGY::Errors;
using InternalZGY::LoggerBase;

// The code that tries to MT-fetch one brick at a time has been removed.
// TODO-Low: converters should have a round=False argument.
// TODO-Low: convert arrays instead of just single values.

static PyObject* _zgy_error = nullptr;
static PyObject* _zgy_error_FormatError = nullptr;
static PyObject* _zgy_error_NeedOldLibrary = nullptr;
static PyObject* _zgy_error_UpdateRules = nullptr;
static PyObject* _zgy_error_CorruptedFile = nullptr;
static PyObject* _zgy_error_UserError = nullptr;
static PyObject* _zgy_error_InternalError = nullptr;
static PyObject* _zgy_error_EndOfFile = nullptr;
static PyObject* _zgy_error_SegmentIsClosed = nullptr;
static PyObject* _zgy_error_Aborted = nullptr;
static PyObject* _zgy_error_MissingFeature = nullptr;
static PyObject* _zgy_error_IoError = nullptr;
static PyObject* _enum_SampleDataType = nullptr;
static PyObject* _enum_UnitDimension = nullptr;
static PyObject* _enum_DecimationType = nullptr;
static PyObject* _enum_FinalizeAction = nullptr;
static PyObject* _enum_LodMode = nullptr;
static PyObject* _namedtuple_Statistics = nullptr;
static PyObject* _namedtuple_Histogram = nullptr;
static PyObject* _namedtuple_FileStats = nullptr;

static PyObject*
_raise_simple_error(PyObject* type, const char *message)
{
  PyErr_SetString(type==NULL ? _zgy_error : type,
                  message ? message : "error");
  return NULL;
}

static PyObject*
_raise_simple_error(PyObject* type, const std::string& message)
{
  return _raise_simple_error(type, message.c_str());
}

static PyObject*
_raise_ex()
{
  try {
    throw;
  }
  catch (const ZgyNeedOldLibrary& ex)  { return _raise_simple_error(_zgy_error_NeedOldLibrary,  ex.what()); }
  catch (const ZgyUpdateRules& ex)     { return _raise_simple_error(_zgy_error_UpdateRules,     ex.what()); }
  catch (const ZgyFormatError& ex)     { return _raise_simple_error(_zgy_error_FormatError,     ex.what()); }
  catch (const ZgyCorruptedFile& ex)   { return _raise_simple_error(_zgy_error_CorruptedFile,   ex.what()); }
  catch (const ZgyUserError& ex)       { return _raise_simple_error(_zgy_error_UserError,       ex.what()); }
  catch (const ZgyInternalError& ex)   { return _raise_simple_error(_zgy_error_InternalError,   ex.what()); }
  catch (const ZgyEndOfFile& ex)       { return _raise_simple_error(_zgy_error_EndOfFile,       ex.what()); }
  catch (const ZgySegmentIsClosed& ex) { return _raise_simple_error(_zgy_error_SegmentIsClosed, ex.what()); }
  catch (const ZgyAborted& ex)         { return _raise_simple_error(_zgy_error_Aborted,         ex.what()); }
  catch (const ZgyMissingFeature& ex)  { return _raise_simple_error(_zgy_error_MissingFeature,  ex.what()); }
  catch (const ZgyIoError& ex)         { return _raise_simple_error(_zgy_error_IoError,         ex.what()); }
  catch (const ZgyError& ex)           { return _raise_simple_error(_zgy_error,                 ex.what()); }
  catch (const std::exception& ex)     { return _raise_simple_error(_zgy_error,                 ex.what()); }
  catch (...)                          { return _raise_simple_error(_zgy_error, "(unknown exception"); }
}

static std::string
_toString(PyObject* obj)
{
  if (!obj)
    return "(null)";
  else if (obj == Py_None)
    return "None";
  PyObject* strobj = PyObject_Str(obj);
  if (!strobj) {
    PyErr_Clear();
    return "(error)";
  }
  const char *cstring = PyUnicode_AsUTF8(strobj);
  if (!cstring) {
    PyErr_Clear();
    return "(error)";
  }
  Py_XDECREF(strobj);
  return std::string(cstring);
}

/**
 * Return (long)obj.value.
 * On error set an exception and return -1.
 * A "value" < 0 is also flagged as an error.
 */
static long
_decodeEnum(PyObject *obj, PyObject* typeobj)
{
  if (!PyObject_TypeCheck(obj, (PyTypeObject*)typeobj)) {
    _raise_simple_error(PyExc_TypeError, "Expected " + _toString(typeobj));
    return -1;
  }
  PyObject *enum_value = PyObject_GetAttrString(obj, "value");
  if (!enum_value)
    return -1;
  long result = PyLong_AsLong(enum_value);
  Py_XDECREF(enum_value);
  if (result < 0) {
    _raise_simple_error(_zgy_error_InternalError, "Cannot get enum tag");
    return -1;
  }
  return result;
}

/////////////////////////////////////////////////////////////////////////////
///   ZgyReader class   /////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/**
 * Using the pimpl pattern, among other things because some of the members
 * are complex C++ types which don't play well with tp_init, tp_dealloc
 */
struct ZgyClassImpl
{
  typedef std::function<bool(int, const std::string&)> LoggerFn;
  std::string filename;
  int numthreads;

  // set to true if we want to unlock the python GIL every time we read / write.
  bool unlockgil;

  // The logger verbosity is captured in the lambda assigned to logger_.
  // The logger_verbose_ member is only there so the current "verbose"
  // level can be read from application code.
  int logger_verbose_;
  LoggerFn logger_;

  //MyMetaData meta;
  //ZgyApi::ReaderPtr read_accessor;
  //ZgyApi::WriterPtr write_accessor;
  std::shared_ptr<IZgyReader> reader_;
  std::shared_ptr<IZgyWriter> writer_;
  std::shared_ptr<IZgyTools> meta_;
  std::shared_ptr<IZgyUtils> utils_;

  PyObject *debug_trace_;

  ZgyClassImpl()
    : numthreads(1)
    , unlockgil(false)
    , logger_verbose_(-1)
    , logger_(LoggerBase::emptyCallback())
    , debug_trace_(nullptr)
  {
  }
  ~ZgyClassImpl()
  {
    Py_XDECREF(debug_trace_);
    debug_trace_ = nullptr;
  }
};

struct ZgyClass
{
  PyObject_HEAD
  ZgyClassImpl *pimpl_;
};

static bool pimplcheck(ZgyClass* self)
{
  if (!self || !self->pimpl_ || !self->pimpl_->meta_) {
    _raise_simple_error(_zgy_error_InternalError,
      "Trying to access a closed ZgyReader or ZgyWriter");
    return false;
  }
  return true;
}

static PyObject * ZgyReader_open(ZgyClass* self, PyObject* args, PyObject* keywds);
static PyObject * ZgyWriter_create(ZgyClass* self, PyObject* args, PyObject* keywds);
static PyObject * ZgyWriter_clone(ZgyClass* self, PyObject* args, PyObject* keywds);
static PyObject * ZgyUtils_open(ZgyClass* self, PyObject* args, PyObject* keywds);

static void
ZgyReader_dealloc(ZgyClass *self)
{
  ZgyClassImpl *pimpl = self->pimpl_;
  self->pimpl_ = NULL;
  pimpl->logger_(1, "Destructed a ZgyReader\n");
  delete pimpl;
  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
ZgyReader_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ZgyClass *self = reinterpret_cast<ZgyClass*>(type->tp_alloc(type, 0));
  self->pimpl_ = new ZgyClassImpl();
  self->pimpl_->logger_verbose_ = LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE");
  self->pimpl_->logger_ = LoggerBase::standardCallback(self->pimpl_->logger_verbose_, "openzgy-pyreader: ", "");
  self->pimpl_->logger_(1, "Created a new ZgyReader\n");
  return (PyObject *)self;
}

static int
ZgyReader_init(ZgyClass *self, PyObject *args, PyObject *kwds)
{
  // TODO-Low should I call the base class init now? But the base is just Object.
  PyObject* dummy = ZgyReader_open(self, args, kwds);
  int ret = dummy == NULL ? -1 : 0;
  Py_XDECREF(dummy);
  self->pimpl_->logger_(1, "Initialized a new ZgyReader\n");
  return ret;
}

static void
ZgyWriter_dealloc(ZgyClass *self)
{
  ZgyClassImpl *pimpl = self->pimpl_;
  self->pimpl_ = NULL;
  pimpl->logger_(1, "Destructed a ZgyWriter\n");
  delete pimpl;
  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
ZgyWriter_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ZgyClass *self = reinterpret_cast<ZgyClass*>(type->tp_alloc(type, 0));
  self->pimpl_ = new ZgyClassImpl();

  self->pimpl_->logger_verbose_ = LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE");
  self->pimpl_->logger_ = LoggerBase::standardCallback(self->pimpl_->logger_verbose_, "openzgy-pywriter: ", "");
  self->pimpl_->logger_(1, "Created a new ZgyWriter\n");
  return (PyObject *)self;
}

/**
 * Check if the arguments to the ZgyWriter constructor looks like
 * the "clone" version. If not, we will assume it is the regular create.
 * This should be unambiguous: A clone requires either the "template"
 * keyword argument or a string as the second arg. A regular create
 * cannot have the "template" keyword and the second arg should be a
 * tuple.
 */
static bool
callerWantsClone(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  // Look for a "template" keyword argument.
  if (keywds && PyMapping_HasKeyString(keywds, "templatename")) {
    //std::cout << "Caller wants clone (1)" << std::endl;
    return true;
  }

  // Look for a string in the second argument.
  if (args) {
    PyObject* item = PySequence_GetItem(args, 1);
    if (item) {
      if (PyUnicode_Check(item)) {
        Py_XDECREF(item);
        //std::cout << "Caller wants clone (2)" << std::endl;
        return true;
      }
      else {
        // Not a string, probably a 3-tuple.
        //std::cout << "Caller doesn't want clone (3)" << std::endl;
      }
      Py_XDECREF(item);
    }
    else {
      // Probably less than two positional arguments.
      //std::cout << "Caller doesn't want clone (4)" << std::endl;
      PyErr_Clear();
    }
  }
  return false;
}

static int
ZgyWriter_init(ZgyClass *self, PyObject *args, PyObject *kwds)
{
  // TODO-Low should I call the base class init now? But the base is just Object.
  PyObject* dummy = 0;
  if (!callerWantsClone(self, args, kwds))
    dummy = ZgyWriter_create(self, args, kwds);
  else
    dummy = ZgyWriter_clone(self, args, kwds);

  int ret = dummy == NULL ? -1 : 0;
  Py_XDECREF(dummy);
  self->pimpl_->logger_(1, "Initialized a new ZgyWriter\n");
  return ret;
}

static void
ZgyUtils_dealloc(ZgyClass *self)
{
  ZgyClassImpl *pimpl = self->pimpl_;
  self->pimpl_ = NULL;
  pimpl->logger_(1, "Destructed a ZgyUtils\n");
  delete pimpl;
  Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject *
ZgyUtils_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  ZgyClass *self = reinterpret_cast<ZgyClass*>(type->tp_alloc(type, 0));
  self->pimpl_ = new ZgyClassImpl();
  self->pimpl_->logger_verbose_ = LoggerBase::getVerboseFromEnv("OPENZGY_VERBOSE");
  self->pimpl_->logger_ = LoggerBase::standardCallback(self->pimpl_->logger_verbose_, "openzgy-pyutils: ", "");
  self->pimpl_->logger_(1, "Created a new ZgyUtils\n");
  return (PyObject *)self;
}

static int
ZgyUtils_init(ZgyClass *self, PyObject *args, PyObject *kwds)
{
  PyObject* dummy = ZgyUtils_open(self, args, kwds);
  int ret = dummy == NULL ? -1 : 0;
  Py_XDECREF(dummy);
  self->pimpl_->logger_(1, "Initialized a new ZgyUtils\n");
  return ret;
}

/////////////////////////////////////////////////////////////////////////////
///   ZgyReader methods   ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

static PyObject *
ZgyCommon_enter(ZgyClass* self, PyObject* args)
{
  Py_INCREF(self);
  return (PyObject*)self;
}

static bool
ZgyCommon_close_all(ZgyClass *self, OpenZGY::FinalizeAction action = OpenZGY::FinalizeAction::BuildDefault, bool force = false)
{
  if (self->pimpl_->writer_) {
    self->pimpl_->logger_(1, "Closing and finalizing \"" + self->pimpl_->filename + "\"\n");
    try {
      // Default decimation and no progress bar.
      // Application should call finalize() directly for better control.
      self->pimpl_->writer_->finalize
        (std::vector<DecimationType>{DecimationType::LowPass,
                                     DecimationType::WeightedAverage},
          [](std::int64_t,std::int64_t){ return true; },
          action, force);
      self->pimpl_->writer_->close();
      self->pimpl_->writer_.reset();
    }
    catch (...) {
      _raise_ex();
      return false;
    }
  }

  try {
    if (self->pimpl_->reader_) {
      self->pimpl_->logger_(1, "Closing read-only file \"" + self->pimpl_->filename + "\"\n");
      self->pimpl_->reader_->close();
      self->pimpl_->reader_.reset();
    }
  }
  catch (...) {
    _raise_ex();
    return false;
  }

  // It is not really a good idea to allow reading metadata from a closed
  // accessor. This can cause other instances to be held on to too long.
  // Reasons to do so anyway: (1) there might be applicatoin code that
  // depends on this behavior, and (2) resetting meta_ here means I need
  // a null pointer check every place meta_ is used. I might have missed some.
  self->pimpl_->meta_.reset();
  self->pimpl_->filename = "";

  self->pimpl_->utils_.reset();

  return true;
}

static PyObject *
ZgyCommon_exit(ZgyClass* self, PyObject* args)
{
  if (!ZgyCommon_close_all(self))
    return NULL;
  Py_INCREF(Py_False);
  return Py_False; // If exiting due to an exception, do not swallow it.
}

/**
 * Get something that is either an attribute or a dict-like entry, i.e.
 * o["attrname"] or o.attrname. Reurns nullptr WITHOUT an exception set
 * if not found. Returns a new reference otherwise.
 */
static PyObject*
getOneAttr(PyObject* obj, const char *attrname)
{
  if (!obj || obj == Py_None || !attrname || !*attrname)
    return nullptr;
  PyObject* attr = PyUnicode_FromString(attrname);
  if (attr == NULL) {
    PyErr_Clear();
    return nullptr;
  }
  PyObject* val = PyObject_GetItem(obj, attr);
  //if (val) printf("Found o[\"%s\"]\n", attrname);
  if (val == NULL) {
    PyErr_Clear();
    val = PyObject_GetAttr(obj, attr);
    //if (val) printf("Found o.%s\n", attrname);
  }
  Py_DECREF(attr);
  if (val == NULL) {
    PyErr_Clear();
    return nullptr;
  }
  return val;
}

static std::string
getStringValuedAttr(PyObject* obj, const char *attrname)
{
  if (!obj || obj == Py_None || !attrname || !*attrname)
    return std::string();
  PyObject* val = getOneAttr(obj, attrname);
  if (val == NULL) {
    // attribute missing, I default it to the empty string.
    PyErr_Clear();
    return std::string();
  }
  const char *str = PyUnicode_AsUTF8(val);
  if (!str) {
    // Attribute might have been the wrong type?
    // Ideally I should propagate this error.
    PyErr_Clear();
    Py_DECREF(val);
    return std::string();
  }
  std::string result(str);
  Py_DECREF(val);
  //printf("Found %s = \"%s\"\n", attrname, result.c_str());
  return result;
}

static std::int64_t
getLongValuedAttr(PyObject* obj, const char *attrname, std::int64_t dflt)
{
  if (!obj || obj == Py_None || !attrname || !*attrname)
    return dflt;
  PyObject* val = getOneAttr(obj, attrname);
  if (val == NULL) {
    // Probably the attribute did not exist.
    PyErr_Clear();
    return dflt;
  }
  long long result = PyLong_AsLongLong(val);
  Py_DECREF(val);
  if (result == -1 && PyErr_Occurred()) {
    // Actually an error; should I raise it?
    PyErr_Clear();
    return dflt;
  }
  return static_cast<std::int64_t>(result);
}

static std::shared_ptr<IOContext>
ZgyCommon_getIOContext(ZgyClass *self, const char *filename, PyObject* obj)
{
  // iocontext is now expected to be a dict,
  // but the code supports both attribute and dict lookup.
  // TODO-Low: Warning if unrecognized attributes are present.
  if (filename && strncmp(filename, "sd://", 5) == 0) {
    auto result = std::make_shared<SeismicStoreIOContext>();
    std::string value;
    std::int64_t ivalue;
    if (!(value = getStringValuedAttr(obj, "sdurl")).empty())
      result->sdurl(value);
    if (!(value = getStringValuedAttr(obj, "sdapikey")).empty())
      result->sdapikey(value);
    if (!(value = getStringValuedAttr(obj, "sdtoken")).empty())
      result->sdtoken(value, std::string());
    if ((ivalue = getLongValuedAttr(obj, "maxsize", -1)) != -1)
      result->maxsize((int)ivalue);
    if ((ivalue = getLongValuedAttr(obj, "maxhole", -1)) != -1)
      result->maxhole((int)ivalue);
    if ((ivalue = getLongValuedAttr(obj, "aligned", -1)) != -1)
      result->aligned((int)ivalue);
    if ((ivalue = getLongValuedAttr(obj, "segsize", -1)) != -1)
      result->segsize((int)ivalue);
    if ((ivalue = getLongValuedAttr(obj, "segsplit", -1)) != -1)
      result->segsplit((int)ivalue);
    if ((ivalue = getLongValuedAttr(obj, "iothreads", -1)) != -1)
      result->iothreads((int)ivalue);
    if ((ivalue = getLongValuedAttr(obj, "cputhreads", -1)) != -1)
      result->cputhreads((int)ivalue);
    if (!(value = getStringValuedAttr(obj, "legaltag")).empty())
      result->legaltag(value);
    if (!(value = getStringValuedAttr(obj, "writeid")).empty())
      result->writeid(value);
    if (!(value = getStringValuedAttr(obj, "seismicmeta")).empty())
      result->seismicmeta(value);
    // TODO-Low: The pure Python api has seismicmeta as a dict. Do likewise.
    if (obj) {
      PyObject *callback = getOneAttr(obj, "_debug_trace");
      // TODO-Low: The callback will be leaked because I am not storing
      // it in the pimpl struct. No big deal since this is only for
      // debugging.
      if (!callback) {
        self->pimpl_->logger_(9, "debug_trace callback was not provided.\n");
        PyErr_Clear();
      }
      else {
        self->pimpl_->logger_(2, "debug_trace callback is being registered.\n");
        // Note that the saved self->pimpl_->debug_trace_ is only used
        // to make sure the python callable is released when the file
        // is closed. The function captures the callable itself.
        // TODO-Worry: Consider leaking the callable just in case
        // I messed up and end up calling it after going out of scope,
        PyObject* old_callback = self->pimpl_->debug_trace_;
        self->pimpl_->debug_trace_ = callback;
        Py_XDECREF(old_callback);
        auto logger = self->pimpl_->logger_;
        result->debug_trace
          ([callback, logger]
           (const std::string& operation,
            std::int64_t need,
            std::int64_t want,
            std::int64_t parts,
            const std::vector<std::int64_t>& sizes)
           {
             logger(2, "debug_trace callback is being invoked.\n");
             // Ignore errors creating and populating the tuple and
             // working with numbers. How could those fail?
             // Note that SetItem steals the reference it is given.
             PyObject* sizes_obj = PyTuple_New(sizes.size());
             for (std::size_t ii = 0; ii < sizes.size(); ++ii)
               PyTuple_SetItem(sizes_obj, ii,
                               PyLong_FromLongLong((long long)sizes[ii]));
             PyObject* rc = PyObject_CallFunction
               (callback, "sLLLN", operation.c_str(),
                (long long)need, (long long)want, (long long)parts, sizes_obj);
             if (!rc) {
               logger(0, "Warning: debug_trace callback raised an error.\n");
               PyErr_Clear(); // ignore errors raised by callback.
             }
             else {
               Py_DECREF(rc);
             }
           });
      }
    }
    return result;
  }
  return nullptr;
}

static PyObject *
ZgyReader_open(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  static char *kwlist[] = {
    const_cast<char*>("filename"),
    const_cast<char*>("iocontext"),
    const_cast<char*>("unlockgil"),
    NULL
  };

  char *filename{nullptr};
  PyObject *iocontext_obj{nullptr};
  int unlockgil{0};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|O$p", kwlist,
                                   &filename, &iocontext_obj, &unlockgil,
                                   NULL))
    return NULL;

  try {
    std::shared_ptr<IOContext> iocontext =
      ZgyCommon_getIOContext(self, filename, iocontext_obj);

    self->pimpl_->logger_(1, "Open \"" + std::string(filename) + "\"\n");
    if (iocontext)
      self->pimpl_->logger_(1, iocontext->toString());
    if (unlockgil > 0)
      self->pimpl_->logger_(1, "UnlockGIL set to True\n");

    self->pimpl_->reader_ = IZgyReader::open(std::string(filename), iocontext.get());
    self->pimpl_->meta_ = self->pimpl_->reader_;
    self->pimpl_->filename = filename;
    self->pimpl_->numthreads = 1;
    self->pimpl_->unlockgil = (unlockgil > 0);
  }
  catch (...) {
    return _raise_ex();
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyUtils_open(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  static char *kwlist[] = {
    const_cast<char*>("prefix"),
    const_cast<char*>("iocontext"),
    NULL
  };

  char *prefix{nullptr};
  PyObject *iocontext_obj{nullptr};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|O", kwlist,
                                   &prefix, &iocontext_obj,
                                   NULL))
    return NULL;

  try {
    std::shared_ptr<IOContext> iocontext =
      ZgyCommon_getIOContext(self, prefix, iocontext_obj);

    self->pimpl_->logger_(1, "Utils \"" + std::string(prefix) + "\"\n");
    if (iocontext)
      self->pimpl_->logger_(1, iocontext->toString());

    self->pimpl_->utils_ = IZgyUtils::utils(prefix, iocontext.get());
    self->pimpl_->filename = prefix;
    self->pimpl_->numthreads = 1;
  }
  catch (...) {
    return _raise_ex();
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyReader_close(ZgyClass* self, PyObject* args)
{
  if (!ZgyCommon_close_all(self))
    return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyWriter_close(ZgyClass* self, PyObject* args)
{
  if (!ZgyCommon_close_all(self))
    return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyWriter_close_incomplete(ZgyClass* self, PyObject* args)
{
  if (!ZgyCommon_close_all(self, OpenZGY::FinalizeAction::Delete))
    return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyUtils_close(ZgyClass* self, PyObject* args)
{
  if (!ZgyCommon_close_all(self))
    return NULL;
  Py_INCREF(Py_None);
  return Py_None;
}

static void dumpBuffer(const Py_buffer& buf)
{
#if 0
  // https://docs.python.org/3/c-api/buffer.html
  const char *ro =
    buf.readonly == 0 ? "writable" :
    buf.readonly == 1 ? "readonly" : "unknown";
  printf("Python %s buffer ", ro);
  //printf("at 0x%llx ", (long long)buf.buf);
  //printf("obj 0x%llx ", (long long)buf.obj);
  printf("length = %lld * %lld ",
         (long long)buf.itemsize,
         (long long)(buf.len / buf.itemsize));
  if (buf.shape) {
    switch (buf.ndim) {
    case 0: printf("scalar "); break;
    case 1: printf("size = [%d] ", (int)buf.shape[0]); break;
    case 2: printf("size = [%d,%d] ", (int)buf.shape[0], (int)buf.shape[1]); break;
    case 3: printf("size = [%d,%d,%d] ", (int)buf.shape[0], (int)buf.shape[1], (int)buf.shape[2]); break;
    default: printf("ndim = %d ", buf.ndim); break;
    }
  }
  else {
    printf("ndim = %d. shape = (null) ", buf.ndim);
  }
  printf("format = %s\n", buf.format ? buf.format : "(null)");
#endif
}

static PyObject *
ZgyReader_read(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  int i0=0, j0=0, k0=0;
  Py_buffer bulk = {0};
  PyObject* bulkobj = NULL;
  int lod=0;
  PyObject* verbose_obj{nullptr};

  std::shared_ptr<IZgyReader> accessor = self->pimpl_->reader_;
  if (!accessor) {
    return _raise_simple_error(_zgy_error_UserError, "Not open for read");
  }

  static char *kwlist[] = {
    const_cast<char*>(""), // positional-only
    const_cast<char*>(""), // introduced in Python 3.6
    const_cast<char*>("lod"),
    const_cast<char*>("verbose"), // ignored
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "(iii)O|i$O", kwlist,
                                   &i0, &j0, &k0, &bulkobj, &lod, &verbose_obj,
                                   NULL))
    return NULL;

  if (PyObject_GetBuffer(bulkobj, &bulk, PyBUF_FORMAT|PyBUF_ND|PyBUF_WRITABLE) < 0)
    return NULL;

  dumpBuffer(bulk);

  if (bulk.ndim != 3 || bulk.shape == NULL || bulk.readonly) {
    PyBuffer_Release(&bulk);
    return _raise_simple_error(PyExc_TypeError, "Expected a 3-dimensional writable array");
  }

  const IZgyWriter::size3i_t start{i0, j0, k0};
  const IZgyWriter::size3i_t count{bulk.shape[0], bulk.shape[1], bulk.shape[2]};
  const char format =
    bulk.format == NULL || bulk.format[0] == '\0' || bulk.format[1] != '\0' ?
    '?' : bulk.format[0];

  PyThreadState *_save = nullptr;
  bool unlocked_gil = (self->pimpl_->unlockgil && self->pimpl_->debug_trace_ == nullptr);
  if (unlocked_gil) {
    // Release the Python GIL if unlockgil is specified and we do not have a debug_trace_ callback declared
    // The GIL is reacquired immediately before any Py* function calls / returns.
    Py_UNBLOCK_THREADS
  }

  try {
    if (count[0]==0 || count[1] == 0 || count[2]==0) {
      // nothing to be done.
    }
    else if (bulk.itemsize == sizeof(float) && format == 'f') {
      accessor->read(start, count, (float*)bulk.buf, lod);
    }
    else if (bulk.itemsize == sizeof(signed short) && format == 'h') {
      accessor->read(start, count, (signed short*)bulk.buf, lod);
    }
    else if (bulk.itemsize == sizeof(signed char) && format == 'b') {
      accessor->read(start, count, (signed char*)bulk.buf, lod);
    }
    else {
      if (unlocked_gil) {
        Py_BLOCK_THREADS
      }
      PyBuffer_Release(&bulk);
      return _raise_simple_error(PyExc_TypeError, "Expected array of float, short, or byte");
    }
  }
  catch (...) {
    if (unlocked_gil) {
      Py_BLOCK_THREADS
    }
    PyBuffer_Release(&bulk);
    return _raise_ex();
  }

  if (unlocked_gil) {
    Py_BLOCK_THREADS
  }
  PyBuffer_Release(&bulk);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyReader_readconst(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  int i0{0}, j0{0}, k0{0};
  int ni{0}, nj{0}, nk{0};
  int lod{0};
  int as_float{1}; // PyArg_Parse.. returns bools as "int".
  PyObject* verbose_obj{nullptr};

  if (!self->pimpl_->reader_) {
    return _raise_simple_error(_zgy_error_UserError, "Not open for read");
  }

  //def readconst(self, start, count, *, lod=0, as_float=True, verbose=None)

  static char *kwlist[] = {
    const_cast<char*>("start"),
    const_cast<char*>("size"),
    const_cast<char*>("lod"),
    const_cast<char*>("as_float"),
    const_cast<char*>("verbose"), // ignored
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "(iii)(iii)|ip$O", kwlist,
                                   &i0, &j0, &k0,
                                   &ni, &nj, &nk,
                                   &lod, &as_float, &verbose_obj,
                                   NULL))
    return NULL;

  try {
    const IZgyWriter::size3i_t start{i0, j0, k0};
    const IZgyWriter::size3i_t count{ni, nj, nk};
    std::shared_ptr<IZgyReader> accessor = self->pimpl_->reader_;
    const std::pair<bool,double> constant =
      accessor->readconst(start, count, lod, as_float!=0);
    if (!constant.first) {
      Py_INCREF(Py_None);
      return Py_None;
    }
    else {
      // TODO-Low should I convert to int if file is integral and !as_float?
      return Py_BuildValue("d", constant.second);
    }
  }
  catch (...) {
    return _raise_ex();
  }
}

/**
 * copy/paste the version in ZgyReader except:
 *  - no "lod" parameter.
 *  - uses writer_ instead of reader_.
 * It looks like making a shared version would be too messy.
 */
static PyObject *
ZgyWriter_read(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  int i0=0, j0=0, k0=0;
  Py_buffer bulk = {0};
  PyObject* bulkobj = NULL;
  PyObject* verbose_obj{nullptr};

  std::shared_ptr<IZgyWriter> accessor = self->pimpl_->writer_;
  if (!accessor) {
    return _raise_simple_error(_zgy_error_UserError, "Not open for read/write");
  }

  static char *kwlist[] = {
    const_cast<char*>(""), // positional-only
    const_cast<char*>(""), // introduced in Python 3.6
    const_cast<char*>("verbose"), // ignored
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "(iii)O|$O", kwlist,
                                   &i0, &j0, &k0, &bulkobj, &verbose_obj,
                                   NULL))
    return NULL;

  if (PyObject_GetBuffer(bulkobj, &bulk, PyBUF_FORMAT|PyBUF_ND|PyBUF_WRITABLE) < 0)
    return NULL;

  dumpBuffer(bulk);

  if (bulk.ndim != 3 || bulk.shape == NULL || bulk.readonly) {
    PyBuffer_Release(&bulk);
    return _raise_simple_error(PyExc_TypeError, "Expected a 3-dimensional writable array");
  }

  const IZgyWriter::size3i_t start{i0, j0, k0};
  const IZgyWriter::size3i_t count{bulk.shape[0], bulk.shape[1], bulk.shape[2]};
  const char format =
    bulk.format == NULL || bulk.format[0] == '\0' || bulk.format[1] != '\0' ?
    '?' : bulk.format[0];

  try {
    if (count[0]==0 || count[1] == 0 || count[2]==0) {
      // nothing to be done.
    }
    else if (bulk.itemsize == sizeof(float) && format == 'f') {
      accessor->read(start, count, (float*)bulk.buf);
    }
    else if (bulk.itemsize == sizeof(signed short) && format == 'h') {
      accessor->read(start, count, (signed short*)bulk.buf);
    }
    else if (bulk.itemsize == sizeof(signed char) && format == 'b') {
      accessor->read(start, count, (signed char*)bulk.buf);
    }
    else {
      PyBuffer_Release(&bulk);
      return _raise_simple_error(PyExc_TypeError, "Expected array of float, short, or byte");
    }
  }
  catch (...) {
    PyBuffer_Release(&bulk);
    return _raise_ex();
  }

  PyBuffer_Release(&bulk);
  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * copy/paste the version in ZgyReader except:
 *  - no "lod" parameter.
 *  - uses writer_ instead of reader_.
 *  - as_float made keyword_only to avoid confusion.
 * It looks like making a shared version would be too messy.
 */
static PyObject *
ZgyWriter_readconst(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  int i0{0}, j0{0}, k0{0};
  int ni{0}, nj{0}, nk{0};
  int as_float{1}; // PyArg_Parse.. returns bools as "int".
  PyObject* verbose_obj{nullptr};

  std::shared_ptr<IZgyWriter> accessor = self->pimpl_->writer_;
  if (!accessor) {
    return _raise_simple_error(_zgy_error_UserError, "Not open for read/write");
  }

  //def readconst(self, start, count, *, lod=0, as_float=True, verbose=None)

  static char *kwlist[] = {
    const_cast<char*>("start"),
    const_cast<char*>("size"),
    const_cast<char*>("as_float"),
    const_cast<char*>("verbose"), // ignored
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "(iii)(iii)|$pO", kwlist,
                                   &i0, &j0, &k0,
                                   &ni, &nj, &nk,
                                   &as_float, &verbose_obj,
                                   NULL))
    return NULL;

  try {
    const IZgyWriter::size3i_t start{i0, j0, k0};
    const IZgyWriter::size3i_t count{ni, nj, nk};
    const std::pair<bool,double> constant =
      accessor->readconst(start, count, as_float!=0);
    if (!constant.first) {
      Py_INCREF(Py_None);
      return Py_None;
    }
    else {
      // TODO-Low should I convert to int if file is integral and !as_float?
      return Py_BuildValue("d", constant.second);
    }
  }
  catch (...) {
    return _raise_ex();
  }
}

static PyObject *
ZgyWriter_write(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  int i0=0, j0=0, k0=0;
  Py_buffer bulk = {0};
  PyObject* bulkobj = NULL;
  if (!self->pimpl_->writer_) {
    return _raise_simple_error(_zgy_error_UserError, "Not open for write");
  }

  static char *kwlist[] = {
#if PY_VERSION_HEX >= 0x03060000
    const_cast<char*>(""), // positional-only
    const_cast<char*>(""), // introduced in Python 3.6
#else
    const_cast<char*>("start"),
    const_cast<char*>("data"),
#endif
    const_cast<char*>("verbose"), // ignored
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "(iii)O|$O", kwlist,
                                   &i0, &j0, &k0, &bulkobj))
    return NULL;

  if (PyObject_GetBuffer(bulkobj, &bulk, PyBUF_FORMAT|PyBUF_ND) < 0)
    return NULL;

  dumpBuffer(bulk);

  if (bulk.ndim != 3 || bulk.shape == NULL) {
    PyBuffer_Release(&bulk);
    return _raise_simple_error(PyExc_TypeError, "Expected a 3-dimensional array");
  }

  const IZgyWriter::size3i_t start{i0, j0, k0};
  const IZgyWriter::size3i_t count{bulk.shape[0], bulk.shape[1], bulk.shape[2]};
  const char format =
    bulk.format == NULL || bulk.format[0] == '\0' || bulk.format[1] != '\0' ?
    '?' : bulk.format[0];
  std::shared_ptr<IZgyWriter> writer = self->pimpl_->writer_;

  PyThreadState *_save = nullptr;
  bool unlocked_gil = (self->pimpl_->unlockgil && self->pimpl_->debug_trace_ == nullptr);
  if (unlocked_gil) {
    // Release the Python GIL if unlockgil is specified and we do not have a debug_trace_ callback declared
    // The GIL is reacquired immediately before any Py* function calls / returns.
    Py_UNBLOCK_THREADS
  }

  try {
    if (bulk.itemsize == sizeof(float) && format == 'f') {
      writer->write(start, count, (const float*)bulk.buf);
    }
    else if (bulk.itemsize == sizeof(signed short) && format == 'h') {
      writer->write(start, count, (const signed short*)bulk.buf);
    }
    else if (bulk.itemsize == sizeof(signed char) && format == 'b') {
      writer->write(start, count, (const signed char*)bulk.buf);
    }
    else {
      if (unlocked_gil) {
        Py_BLOCK_THREADS
      }
      PyBuffer_Release(&bulk);
      return _raise_simple_error(PyExc_TypeError, "Expected array of float, short, or byte");
    }
  }
  catch (...) {
    if (unlocked_gil) {
      Py_BLOCK_THREADS
    }
    PyBuffer_Release(&bulk);
    return _raise_ex();
  }
  if (unlocked_gil) {
    Py_BLOCK_THREADS
  }
  PyBuffer_Release(&bulk);
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyWriter_writeconst(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  int i0=0, j0=0, k0=0;
  int ni=0, nj=0, nk=0;
  double value{0};
  int is_storage{0};
  int verbose{0};
  if (!self->pimpl_->writer_) {
    return _raise_simple_error(_zgy_error_UserError, "Not open for write");
  }

  // def writeconst(self, start, value, size, is_storage, *, verbose = None):

  static char *kwlist[] = {
    const_cast<char*>("start"),
    const_cast<char*>("value"),
    const_cast<char*>("size"),
    const_cast<char*>("is_storage"),
    const_cast<char*>("verbose"), // ignored
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "(iii)d(iii)p|$p", kwlist,
                                   &i0, &j0, &k0,
                                   &value,
                                   &ni, &nj, &nk,
                                   &is_storage,
                                   &verbose))
    return NULL;

  const IZgyWriter::size3i_t start{i0, j0, k0};
  const IZgyWriter::size3i_t count{ni, nj, nk};
  const SampleDataType dt = is_storage ?
    self->pimpl_->meta_->datatype() : SampleDataType::float32;
  std::shared_ptr<IZgyWriter> writer = self->pimpl_->writer_;

  // The Python API writes "float" (converting as needed) or "storage"
  // (always matching the type on file, might be "float").
  // The C++ API decides what to write based on the type of the
  // value that is passed in. Which may cause a runtime error
  // if trying to use the int8 overload to write to an int16 file
  // or vice versa.
  // TODO-Low throw if outside integer range.
  try {
    switch (dt) {
    case SampleDataType::float32: {
      const float typed_value{(float)value};
      writer->writeconst(start, count, &typed_value);
      break;
    }
    case SampleDataType::int16: {
      const std::int16_t typed_value{(std::int16_t)value};
      writer->writeconst(start, count, &typed_value);
      break;
    }
    case SampleDataType::int8: {
      const std::int8_t typed_value{(std::int8_t)value};
      writer->writeconst(start, count, &typed_value);
      break;
    }
    default: {
      return _raise_simple_error(_zgy_error_InternalError,
                                 "Unsupported data type on file");
    }
    }
  }
  catch (...) {
    return _raise_ex();
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyWriter_finalize(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  if (!self->pimpl_->writer_) {
    return _raise_simple_error(_zgy_error_UserError, "Not open for write");
  }

  static char *kwlist[] = {
    const_cast<char*>("decimation"),
    const_cast<char*>("progress"),
    const_cast<char*>("action"),
    const_cast<char*>("force"),
    const_cast<char*>("verbose"),
    NULL
  };

  PyObject *decimation_obj{nullptr};
  PyObject *progress_obj{nullptr};
  PyObject *action_obj{nullptr};
  int force{0}; // PyArg_Parse.. returns bools as "int".
  PyObject *verbose_obj{nullptr}; // currently ignored.
  if (!PyArg_ParseTupleAndKeywords(args, keywds, "|$OOOpO", kwlist,
                                   &decimation_obj,
                                   &progress_obj,
                                   &action_obj,
                                   &force,
                                   &verbose_obj,
                                   NULL))
    return NULL;

  if (decimation_obj || action_obj || force != 0) {
    // Warning, not Deprecation warning. This is already a behavior change.
    if (PyErr_WarnEx(PyExc_Warning,
                     "The \"decimation\", \"action\", and \"force\""
                     " arguments to finalize() are now ignored."
                     " Please specify that information in the ZgyWriter"
                     " constructor instead", /*stack_level=*/1) < 0)
    {
      return NULL;
    }
  }

#if 1
  // TODO-WIP-BrickedAPI: "decimation", "action", and "force" are all
  // deprecated, and if set in the call to finalize they will
  // (soon) be ignored. After things settle down, keep the warning but
  // remove the setting.
  std::vector<DecimationType> decimation;
  if (decimation_obj && decimation_obj != Py_None) {
    const std::size_t len = PyList_Size(decimation_obj);
    for (std::size_t ii = 0; ii < len; ++ii) {
      long value = _decodeEnum(PyList_GetItem(decimation_obj, ii),
                               _enum_DecimationType);
      if (value < 0)
        return NULL;
      decimation.push_back(static_cast<DecimationType>(value));
    }
  }
#endif

  std::function<bool(std::int64_t,std::int64_t)> progress;
  int mainpid = getpid();
  if (progress_obj && progress_obj != Py_None) {
    progress = [progress_obj,mainpid](std::int64_t done, std::int64_t total) {
                 const bool asynchronous = (mainpid != getpid());
                 //printf("** progress %d / %d %s\n",
                 //     (int)done, (int)total, asynchronous ? "(async)" : "");
                 // TODO-Worry: If the progress report can be invoked
                 // asynchronously then this will need a major rethink.
                 // For now just pretend the user didn't provide a callback.
                 if (asynchronous)
                   return true;
                 PyObject* ret = PyObject_CallFunction
                   (progress_obj, "ll", (long long)done, (long long)total);
                 if (!ret) {
                   // TODO-Low: propagate the error to the caller: Return
                   // false so the computation aborts, then cancel the
                   // ZgyAborted exception and somehow put back the current
                   // exception.
                   PyErr_Clear();
                   printf("Exception in progress callback\n");
                   return false;
                 }
                 int result = PyObject_IsTrue(ret);
                 Py_DECREF(ret);
                 return result > 0; // "-1" errors from IsTrue ignored.
               };
  }

#if 1
  // TODO-WIP-BrickedAPI: "decimation", "action", and "force" are all deprecated.
  // After things settle down, keep the warning but remove the setting.
  FinalizeAction action{FinalizeAction::BuildDefault};
  if (action_obj) {
    long value = _decodeEnum(action_obj, _enum_FinalizeAction);
    if (value < 0)
      return NULL;
    action = static_cast<FinalizeAction>(value);
  }
#endif

  self->pimpl_->logger_(1, "Finalize \"" + self->pimpl_->filename + "\"\n");
  //PyObject_Print(decimation_obj, stdout, 0); printf("\n");

  try {
    self->pimpl_->writer_->finalize(decimation, progress, action, force!=0);
  }
  catch (...) {
    return _raise_ex();
  }

  Py_INCREF(Py_None);
  return Py_None;
}

/**
 * Make a shallow copy of a dict-like object into a new dict.
 * A nullptr input (rare in this module) returns an empty dict.
 */
static PyObject*
ZgyWriter_shallowCopyDict(PyObject* dict)
{
  // The docs I read don't explicitly state that this works for any
  // dict-like objects. It does state this for PyDict_Merge.
  //return PyDict_Copy(dict);
  PyObject* result = PyDict_New();
  if (!result) {
    return NULL;
  }
  if (dict) {
    if (PyDict_Merge(result, dict, false) < 0) {
      Py_XDECREF(result);
      return NULL;
    }
  }
  return result;
}

/**
 * Split a dict of keyword arguments into one we recognize and one we don't.
 * Note the input kwargs WILL BE MODIFIED. Use ZgyWriter_shallowCopyDict()
 * if it is not clear whether this is ok.
 *
 * Returns a new reference that needs to be freed by the caller.
 */
static PyObject*
ZgyWriter_splitArguments(char* kwlist[], PyObject* kwargs)
{
  PyObject* my_kwargs = PyDict_New();
  if (!my_kwargs)
    return NULL;
  for (char** kw = kwlist; *kw != nullptr; ++kw) {
    if (PyMapping_HasKeyString(kwargs, *kw)) {
      PyObject* oo = PyMapping_GetItemString(kwargs, *kw); // oo is borrowed ref
      if (!oo) {
        Py_XDECREF(my_kwargs);
        return NULL;
      }
      if (PyDict_SetItemString(my_kwargs, *kw, oo) < 0) { // list gets ref to oo.
        Py_XDECREF(my_kwargs);
        return NULL;
      }
      PyMapping_DelItemString(kwargs, *kw); // removes ref owned by input list.
    }
  }
  return my_kwargs;
}

/**
 * Parse some of the keyword_only parameters from write, update, or
 * clone and store the result into the provided ZgyWriterArgsV3.
 * This function handle attributes introduced with ZgyWriterArgsV3 (or -2).
 * The result is a simple pass/fail boolean; it does not return a
 * PyObject* as is the convention in Python. It does, however,
 * set the current exception before returning false.
 *
 * The long list of possible attributes has been spit up only to
 * make the code more readable. Which of the ZgyWriter_parseXxxArgs
 * that a particular attribute will be processed by is up to the
 * whims of the developer.
 *
 * Only optional, keyword-only attributes that correspond to a
 * setting in ZgyWriterArgs can be handled using this pattern.
 *
 * The fuction might throw a C++ exception.
 *
 * There should be no net changes in reference counts.
 */
static bool
ZgyWriter_parseV3Args(PyObject* kwargs, ZgyWriterArgsV3& writer_args)
{
  static char *kwlist[] = {
    const_cast<char*>("historange"),  // (ff)
    const_cast<char*>("lodmode"),     // O
    const_cast<char*>("decimation"),  // O
    NULL
  };

  float     histomin{std::numeric_limits<float>::infinity()};
  float     histomax{-std::numeric_limits<float>::infinity()};
  PyObject* lodmode_obj{nullptr};
  PyObject* decimation_obj{nullptr};

  {
    PyObject *my_args = PyTuple_New(0);
    PyObject* my_kwargs = ZgyWriter_splitArguments(kwlist, kwargs);
    if (!my_args || !my_kwargs || !PyArg_ParseTupleAndKeywords
        (my_args, my_kwargs, "|$(ff)OO", kwlist,
         &histomin, &histomax, &lodmode_obj, &decimation_obj, NULL))
    {
      Py_XDECREF(my_kwargs);
      Py_XDECREF(my_args);
      return false;
    }
    Py_XDECREF(my_kwargs);
    Py_XDECREF(my_args);
  }

  // The xxx_obj variables are borrowed references.
  // There should be no more references we must remember to drop.

  // If unset, ignore. Else, even if set to min >= max, pass to C++.
  // The PyArgs function raises an error if not a 2-tuple of scalars.
  if (histomin < histomax)
    writer_args.historange(histomin, histomax);

  // If unset, ignore. If set to Null, error. Else pass to C++
  if (lodmode_obj) {
    std::vector<DecimationType> lodmode;
    long value = _decodeEnum(lodmode_obj, _enum_LodMode);
    if (value < 0) {
      return false;
    }
    writer_args.lodmode(static_cast<LodMode>(value));
  }

  // If unset or set to Null, ignore. If set to {} or {...}, pass to C++
  if (decimation_obj && decimation_obj != Py_None) {
    if (!PySequence_Check(decimation_obj)) {
      return _raise_simple_error(_zgy_error_UserError, "\"decimation\" needs to be a sequence");
    }
    std::vector<DecimationType> decimation;
    const std::size_t len = PySequence_Size(decimation_obj);
    for (std::size_t ii = 0; ii < len; ++ii) {
      long value = _decodeEnum(PySequence_GetItem(decimation_obj, ii),
                               _enum_DecimationType);
      if (value < 0) {
        return false;
      }
      decimation.push_back(static_cast<DecimationType>(value));
    }
    writer_args.decimation(decimation);
  }
  return true;
}

/**
 * Parse some of the keyword_only parameters from write, update,
 * or clone. Store the result into the provided ZgyWriterArgsV3.
 *
 * Returns false, with a Python exception set in the context, on error.
 * The fuction might also throw a C++ exception with no Python context.
 * There should be no changes in net reference counts.
 *
 * The long list of possible attributes has been spit up only to
 * make the code more readable. Which of the ZgyWriter_parseXxxArgs
 * that a particular attribute will be processed by is up to the
 * whims of the developer.
 *
 * Only optional, keyword-only attributes that correspond to a
 * setting in ZgyWriterArgs can be handled using this pattern.
 *
 */
static bool
ZgyWriter_parseAnnotArgs(PyObject* kwargs, ZgyWriterArgsV3& writer_args)
{
  static char *kwlist[] = {
    const_cast<char*>("zstart"),      // f
    const_cast<char*>("zinc"),        // f
    const_cast<char*>("annotstart"),  // (ff)
    const_cast<char*>("annotinc"),    // (ff)
    const_cast<char*>("corners"),     // ((dd)(dd)(dd)(dd))
    NULL
  };

  float zstart{0}, zinc{0};
  float ilstart{0}, xlstart{0};
  float ilinc{0}, xlinc{0};
  ZgyWriterArgsV3::corners_t corners{0,0,0,0,0,0,0,0};

  PyObject *my_args = PyTuple_New(0);
  PyObject* my_kwargs = ZgyWriter_splitArguments(kwlist, kwargs);
  if (!my_args || !my_kwargs || !PyArg_ParseTupleAndKeywords
      (my_args, my_kwargs, "|$ff(ff)(ff)((dd)(dd)(dd)(dd))", kwlist,
       &zstart, &zinc, &ilstart, &xlstart, &ilinc, &xlinc,
       &corners[0][0], &corners[0][1],
       &corners[1][0], &corners[1][1],
       &corners[2][0], &corners[2][1],
       &corners[3][0], &corners[3][1],
       NULL))
  {
    Py_XDECREF(my_kwargs);
    Py_XDECREF(my_args);
    return false;
  }

  if (PyMapping_HasKeyString(my_kwargs, "annotstart"))
    writer_args.ilstart(ilstart).xlstart(xlstart);

  if (PyMapping_HasKeyString(my_kwargs, "annotinc"))
    writer_args.ilinc(ilinc).xlinc(xlinc);

  if (PyMapping_HasKeyString(my_kwargs, "zstart"))
    writer_args.zstart(zstart);

  if (PyMapping_HasKeyString(my_kwargs, "zinc"))
    writer_args.zinc(zinc);

  if (PyMapping_HasKeyString(my_kwargs, "corners"))
    writer_args.corners(corners);

  Py_XDECREF(my_kwargs);
  Py_XDECREF(my_args);
  return true;
}

/**
 * Parse some of the keyword_only parameters from write, update,
 * or clone. Store the result into the provided ZgyWriterArgsV3.
 * See ZgyWriter_parseAnnotArgs() for more details.
 */
static bool
ZgyWriter_parseUnitArgs(PyObject* kwargs, ZgyWriterArgsV3& writer_args)
{
  static char *kwlist[] = {
    const_cast<char*>("zunitdim"),    // O
    const_cast<char*>("hunitdim"),    // O
    const_cast<char*>("zunitname"),   // s
    const_cast<char*>("hunitname"),   // s
    const_cast<char*>("zunitfactor"), // d
    const_cast<char*>("hunitfactor"), // d
    NULL
  };

  PyObject *zunitdim_obj{nullptr};
  PyObject *hunitdim_obj{nullptr};
  char *zunitname_str = const_cast<char*>("");
  char *hunitname_str = const_cast<char*>("");
  double zunitfactor{1.0}, hunitfactor{1.0};

  PyObject *my_args = PyTuple_New(0);
  PyObject* my_kwargs = ZgyWriter_splitArguments(kwlist, kwargs);
  if (!my_args || !my_kwargs || !PyArg_ParseTupleAndKeywords
      (my_args, my_kwargs, "|$OOssdd", kwlist,
       &zunitdim_obj,  &hunitdim_obj,
       &zunitname_str, &hunitname_str,
       &zunitfactor,   &hunitfactor,
       NULL))
  {
    Py_XDECREF(my_kwargs);
    Py_XDECREF(my_args);
    return false;
  }

  int has_z = ((PyMapping_HasKeyString(my_kwargs, "zunitdim")    ? 1 : 0) +
               (PyMapping_HasKeyString(my_kwargs, "zunitname")   ? 1 : 0) +
               (PyMapping_HasKeyString(my_kwargs, "zunitfactor") ? 1 : 0) +
               (zunitdim_obj ? 1 : 0));
  int has_h = ((PyMapping_HasKeyString(my_kwargs, "hunitdim")    ? 1 : 0) +
               (PyMapping_HasKeyString(my_kwargs, "hunitname")   ? 1 : 0) +
               (PyMapping_HasKeyString(my_kwargs, "hunitfactor") ? 1 : 0) +
               (hunitdim_obj ? 1 : 0));

  Py_XDECREF(my_kwargs);
  Py_XDECREF(my_args);
  my_kwargs = my_args = nullptr;

  // TODO-WIP-BrickedAPI: None of these are valid to change in a file
  // being opened for update, because the unit name is a variable
  // length string. Arguably that rarely makes sense anyway.
  // TODO-WIP-BrickedAPI: Maybe raise an exception here?
  // Or at least, make an explicit unit test for error handling.
  if (has_z == 4) {
      long value = _decodeEnum(zunitdim_obj, _enum_UnitDimension);
      if (value < 0)
        return NULL;
      writer_args.zunit(static_cast<UnitDimension>(value),
                        std::string(zunitname_str), zunitfactor);
    }
  else if (has_z != 0) {
    return _raise_simple_error(_zgy_error_UserError, "zunitdim, zunitname, and zunitfactor must all be given if one of them are." + std::to_string(has_z));
  }
  if (has_h == 4) {
    long value = _decodeEnum(hunitdim_obj, _enum_UnitDimension);
    if (value < 0)
      return NULL;
    writer_args.hunit(static_cast<UnitDimension>(value),
                      std::string(hunitname_str), hunitfactor);
  }
  else if (has_h != 0) {
    return _raise_simple_error(_zgy_error_UserError, "hunitdim, hunitname, and hunitfactor must all be given if one of them are." + std::to_string(has_h));
  }
  return true;
}

/**
 * Parse some of the keyword_only parameters from write, update,
 * or clone. Store the result into the provided ZgyWriterArgsV3.
 * See ZgyWriter_parseAnnotArgs() for more details.
 */
static bool
ZgyWriter_parseCompressArgs(PyObject* kwargs, ZgyWriterArgsV3& writer_args)
{
  static char *kwlist[] = {
    const_cast<char*>("compressor"),  // O
    const_cast<char*>("lodcompressor"), // O
    const_cast<char*>("zfp_compressor"), // f
    const_cast<char*>("zfp_lodcompressor"), // f
    NULL
  };

  PyObject* compressor_obj{nullptr};
  PyObject* lodcompressor_obj{nullptr};
  float zfp_snr = 0, zfp_lodsnr = 0;

  {
    PyObject *my_args = PyTuple_New(0);
    PyObject* my_kwargs = ZgyWriter_splitArguments(kwlist, kwargs);
    if (!my_args || !my_kwargs || !PyArg_ParseTupleAndKeywords
        (my_args, my_kwargs, "|$OOff", kwlist,
         &compressor_obj, &lodcompressor_obj, &zfp_snr, &zfp_lodsnr, NULL))
    {
      Py_XDECREF(my_kwargs);
      Py_XDECREF(my_args);
      return false;
    }
    Py_XDECREF(my_kwargs);
    Py_XDECREF(my_args);
  }

  // The pure Python API supports zfp_[lod]compressor and not [lod]compressor,
  // but that code is deprecated. The C++ wrapper *might* have supported
  // [lod]compressor but implementing compression in the Python code
  // doesn't make much sense.
  if (compressor_obj)
    return _raise_simple_error(_zgy_error_UserError, "The wrapper requires zfp_compressor(snr) instead of compressor");
  else if (zfp_snr > 0)
    writer_args.zfp_compressor(zfp_snr);

  if (lodcompressor_obj)
    return _raise_simple_error(_zgy_error_UserError, "The wrapper requires zfp_lodcompressor(snr) instead of lodcompressor");
  else if (zfp_lodsnr > 0)
    writer_args.zfp_lodcompressor(zfp_lodsnr);

  return true;
}

/**
 * Parse some of the keyword_only parameters from write, update,
 * or clone. Store the result into the provided ZgyWriterArgsV3.
 * See ZgyWriter_parseAnnotArgs() for more details.
 */
static bool
ZgyWriter_parseDataTypeArgs(PyObject* kwargs, ZgyWriterArgsV3& writer_args)
{
  static char *kwlist[] = {
    const_cast<char*>("datatype"),    // O
    const_cast<char*>("datarange"),   // (ff)
    NULL
  };

  PyObject *datatype_obj{nullptr};
  float datamin{0}, datamax{-1};

  PyObject *my_args = PyTuple_New(0);
  PyObject* my_kwargs = ZgyWriter_splitArguments(kwlist, kwargs);
  if (!my_args || !my_kwargs || !PyArg_ParseTupleAndKeywords
      (my_args, my_kwargs, "|$O(ff)", kwlist,
       &datatype_obj, &datamin, &datamax, NULL))
  {
    Py_XDECREF(my_kwargs);
    Py_XDECREF(my_args);
    return false;
  }
  if (PyMapping_HasKeyString(my_kwargs, "datarange")) {
    writer_args.datarange(datamin, datamax);
  }

  if (datatype_obj) {
    long value = _decodeEnum(datatype_obj, _enum_SampleDataType);
    if (value < 0)
      return false;
    writer_args.datatype(static_cast<SampleDataType>(value));
  }

  Py_XDECREF(my_kwargs);
  Py_XDECREF(my_args);
  return true;
}

static PyObject *
ZgyWriter_create(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  static char *kwlist[] = {
    const_cast<char*>("filename"),    // s
    const_cast<char*>("size"),        // (iii)
    const_cast<char*>("bricksize"),   // (iii)
    const_cast<char*>("iocontext"),   // O
    const_cast<char*>("update"), // p
    const_cast<char*>("unlockgil"), // p
    NULL
  };

  char *filename = NULL;
  int size[3]{0,0,0};
  int bricksize[3]{64,64,64};
  PyObject* iocontext_obj{nullptr};
  int update{0};
  int unlockgil{0};

  // Nearly duplicated between ZgyWriter_create and ZgyWriter_clone().
  // But it had been a lot worse without the refactored ZgyWriter_parseXxxArgs.
  // Also, the list of parseXxxArgs functions called might differ.
  if (keywds && !PyArg_ValidateKeywordArguments(keywds))
    return NULL;
  ZgyWriterArgsV3 writer_args;
  PyObject* cloned_keywds = ZgyWriter_shallowCopyDict(keywds);
  if (!cloned_keywds)
    return NULL;
  try {
    // Some of the parsing refactored into separate functions.
    if (!ZgyWriter_parseV3Args(cloned_keywds, writer_args) ||
        !ZgyWriter_parseAnnotArgs(cloned_keywds, writer_args) ||
        !ZgyWriter_parseUnitArgs(cloned_keywds, writer_args) ||
        !ZgyWriter_parseCompressArgs(cloned_keywds, writer_args) ||
        !ZgyWriter_parseDataTypeArgs(cloned_keywds, writer_args))
    {
      Py_XDECREF(cloned_keywds);
      return NULL;
    }
  }
  catch (...) {
    Py_XDECREF(cloned_keywds);
    return _raise_ex();
  }
  // End nearly duplicated code.

  // Now handle positional arguments and any keywords not handled above.
  // This will also raise an exception on unrecognized keywords.
  if (!PyArg_ParseTupleAndKeywords(args, cloned_keywds, "s(iii)|$(iii)Opp", kwlist,
                                   &filename,
                                   &size[0], &size[1], &size[2], // (iii)
                                   &bricksize[0], &bricksize[1], &bricksize[2], // (iii)
                                   &iocontext_obj,
                                   &update,
                                   &unlockgil,
                                   NULL))
  {
    Py_XDECREF(cloned_keywds);
    return NULL;
  }
  // The clone is no longer useful. It will contrain the named parameters
  // that were handled locally instead of some ZgyWriter_parseXxxArgs().
  Py_XDECREF(cloned_keywds);
  cloned_keywds = NULL;

  try {
    std::shared_ptr<IOContext> iocontext =
      ZgyCommon_getIOContext(self, filename, iocontext_obj);

    if (filename != nullptr)
      writer_args.filename(std::string(filename));

    // filename and size are declared as mandatory, and can be passed
    // either as positional or keyword arguments. With update=True
    // it is valid to skip the size. Caller can fake a missing size
    // by passing (0,0,0).
    if (size[0] != 0 || size[1] != 0 || size[2] != 0)
      writer_args.size(size[0], size[1], size[2]);

    if (PyMapping_HasKeyString(keywds, "iocontext"))
      writer_args.iocontext(iocontext.get());

    if (PyMapping_HasKeyString(keywds, "bricksize"))
      writer_args.bricksize(bricksize[0], bricksize[1], bricksize[2]);

    self->pimpl_->logger_(1, "Create \"" + std::string(filename) + "\"\n");
    if (self->pimpl_->logger_(2, "")) {
      std::stringstream ss;
      ss << "ZgyWriterArgs built in Python code:";
      writer_args.dump(ss);
      self->pimpl_->logger_(2, ss.str());
    }

    self->pimpl_->writer_ = update ?
      OpenZGY::IZgyWriter::reopen(writer_args) :
      OpenZGY::IZgyWriter::open(writer_args);
    self->pimpl_->meta_ = self->pimpl_->writer_;
    self->pimpl_->filename = std::string(filename);
    self->pimpl_->numthreads = 1;
    self->pimpl_->unlockgil = (unlockgil > 0);
  }
  catch (...) {
    return _raise_ex();
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyWriter_clone(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  static char *kwlist[] = {
    const_cast<char*>("filename"),    // s
    const_cast<char*>("templatename"),// s
    const_cast<char*>("datatype"),    // O
    const_cast<char*>("datarange"),   // (ff)
    const_cast<char*>("iocontext"),   // O
    const_cast<char*>("unlockgil"),   // p
    NULL
  };

  char *filename_str{nullptr};
  char *template_str{nullptr};
  PyObject *datatype_obj{nullptr};
  float datamin = 0;
  float datamax = 0;
  PyObject* iocontext_obj{nullptr};
  int unlockgil{0};

  // Nearly duplicated between ZgyWriter_create and ZgyWriter_clone().
  // But it had been a lot worse without the refactored ZgyWriter_parseXxxArgs.
  // Also, the list of parseXxxArgs functions called might differ.
  if (keywds && !PyArg_ValidateKeywordArguments(keywds))
    return NULL;
  ZgyWriterArgsV3 writer_args;
  PyObject* cloned_keywds = ZgyWriter_shallowCopyDict(keywds);
  if (!cloned_keywds)
    return NULL;
  try {
    // Some of the parsing refactored into separate functions.
    // NOTE! Chicken and egg problem: metafrom() will override
    // settings done here. For now, only handle parameters not
    // available in a ZGY file, and thus not touched by metafrom().
    if (!ZgyWriter_parseV3Args(cloned_keywds, writer_args) ||
        !ZgyWriter_parseCompressArgs(cloned_keywds, writer_args) ||
        false/*!ZgyWriter_parseDataTypeArgs(cloned_keywds, writer_args)*/)
    {
      Py_XDECREF(cloned_keywds);
      return NULL;
    }
  }
  catch (...) {
    Py_XDECREF(cloned_keywds);
    return _raise_ex();
  }
  // End nearly duplicated code.

  // Now handle positional arguments and any keywords not handled above.
  // This will also raise an exception on unrecognized keywords.
  if (!PyArg_ParseTupleAndKeywords(args, cloned_keywds, "ss|$O(ff)Op", kwlist,
                                   &filename_str, &template_str,
                                   &datatype_obj, &datamin, &datamax,
                                   &iocontext_obj,
                                   &unlockgil,
                                   NULL))
  {
    Py_XDECREF(cloned_keywds);
    return NULL;
  }
  // The clone is no longer useful. It will contrain the named parameters
  // that were handled locally instead of some ZgyWriter_parseXxxArgs().
  Py_XDECREF(cloned_keywds);
  cloned_keywds = NULL;

  static const auto iscloud = [](const char* fn) {
                                return fn && strncmp(fn, "sd://", 5) == 0;
                              };
  const bool has_iocontext = PyMapping_HasKeyString(keywds, "iocontext");
  const bool out_sdcontext = has_iocontext && iscloud(filename_str);
  const bool tpl_sdcontext = has_iocontext && iscloud(template_str);
  std::shared_ptr<IOContext> iocontext;
  try {
    if (out_sdcontext || tpl_sdcontext)
      iocontext = ZgyCommon_getIOContext(self, "sd://", iocontext_obj);
  }
  catch (...) {
    return _raise_ex();
  }

  try {
    auto parent = IZgyReader::open(std::string(template_str),
                                   tpl_sdcontext ? iocontext.get() : nullptr);
    writer_args.metafrom(parent);
  }
  catch (...) {
    return _raise_ex();
  }
  writer_args.filename(std::string(filename_str));
  if (out_sdcontext)
    writer_args.iocontext(iocontext.get());
  if (datamin < datamax)
    writer_args.datarange(datamin, datamax);
  if (datatype_obj) {
    long value = _decodeEnum(datatype_obj, _enum_SampleDataType);
    if (value < 0)
      return NULL;
    writer_args.datatype(static_cast<SampleDataType>(value));
  }

  self->pimpl_->logger_(1, "Clone \"" + std::string(filename_str) + "\" from \"" + std::string(template_str) + "\"\n");
  if (self->pimpl_->logger_(2, "")) {
    std::stringstream ss;
    ss << "ZgyWriterArgs built in Python code:";
    writer_args.dump(ss);
    self->pimpl_->logger_(2, ss.str());
  }

  try {
    self->pimpl_->writer_ = OpenZGY::IZgyWriter::open(writer_args);
    self->pimpl_->meta_ = self->pimpl_->writer_;
    self->pimpl_->filename = std::string(filename_str);
    self->pimpl_->numthreads = 1;
    self->pimpl_->unlockgil = (unlockgil > 0);
  }
  catch (...) {
    return _raise_ex();
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyUtils_delete(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  static char *kwlist[] = {
    const_cast<char*>("filename"),
    const_cast<char*>("missing_ok"),
    NULL
  };

  char *filename{nullptr};
  int missing_ok{1}; // PyArg_Parse.. returns bools as "int".

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|$p", kwlist,
                                   &filename, &missing_ok,
                                   NULL))
    return NULL;

  if (self->pimpl_->logger_(1, "")) {
    std::stringstream ss;
    ss << "ZgyUtils::delete("
       << "\"" << std::string(filename) << "\""
       << ", missing_ok=" << (missing_ok?"True":"False")
       << ")\n";
    self->pimpl_->logger_(1, ss.str());
  }

  try {
    self->pimpl_->utils_->deletefile(filename, missing_ok!=0);
  }
  catch (...) {
    return _raise_ex();
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyUtils_alturl(ZgyClass* self, PyObject* args, PyObject* keywds)
{
  static char *kwlist[] = {
    const_cast<char*>("filename"),
    NULL
  };

  char *filename{nullptr};

  if (!PyArg_ParseTupleAndKeywords(args, keywds, "s", kwlist,
                                   &filename,
                                   NULL))
    return NULL;

  if (self->pimpl_->logger_(1, "")) {
    std::stringstream ss;
    ss << "ZgyUtils::alturl("
       << "\"" << std::string(filename) << "\""
       << ")\n";
    self->pimpl_->logger_(1, ss.str());
  }

  std::string result;
  try {
    result = self->pimpl_->utils_->alturl(filename);
  }
  catch (...) {
    return _raise_ex();
  }
  return Py_BuildValue("s", result.c_str());
}

/*********************************************************************
 *
 * Purpose:
 *
 *      We have two coordinate systems (e.g. a grid coordinate
 *      system and a world (projection) coordinate system), but
 *      we don't know the general transformation between them.
 *
 *      What we do have is three arbitrary points, and their
 *      positions in both coordinate systems.
 *
 *      This routine uses these three points to transform
 *      points from one system to the other.
 *
 *
 * Description:
 *
 *      Let p0, p1, p2 be the 3 points, and let a = p1 - p0,
 *      and b = p2 - p0. An arbitrary point q can then be
 *      written as: q = p0 + s*a + t*b; where (s,t) are
 *      the "p"-coordinates of point q.
 *
 *      (s,t) are independant of which coordinate system
 *      we use to express p, so if we solve for s,t in
 *      one system, we can then directly compute q in
 *      the other system.
 *
 *      To see why (s,t) don't depend on system, write the
 *      equation as: q - p0 = s*a + t*b.
 *      Let T(v) be the transformation (some combination of
 *      scaling and rotation) that takes a point from A to B.
 *      A linear transformation is defined by these rules:
 *       1. T(cv) = cT(v)     - invariant under multiplication with a constant
 *       2. T(v+w) = T(v) + T(w) - invariant under vector addition
 *
 *      If we apply this, we get:
 *        T(q-p0) = T(s*a + t*b)
 *                = T(s*a) + T(t*b)   - rule 2
 *                = s*T(a) + t*T(b)   - rule 1
 *      Ie; (s,t) are the "p"-coordinates of q-p0 in both
 *      coordinate systems. QED.
 *
 *      If we expand the above equation, we get:
 *        xq - xp0 = s*xa + t*xb
 *        yq - yp0 = s*ya + t*yb
 *      We then have two equations in two unknowns, which
 *      we solve using Cramer's rule, and then use
 *      (s,t) to compute q in the other system.
 *
 *********************************************************************/
static bool ZgyReader_transform(
     double AX0, double AY0,
     double AX1, double AY1,
     double AX2, double AY2,
     double BX0, double BY0,
     double BX1, double BY1,
     double BX2, double BY2,
     double *X,  double* Y, std::int64_t length)
{
  // Make everything relative to p0
  AX1 -= AX0; AY1 -= AY0;
  AX2 -= AX0; AY2 -= AY0;
  BX1 -= BX0; BY1 -= BY0;
  BX2 -= BX0; BY2 -= BY0;

  double det = AX1*AY2 - AX2*AY1; // The determinant

  if(::fabs(det) < 1.0e-6) {
    _raise_simple_error(PyExc_ValueError, "Transform is not well defined due to colinear or coincident control points");
    return false;
  }

  for (std::int64_t ii = 0; ii < length; ++ii) {
    double xq = X[ii] - AX0 ;
    double yq = Y[ii] - AY0 ;
    double s  = (xq*AY2 - AX2*yq)/det ;
    double t  = (AX1*yq - xq*AY1)/det ;
    X[ii] = BX0 + s*BX1 + t*BX2 ;
    Y[ii] = BY0 + s*BY1 + t*BY2 ;
  }
  return true;
}

static PyObject *
ZgyReader_generalTransform(ZgyClass* self, PyObject* args)
{
  double ax0 = 0, ay0 = 0;
  double ax1 = 0, ay1 = 0;
  double ax2 = 0, ay2 = 0;
  double bx0 = 0, by0 = 0;
  double bx1 = 0, by1 = 0;
  double bx2 = 0, by2 = 0;
  PyObject* bulkobj = NULL;
  Py_buffer bulk = {0};

  if (!PyArg_ParseTuple(args, "((dd)(dd)(dd))((dd)(dd)(dd))O",
                        &ax0, &ay0, &ax1, &ay1, &ax2, &ay2,
                        &bx0, &by0, &bx1, &by1, &bx2, &by2,
                        &bulkobj))
    return NULL;

  if (PyObject_GetBuffer(bulkobj, &bulk, PyBUF_FORMAT|PyBUF_ND|PyBUF_WRITABLE) < 0)
    return NULL;

  dumpBuffer(bulk);

  if (bulk.ndim != 2 || bulk.shape == NULL || bulk.shape[0] != 2 || bulk.readonly) {
    PyBuffer_Release(&bulk);
    return _raise_simple_error(PyExc_TypeError, "Expected a 2d [2,*] writable array");
  }

  if (bulk.itemsize == sizeof(double) && bulk.format[0] == 'd') {
    const Py_ssize_t length = bulk.shape[1];
    double* data = (double*)bulk.buf;
    if (!ZgyReader_transform(
                             ax0, ay0, ax1, ay1, ax2, ay2,
                             bx0, by0, bx1, by1, bx2, by2,
                             data, data+length, length))
      return NULL;
  } else {
    // Could handle this if we wanted to, by copying into a temporary.
    return _raise_simple_error(PyExc_TypeError, "Expected an array of double");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject *
ZgyReader_indexToAnnot(ZgyClass* self, PyObject* args)
{
  if (!pimplcheck(self))
    return nullptr;

  std::array<double,2> point;
  if (!PyArg_ParseTuple(args, "(dd)", &point[0], &point[1]))
    return NULL;
  try {
    point = self->pimpl_->meta_->indexToAnnot(point);
  }
  catch (...) {
    return _raise_ex();
  }
  return Py_BuildValue("(dd)", point[0], point[1]);
}

static PyObject *
ZgyReader_annotToIndex(ZgyClass* self, PyObject* args)
{
  if (!pimplcheck(self))
    return nullptr;

  std::array<double,2> point;
  if (!PyArg_ParseTuple(args, "(dd)", &point[0], &point[1]))
    return NULL;
  try {
    point = self->pimpl_->meta_->annotToIndex(point);
  }
  catch (...) {
    return _raise_ex();
  }
  return Py_BuildValue("(dd)", point[0], point[1]);
}

static PyObject *
ZgyReader_indexToWorld(ZgyClass* self, PyObject* args)
{
  if (!pimplcheck(self))
    return nullptr;

  std::array<double,2> point;
  if (!PyArg_ParseTuple(args, "(dd)", &point[0], &point[1]))
    return NULL;
  try {
    point = self->pimpl_->meta_->indexToWorld(point);
  }
  catch (...) {
    return _raise_ex();
  }
  return Py_BuildValue("(dd)", point[0], point[1]);
}

static PyObject *
ZgyReader_worldToIndex(ZgyClass* self, PyObject* args)
{
  if (!pimplcheck(self))
    return nullptr;

  std::array<double,2> point;
  if (!PyArg_ParseTuple(args, "(dd)", &point[0], &point[1]))
    return NULL;
  try {
    point = self->pimpl_->meta_->worldToIndex(point);
  }
  catch (...) {
    return _raise_ex();
  }
  return Py_BuildValue("(dd)", point[0], point[1]);
}

static PyObject *
ZgyReader_annotToWorld(ZgyClass* self, PyObject* args)
{
  if (!pimplcheck(self))
    return nullptr;

  std::array<double,2> point;
  if (!PyArg_ParseTuple(args, "(dd)", &point[0], &point[1]))
    return NULL;
  try {
    point = self->pimpl_->meta_->annotToWorld(point);
  }
  catch (...) {
    return _raise_ex();
  }
  return Py_BuildValue("(dd)", point[0], point[1]);
}

static PyObject *
ZgyReader_worldToAnnot(ZgyClass* self, PyObject* args)
{
  if (!pimplcheck(self))
    return nullptr;

  std::array<double,2> point;
  if (!PyArg_ParseTuple(args, "(dd)", &point[0], &point[1]))
    return NULL;
  try {
    point = self->pimpl_->meta_->worldToAnnot(point);
  }
  catch (...) {
    return _raise_ex();
  }
  return Py_BuildValue("(dd)", point[0], point[1]);
}

/////////////////////////////////////////////////////////////////////////////
///   ZgyReader attributes   ////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// TODO-Low minor performance: Most attributes are immutable so we could
// have built the python values just once after reading the metadata.
// numthreads and verbose are notable exceptions.

static PyObject *
ZgyReader_getnumthreads(ZgyClass *self, void *closure)
{
  return Py_BuildValue("i", self->pimpl_->numthreads);
}

static int
ZgyReader_setnumthreads(ZgyClass *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError,
                    "The numthreads attribute cannot be deleted.");
    return -1;
  }
  if (!PyLong_Check(value)) {
    PyErr_SetString(PyExc_TypeError,
                    "The numthreads attribute must be an integer.");
    return -1;
  }
  long newvalue = PyLong_AsLong(value);
  if (newvalue < 1 || newvalue >= 1024) {
    PyErr_SetString(PyExc_AttributeError,
                    "The numthreads attribute must be between 1 and 1024.");
    return -1;
  }
  self->pimpl_->numthreads = static_cast<int>(newvalue);
  return 0;
}

static PyObject *
ZgyReader_getverbose(ZgyClass *self, void *closure)
{
  return Py_BuildValue("i", self->pimpl_->logger_verbose_);
}

static int
ZgyReader_setverbose(ZgyClass *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError,
                    "The verbose attribute cannot be deleted.");
    return -1;
  }
  if (!PyLong_Check(value)) {
    PyErr_SetString(PyExc_TypeError,
                    "The verbose attribute must be an integer.");
    return -1;
  }
  long newvalue = PyLong_AsLong(value);
  self->pimpl_->logger_verbose_ = static_cast<int>(newvalue);
  self->pimpl_->logger_ = LoggerBase::standardCallback(self->pimpl_->logger_verbose_, "openzgy-python: ", "");
  // Nitpick: not distinguishing between reader and writer in prefix,
  // as happens when logger is set from environment variable.
  // Because the same setter is used in both classes.
  return 0;
}

static PyObject *
ZgyReader_getfilename(ZgyClass *self, void *closure)
{
  return Py_BuildValue("s", self->pimpl_->filename.c_str());
}

static PyObject *
ZgyReader_getsize(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::array<std::int64_t,3> size = self->pimpl_->meta_->size();
  return Py_BuildValue("iii", size[0], size[1], size[2]);
}

static PyObject *
ZgyReader_getbricksize(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::array<std::int64_t,3> size = self->pimpl_->meta_->bricksize();
  return Py_BuildValue("iii", size[0], size[1], size[2]);
}

static PyObject *
ZgyReader_getbrickcount(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::vector<std::array<int64_t,3>> count = self->pimpl_->meta_->brickcount();
  PyObject* result = PyTuple_New(count.size());
  for (std::size_t lod=0; lod<count.size(); ++lod) {
    PyObject* one = Py_BuildValue("(iii)", count[lod][0], count[lod][1], count[lod][2]);
    if (!result || !one || PyTuple_SetItem(result, lod, one) != 0) {
      Py_XDECREF(one);
      Py_XDECREF(result);
      return nullptr;
    }
  }
  return result;
}

static PyObject *
ZgyReader_getdatatype(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  SampleDataType dt = self->pimpl_->meta_->datatype();
  return PyObject_CallFunction(_enum_SampleDataType, "i", (int)dt);
}

static PyObject *
ZgyReader_getdatarange(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::array<float,2> datarange = self->pimpl_->meta_->datarange();
  return Py_BuildValue("(ff)", datarange[0], datarange[1]);
}

static PyObject *
ZgyReader_getrawdatarange(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::array<float,2> datarange = self->pimpl_->meta_->raw_datarange();
  return Py_BuildValue("(ff)", datarange[0], datarange[1]);
}

static PyObject *
ZgyReader_getzunitdim(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  UnitDimension dim = self->pimpl_->meta_->zunitdim();
  return PyObject_CallFunction(_enum_UnitDimension, "i", (int)dim);
}

static PyObject *
ZgyReader_gethunitdim(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  UnitDimension dim = self->pimpl_->meta_->hunitdim();
  return PyObject_CallFunction(_enum_UnitDimension, "i", (int)dim);
}

static PyObject *
ZgyReader_getzunitname(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::string zunitname = self->pimpl_->meta_->zunitname();
  return Py_BuildValue("s", zunitname.c_str());
}

static PyObject *
ZgyReader_gethunitname(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::string hunitname = self->pimpl_->meta_->hunitname();
  return Py_BuildValue("s", hunitname.c_str());
}

static PyObject *
ZgyReader_getzunitfactor(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  return Py_BuildValue("d", self->pimpl_->meta_->zunitfactor());
}

static PyObject *
ZgyReader_gethunitfactor(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  return Py_BuildValue("d", self->pimpl_->meta_->hunitfactor());
}

static PyObject *
ZgyReader_getzstart(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  return Py_BuildValue("d", self->pimpl_->meta_->zstart());
}

static PyObject *
ZgyReader_getzinc(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  return Py_BuildValue("d", self->pimpl_->meta_->zinc());
}

static PyObject *
ZgyReader_getannotstart(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::array<float,2> annotstart = self->pimpl_->meta_->annotstart();
  return Py_BuildValue("(ff)", annotstart[0], annotstart[1]);
}

static PyObject *
ZgyReader_getannotinc(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::array<float,2> annotinc = self->pimpl_->meta_->annotinc();
  return Py_BuildValue("(ff)", annotinc[0], annotinc[1]);
}

static PyObject *
ZgyReader_getcorners(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  const IZgyMeta::corners_t& ocp = self->pimpl_->meta_->corners();
  return Py_BuildValue("((dd)(dd)(dd)(dd))",
                       ocp[0][0], ocp[0][1],
                       ocp[1][0], ocp[1][1],
                       ocp[2][0], ocp[2][1],
                       ocp[3][0], ocp[3][1]);
}

static PyObject *
ZgyReader_getindexcorners(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  const IZgyMeta::corners_t& ocp = self->pimpl_->meta_->indexcorners();
  return Py_BuildValue("((LL)(LL)(LL)(LL))",
                       (long long)ocp[0][0], (long long)ocp[0][1],
                       (long long)ocp[1][0], (long long)ocp[1][1],
                       (long long)ocp[2][0], (long long)ocp[2][1],
                       (long long)ocp[3][0], (long long)ocp[3][1]);
}

static PyObject *
ZgyReader_getannotcorners(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  const IZgyMeta::corners_t& ocp = self->pimpl_->meta_->annotcorners();
  return Py_BuildValue("((dd)(dd)(dd)(dd))",
                       ocp[0][0], ocp[0][1],
                       ocp[1][0], ocp[1][1],
                       ocp[2][0], ocp[2][1],
                       ocp[3][0], ocp[3][1]);
}

static PyObject *
ZgyReader_getnlods(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  return Py_BuildValue("i", self->pimpl_->meta_->nlods());
}

//static PyObject *
//ZgyReader_getdataid(ZgyClass *self, void *closure)
//{
//  return Py_BuildValue("s", self->pimpl_->meta_->dataid().c_str());
//}

static PyObject *
ZgyReader_getverid(ZgyClass *self, void *closure)
{
  return Py_BuildValue("s", self->pimpl_->meta_->verid().c_str());
}

//static PyObject *
//ZgyReader_getprevid(ZgyClass *self, void *closure)
//{
//  return Py_BuildValue("s", self->pimpl_->meta_->previd().c_str());
//}

static PyObject *
ZgyReader_getmeta(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
   return Py_BuildValue(
     "{s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N,s:N}",
     "filename",    ZgyReader_getfilename(self, closure),
     "size",        ZgyReader_getsize(self, closure),
     "bricksize",   ZgyReader_getbricksize(self, closure),
     "datatype",    ZgyReader_getdatatype(self, closure),
     "datarange",   ZgyReader_getdatarange(self, closure),
     "rawdatarange",ZgyReader_getrawdatarange(self, closure),
     "zunitdim",    ZgyReader_getzunitdim(self, closure),
     "hunitdim",    ZgyReader_gethunitdim(self, closure),
     "zunitname",   ZgyReader_getzunitname(self, closure),
     "hunitname",   ZgyReader_gethunitname(self, closure),
     "zunitfactor", ZgyReader_getzunitfactor(self, closure),
     "hunitfactor", ZgyReader_gethunitfactor(self, closure),
     "zstart",      ZgyReader_getzstart(self, closure),
     "zinc",        ZgyReader_getzinc(self, closure),
     "annotstart",  ZgyReader_getannotstart(self, closure),
     "annotinc",    ZgyReader_getannotinc(self, closure),
     "corners",     ZgyReader_getcorners(self, closure));
}

static PyObject *
ZgyReader_getstatistics(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  const OpenZGY::SampleStatistics s = self->pimpl_->meta_->statistics();
    PyObject* args = Py_BuildValue
    ("Ldddd", (long long)s.cnt,
     (double)s.sum, (double)s.ssq,
     (double)s.min, (double)s.max);
  if (args == nullptr)
    return nullptr;
  PyObject* result = PyObject_Call(_namedtuple_Statistics, args, nullptr);
  Py_DECREF(args);
  return result;
}

static PyObject *
ZgyReader_gethistogram(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  const OpenZGY::SampleHistogram h = self->pimpl_->meta_->histogram();
  PyObject *bins_obj = PyList_New(h.bins.size());
  if (bins_obj == nullptr)
    return nullptr;
  for (std::size_t ii=0; ii<h.bins.size(); ++ii)
    PyList_SetItem(bins_obj, ii, PyLong_FromLongLong(h.bins[ii]));
  PyObject* args = Py_BuildValue
    ("LddO", (long long)h.samplecount,
     (double)h.minvalue, (double)h.maxvalue,
     bins_obj);
  Py_XDECREF(bins_obj); // Because "args" has a ref to it, or args is dead.
  if (args == nullptr)
    return nullptr;
  PyObject* result = PyObject_Call(_namedtuple_Histogram, args, nullptr);
  Py_XDECREF(args);
  return result;
}

/*
 * FileStatistics. All members are std::int64_t except the last three
 * which are double, bool, and vector<int64_t> respectively.
 *
 */
static PyObject *
ZgyReader_getfilestats(ZgyClass *self, void *closure)
{
  if (!pimplcheck(self))
    return nullptr;
  std::shared_ptr<const OpenZGY::FileStatistics> filestats =
    self->pimpl_->meta_->filestats();
  if (!filestats)
    return _raise_simple_error(_zgy_error_InternalError, "No filestats");
  const OpenZGY::FileStatistics& s = *filestats;
  PyObject* iscompressed_obj = PyBool_FromLong(s.isCompressed() ? 1 : 0);
  const std::vector<std::int64_t> sizes = s.segmentSizes();
  PyObject* sizes_obj = PyTuple_New(sizes.size());
  for (std::size_t ii = 0; ii < sizes.size(); ++ii)
    PyTuple_SetItem(sizes_obj, ii,
                    PyLong_FromLongLong((long long)sizes[ii]));
  PyObject* args = Py_BuildValue
    ("LLLLLLLLLLLLLLLLLLdOO",
     (long long)s.fileVersion(),
     (long long)s.fileSize(),
     (long long)s.headerSize(),
     (long long)s.dataStart(),
     (long long)s.alphaNormalCount(),
     (long long)s.alphaNormalSizePerEntry(),
     (long long)s.alphaCompressedCount(),
     (long long)s.alphaCcompressedSize(),
     (long long)s.alphaMissingCount(),
     (long long)s.alphaConstantCount(),
     (long long)s.brickNormalCount(),
     (long long)s.brickNormalSizePerEntry(),
     (long long)s.brickCompressedCount(),
     (long long)s.brickCompressedSize(),
     (long long)s.brickMissingCount(),
     (long long)s.brickConstantCount(),
     (long long)s.usedSize(),
     (long long)s.usedIfUncompressed(),
     (double)s.compressionFactor(),
     iscompressed_obj,
     sizes_obj);
  if (args == nullptr)
    return nullptr;
  PyObject* result = PyObject_Call(_namedtuple_FileStats, args, nullptr);
  Py_DECREF(args);
  Py_DECREF(iscompressed_obj);
  Py_DECREF(sizes_obj);
  return result;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

static PyMethodDef ZgyReader_methods[] = {
    {"__enter__", (PyCFunction)ZgyCommon_enter, METH_VARARGS,
     "Enter a \"with\" context"
    },
    {"__exit__", (PyCFunction)ZgyCommon_exit, METH_VARARGS,
     "Leave a \"with\" context"
    },
    // Hidden to force users to use the with statement
    //{"open", (PyCFunction)ZgyReader_open, METH_VARARGS | METH_KEYWORDS,
    // "Open an existing ZGY file for read or write."
    //},
    {"close", (PyCFunction)ZgyReader_close, METH_NOARGS,
     "Close the currently open file."
    },
    {"read", (PyCFunction)ZgyReader_read, METH_VARARGS | METH_KEYWORDS,
     "Read bulk data into a caller specified buffer.\n"
     "The buffer's type must be float, short, or char.\n"
     "Any file may be read as float. short and char\n"
     "require the file to be of exactly that type.\n"
     "Arguments: (i0,j0,k0), buffer, lod=0"
    },
    {"readconst", (PyCFunction)ZgyReader_readconst, METH_VARARGS | METH_KEYWORDS,
     "Return the scalar value that all samples in the range are set to,\n"
     "or None if the range is not known to have all-constant values.\n"
     "Note that there are some cases where the function returns None\n"
     "even if all samples are in fact equal.\n"
     "Arguments: (i0,j0,k0), (ni, nj, nk), lod=0, as_float=True"
    },
    {"transform", (PyCFunction)ZgyReader_generalTransform, METH_VARARGS,
     "Linear transformation of an array of double-precision coordinates.\n"
     "The coordinate systems to convert between are defined by\n"
     "three arbitrary points in the source system and the target.\n"
     "Arguments: ((ax0,ay0), (ax1,ay1), (ax2,ay2)),\n"
     "           ((bx0,by0), (bx1,by1), (bx2,by2)),\n"
     "           data\n"
     "where data is a 2d array of size (length, 2)"
    },
    {"indexToAnnot", (PyCFunction)ZgyReader_indexToAnnot, METH_VARARGS,
     "Convert ordinal to inline, crossline"
    },
    {"annotToIndex", (PyCFunction)ZgyReader_annotToIndex, METH_VARARGS,
     "Convert inline, crossline to ordinal"
    },
    {"indexToWorld", (PyCFunction)ZgyReader_indexToWorld, METH_VARARGS,
     "Convert ordinal to world X,Y"
    },
    {"worldToIndex", (PyCFunction)ZgyReader_worldToIndex, METH_VARARGS,
     "Convert world X,Y to ordinal"
    },
    {"annotToWorld", (PyCFunction)ZgyReader_annotToWorld, METH_VARARGS,
     "Convert inline, crossline to world X,Y"
    },
    {"worldToAnnot", (PyCFunction)ZgyReader_worldToAnnot, METH_VARARGS,
     "Convert world X,Y to inline, crossline"
    },
    {NULL}  /* Sentinel */
};

static PyMethodDef ZgyWriter_methods[] = {
    {"__enter__", (PyCFunction)ZgyCommon_enter, METH_VARARGS,
     "Enter a \"with\" context"
    },
    {"__exit__", (PyCFunction)ZgyCommon_exit, METH_VARARGS,
     "Leave a \"with\" context"
    },
    {"finalize", (PyCFunction)ZgyWriter_finalize, METH_VARARGS | METH_KEYWORDS,
     "Generate low resolution data and compute histogram and statistics."
    },
    {"close", (PyCFunction)ZgyWriter_close, METH_NOARGS,
     "Close the currently open file after calculating statistics and low resolution data."
    },
    {"close_incomplete", (PyCFunction)ZgyWriter_close_incomplete, METH_NOARGS,
     "Close the currently open file without statistics and low resolution data."
    },
    // Hidden to force users to use the with statement
    //{"create", (PyCFunction)ZgyWriter_create, METH_VARARGS|METH_KEYWORDS,
    // "Create a new ZGY file."
    //},
    // Hidden to force users to use the with statement
    //{"clone", (PyCFunction)ZgyWriter_clone, METH_VARARGS|METH_KEYWORDS,
    // "Create a new ZGY file similar to an existing one."
    //},
    {"read", (PyCFunction)ZgyWriter_read, METH_VARARGS | METH_KEYWORDS,
     "Read bulk data into a caller specified buffer.\n"
     "The buffer's type must be float, short, or char.\n"
     "Any file may be read as float. short and char\n"
     "require the file to be of exactly that type.\n"
     "Arguments: (i0,j0,k0), buffer"
    },
    {"readconst", (PyCFunction)ZgyWriter_readconst, METH_VARARGS | METH_KEYWORDS,
     "Return the scalar value that all samples in the range are set to,\n"
     "or None if the range is not known to have all-constant values.\n"
     "Note that there are some cases where the function returns None\n"
     "even if all samples are in fact equal.\n"
     "Arguments: (i0,j0,k0), (ni, nj, nk), as_float=True"
    },
    {"write", (PyCFunction)ZgyWriter_write, METH_VARARGS | METH_KEYWORDS,
     "Write bulk data. Type must be float, short, or char.\n"
     "float may be written to any file. short and char require\n"
     "the file to be of exactly that type.\n"
     "Arguments: (i0,j0,k0), buffer"
    },
    {"writeconst", (PyCFunction)ZgyWriter_writeconst, METH_VARARGS | METH_KEYWORDS,
     "Write all-constant bulk data. Input should be an array of\n"
     "length 1. Type must be float, short, or char.\n"
     "float may be written to any file. short and char require\n"
     "the file to be of exactly that type.\n"
     "Arguments: (i0,j0,k0), value, (ni, nj, nk), is_storage"
    },
    {NULL}  /* Sentinel */
};

static PyMethodDef ZgyUtils_methods[] = {
    {"__enter__", (PyCFunction)ZgyCommon_enter, METH_VARARGS,
     "Enter a \"with\" context"
    },
    {"__exit__", (PyCFunction)ZgyCommon_exit, METH_VARARGS,
     "Leave a \"with\" context"
    },
    {"close", (PyCFunction)ZgyUtils_close, METH_NOARGS,
     "Close the currently open file."
    },
    {"delete", (PyCFunction)ZgyUtils_delete, METH_VARARGS | METH_KEYWORDS,
     "Delete a file.\n"
    },
    {"alturl", (PyCFunction)ZgyUtils_alturl, METH_VARARGS | METH_KEYWORDS,
     "Return a url that might be opened faster.\n"
    },
    {NULL}  /* Sentinel */
};

// The reader and writer share the same C++ class implementation
// and the getters/setters are all the same. So this tablle
// can be shared between the two.

static PyGetSetDef ZgyCommon_getseters[] = {
  {const_cast<char*>("numthreads"),
     (getter)ZgyReader_getnumthreads, (setter)ZgyReader_setnumthreads,
     const_cast<char*>("How many threads to use when reading. Currently ignored. Use iocontext instead."),
     NULL},
  {const_cast<char*>("verbose"),
     (getter)ZgyReader_getverbose, (setter)ZgyReader_setverbose,
     const_cast<char*>("How many debug messages to output to std::cerr."),
     NULL},
  {const_cast<char*>("filename"),
     (getter)ZgyReader_getfilename, (setter)NULL,
     const_cast<char*>("Name of the currently open file."),
     NULL},
  {const_cast<char*>("size"),
     (getter)ZgyReader_getsize, (setter)NULL,
     const_cast<char*>("Size in (inline, crossline, vertical) dimension."),
     NULL},
  {const_cast<char*>("bricksize"),
     (getter)ZgyReader_getbricksize, (setter)NULL,
     const_cast<char*>("Brick size in (inline, crossline, vertical) dimension."),
     NULL},
  {const_cast<char*>("brickcount"),
     (getter)ZgyReader_getbrickcount, (setter)NULL,
     const_cast<char*>("Brick size for each LOD level."),
     NULL},
  {const_cast<char*>("datatype"),
     (getter)ZgyReader_getdatatype, (setter)NULL,
     const_cast<char*>("Stringly typed enum. One of 'int8', 'int16', 'float'."),
     NULL},
  {const_cast<char*>("datarange"),
     (getter)ZgyReader_getdatarange, (setter)NULL,
     const_cast<char*>("Minimum and maximum data value.\nServes as clipping values if datatype is integral."),
     NULL},
  {const_cast<char*>("raw_datarange"),
     (getter)ZgyReader_getrawdatarange, (setter)NULL,
     const_cast<char*>("Minimum and maximum data value before adjustment."),
     NULL},
  {const_cast<char*>("zunitdim"),
     (getter)ZgyReader_getzunitdim, (setter)NULL,
     const_cast<char*>("Stringly typed enum. One of 'time', 'length', 'arcangle', 'unknown'."),
     NULL},
  {const_cast<char*>("hunitdim"),
     (getter)ZgyReader_gethunitdim, (setter)NULL,
     const_cast<char*>("Stringly typed enum. One of 'time', 'length', 'arcangle', 'unknown'."),
     NULL},
  {const_cast<char*>("zunitname"),
     (getter)ZgyReader_getzunitname, (setter)NULL,
     const_cast<char*>("Name of vertical unit, e.g. \"ms\", \"m\" or \"ft\"."),
     NULL},
  {const_cast<char*>("hunitname"),
     (getter)ZgyReader_gethunitname, (setter)NULL,
     const_cast<char*>("Name of horizontal unit, e.g. \"m\" or \"ft\"."),
     NULL},
  {const_cast<char*>("zunitfactor"),
     (getter)ZgyReader_getzunitfactor, (setter)NULL,
     const_cast<char*>("Scaling factor of vertical unit, e.g. 0.001 for ms, 1.0 for m or 0.3048 for ft."),
     NULL},
  {const_cast<char*>("hunitfactor"),
     (getter)ZgyReader_gethunitfactor, (setter)NULL,
     const_cast<char*>("Scaling factor of horizontal unit, e.g. 1.0 for m or 0.3048 for ft."),
     NULL},
  {const_cast<char*>("zstart"),
     (getter)ZgyReader_getzstart, (setter)NULL,
     const_cast<char*>("Distance from surface/MSL to first sample, given in the vertical unit."),
     NULL},
  {const_cast<char*>("zinc"),
     (getter)ZgyReader_getzinc, (setter)NULL,
     const_cast<char*>("Sample interval, given in the vertical unit."),
     NULL},
  {const_cast<char*>("annotstart"),
     (getter)ZgyReader_getannotstart, (setter)NULL,
     const_cast<char*>("First inline and crossline numbers."),
     NULL},
  {const_cast<char*>("annotinc"),
     (getter)ZgyReader_getannotinc, (setter)NULL,
     const_cast<char*>("Inline and crossline number increments between adjacent sections of the cube."),
     NULL},
  {const_cast<char*>("corners"),
     (getter)ZgyReader_getcorners, (setter)NULL,
     const_cast<char*>(
     "Corner coordinates (i0,x0) (i1,x0), (i0,x1), (i1,x1)\n"
     "Where i0 and x0 are the lowest numbered inline and\n"
     "crossline respectively and i1 and x1 are the highest"),
     NULL},
  {const_cast<char*>("indexcorners"),
     (getter)ZgyReader_getindexcorners, (setter)NULL,
     const_cast<char*>(
     "Corner coordinates in index / ordinal space, ordered the\n"
     "same way as in corners. The same information could be derived\n"
     "from just size."),
     NULL},
  {const_cast<char*>("annotcorners"),
     (getter)ZgyReader_getannotcorners, (setter)NULL,
     const_cast<char*>(
     "Corner coordinates in annotation space, ordered the same way\n"
     "as in corners. The same information could be derived from size,\n"
     "annotstart, and annotinc."),
     NULL},
  {const_cast<char*>("nlods"),
     (getter)ZgyReader_getnlods, (setter)NULL,
     const_cast<char*>("Number of level-of-detail layers, including lod 0 a.k.a. full resolution."),
     NULL},
  //{const_cast<char*>("dataid"),
  //   (getter)ZgyReader_getdataid, (setter)NULL,
  //   const_cast<char*>("GUID set on file creation."),
  //   NULL},
  {const_cast<char*>("verid"),
     (getter)ZgyReader_getverid, (setter)NULL,
     const_cast<char*>("GUID set each time the file is changed."),
     NULL},
  //{const_cast<char*>("previd"),
  //   (getter)ZgyReader_getprevid, (setter)NULL,
  //   const_cast<char*>("GUID before last change."),
  //   NULL},
  {const_cast<char*>("meta"),
     (getter)ZgyReader_getmeta, (setter)NULL,
     const_cast<char*>(
     "A dictionary of all the meta information, which can\n"
     "later be passed as **kwargs to the ZgyWriter constructor.\n"
     "Does not include \"nlods\" or guids because this is derived\n"
     "information which is not acceptable to pass to the constructor."),
     NULL},
  {const_cast<char*>("statistics"),
     (getter)ZgyReader_getstatistics, (setter)NULL,
     const_cast<char*>("Statistics for all samples"),
     NULL},
  {const_cast<char*>("histogram"),
     (getter)ZgyReader_gethistogram, (setter)NULL,
     const_cast<char*>("Histogram for all samples"),
     NULL},
  {const_cast<char*>("filestats"),
     (getter)ZgyReader_getfilestats, (setter)NULL,
     const_cast<char*>("Information for display purposes only."),
     NULL},
  {NULL}  /* Sentinel */
};

/**
 * Some of these methos are shared between all three classes
 * so the functions ought to be prefixed ZgyCommon_ not Zgyreader_
 */
static PyGetSetDef ZgyUtils_getseters[] = {
  {const_cast<char*>("numthreads"),
     (getter)ZgyReader_getnumthreads, (setter)ZgyReader_setnumthreads,
     const_cast<char*>("How many threads to use. Currently ignored. Use iocontext instead."),
     NULL},
  {const_cast<char*>("verbose"),
     (getter)ZgyReader_getverbose, (setter)ZgyReader_setverbose,
     const_cast<char*>("How many debug messages to output to std::cerr."),
     NULL},
  {NULL}  /* Sentinel */
};

static PyTypeObject ZgyReaderClassType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    FULL_MODULE_NAME ".ZgyReader", /* tp_name */
    sizeof(ZgyClass),    /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)ZgyReader_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "ZGY reader",              /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    ZgyReader_methods,         /* tp_methods */
    0,                         /* tp_members */
    ZgyCommon_getseters,       /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)ZgyReader_init,  /* tp_init */
    0,                         /* tp_alloc */
    ZgyReader_new,             /* tp_new */
};

static PyTypeObject ZgyWriterClassType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    FULL_MODULE_NAME ".ZgyWriter", /* tp_name */
    sizeof(ZgyClass),    /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)ZgyWriter_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "ZGY writer",              /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    ZgyWriter_methods,         /* tp_methods */
    0,                         /* tp_members */
    ZgyCommon_getseters,       /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)ZgyWriter_init,  /* tp_init */
    0,                         /* tp_alloc */
    ZgyWriter_new,             /* tp_new */
};

static PyTypeObject ZgyUtilsClassType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    FULL_MODULE_NAME ".ZgyUtils", /* tp_name */
    sizeof(ZgyClass),    /* tp_basicsize */
    0,                         /* tp_itemsize */
    (destructor)ZgyUtils_dealloc, /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /* tp_flags */
    "ZGY utils",               /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    ZgyUtils_methods,          /* tp_methods */
    0,                         /* tp_members */
    ZgyUtils_getseters,        /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)ZgyUtils_init,   /* tp_init */
    0,                         /* tp_alloc */
    ZgyUtils_new,              /* tp_new */
};

#if 0

// Alternative approach: Define free functions
// that creates and returns a new instance of
// the appropriate reader or writer class.
// This allows having distinct names for the
// create and clone operations which still
// return the same class instance.

static PyObject *
ZgyModule_reader(PyObject* self, PyObject* args, PyObject* kwargs)
{
  printf("ZGY: Create reader\n");

  // Call the class object in order to create a new instance.
  // If arguments are needed then build using Py_BuildValue
  // and remember to Py_DECREF the args when done.
  // Or, just pass on the arguments we were sent.
  PyObject *instance = PyObject_Call((PyObject *)&ZgyReaderClassType, args, kwargs);
  return instance;
}
#endif

/**
 * In OpenZGY this is a no-op, because there is no set/get current error.
 * All errors should be thrown as exceptions.
 */
static PyObject *
ZgyModule_setErrorHooks(PyObject* self, PyObject* args)
{
  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef ZgyModuleMethods[] = {
  {"setErrorHooks", (PyCFunction)ZgyModule_setErrorHooks, METH_VARARGS,
   "This is a no-op for backwards compatibility\n"},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

static const char *zgy_doc = "Python bindings for the public ZGY API.";

static struct PyModuleDef zgymodule = {
  PyModuleDef_HEAD_INIT,
  MODULE_NAME, /* m_name, sets __name__, should probably be unqualified. */
  zgy_doc, /* module documentation, may be NULL */
  -1,      /* size of per-interpreter state of the module,
              or -1 if the module keeps state in global variables. */
  ZgyModuleMethods
};

static bool
createEnums()
{
  // Create enums used in the API using the functional inteface
  // https://docs.python.org/3/library/enum.html#functional-api
  // The "module" assignment is to help unpickling; I hope.

  _enum_SampleDataType = NULL;
  _enum_UnitDimension  = NULL;
  _enum_DecimationType = NULL;
  _enum_FinalizeAction = NULL;
  _enum_LodMode = NULL;
  bool ok = true;

  PyObject* enum_module = PyImport_ImportModule("enum");
  if (!enum_module)
    return false;

  PyObject* enum_class = PyObject_GetAttrString(enum_module, "Enum");
  if (!enum_class) {
    Py_DECREF(enum_module);
    return false;
  }

  if (ok) {
    PyObject* tags = PyList_New(0);
    PyList_Append(tags, Py_BuildValue("(si)", "unknown", (int)SampleDataType::unknown));
    PyList_Append(tags, Py_BuildValue("(si)", "int8",    (int)SampleDataType::int8));
    PyList_Append(tags, Py_BuildValue("(si)", "int16",   (int)SampleDataType::int16));
    PyList_Append(tags, Py_BuildValue("(si)", "float",   (int)SampleDataType::float32));
    PyList_Append(tags, Py_BuildValue("(si)", "float32", (int)SampleDataType::float32));
    PyObject* args = Py_BuildValue("sO", "SampleDataType", tags);
    PyObject* kwargs = Py_BuildValue("{s:s}", "module", FULL_MODULE_NAME);
    _enum_SampleDataType = PyObject_Call(enum_class, args, kwargs);
    ok = ok && (_enum_SampleDataType != NULL);
    Py_DECREF(kwargs);
    Py_DECREF(args);
    for (Py_ssize_t ii=0; ii<PyList_Size(tags); ++ii)
      Py_DECREF(PyList_GetItem(tags, ii));
    Py_DECREF(tags);
  }

  if (ok) {
    PyObject* tags = PyList_New(0);
    PyList_Append(tags, Py_BuildValue("(si)", "unknown",  (int)UnitDimension::unknown));
    PyList_Append(tags, Py_BuildValue("(si)", "time",     (int)UnitDimension::time));
    PyList_Append(tags, Py_BuildValue("(si)", "length",   (int)UnitDimension::length));
    PyList_Append(tags, Py_BuildValue("(si)", "arcangle", (int)UnitDimension::arcangle));
    PyObject* args = Py_BuildValue("sO", "UnitDimension", tags);
    PyObject* kwargs = Py_BuildValue("{s:s}", "module", FULL_MODULE_NAME);
    _enum_UnitDimension = PyObject_Call(enum_class, args, kwargs);
    ok = ok && (_enum_UnitDimension != NULL);
    Py_DECREF(kwargs);
    Py_DECREF(args);
    for (Py_ssize_t ii=0; ii<PyList_Size(tags); ++ii)
      Py_DECREF(PyList_GetItem(tags, ii));
    Py_DECREF(tags);
  }

  if (ok) {
    PyObject* tags = PyList_New(0);
    PyList_Append(tags, Py_BuildValue("(si)", "LowPass",  (int)DecimationType::LowPass));
    PyList_Append(tags, Py_BuildValue("(si)", "WeightedAverage",  (int)DecimationType::WeightedAverage));
    PyList_Append(tags, Py_BuildValue("(si)", "Average",  (int)DecimationType::Average));
    PyList_Append(tags, Py_BuildValue("(si)", "Median",  (int)DecimationType::Median));
    PyList_Append(tags, Py_BuildValue("(si)", "Minimum",  (int)DecimationType::Minimum));
    PyList_Append(tags, Py_BuildValue("(si)", "Maximum",  (int)DecimationType::Maximum));
    //PyList_Append(tags, Py_BuildValue("(si)", "MinMax",  (int)DecimationType::MinMax));
    PyList_Append(tags, Py_BuildValue("(si)", "Decimate",  (int)DecimationType::Decimate));
    PyList_Append(tags, Py_BuildValue("(si)", "DecimateSkipNaN",  (int)DecimationType::DecimateSkipNaN));
    //PyList_Append(tags, Py_BuildValue("(si)", "DecimateRandom",  (int)DecimationType::DecimateRandom));
    PyList_Append(tags, Py_BuildValue("(si)", "AllZero",  (int)DecimationType::AllZero));
    //PyList_Append(tags, Py_BuildValue("(si)", "WhiteNoise",  (int)DecimationType::WhiteNoise));
    PyList_Append(tags, Py_BuildValue("(si)", "MostFrequent",  (int)DecimationType::MostFrequent));
    PyList_Append(tags, Py_BuildValue("(si)", "MostFrequentNon0",  (int)DecimationType::MostFrequentNon0));
    PyList_Append(tags, Py_BuildValue("(si)", "AverageNon0",  (int)DecimationType::AverageNon0));
    PyObject* args = Py_BuildValue("sO", "DecimationType", tags);
    PyObject* kwargs = Py_BuildValue("{s:s}", "module", FULL_MODULE_NAME);
    _enum_DecimationType = PyObject_Call(enum_class, args, kwargs);
    ok = ok && (_enum_DecimationType != NULL);
    Py_DECREF(kwargs);
    Py_DECREF(args);
    for (Py_ssize_t ii=0; ii<PyList_Size(tags); ++ii)
      Py_DECREF(PyList_GetItem(tags, ii));
    Py_DECREF(tags);
  }

  if (ok) {
    PyObject* tags = PyList_New(0);
    PyList_Append(tags, Py_BuildValue("(si)", "Delete",  (int)FinalizeAction::Delete));
    PyList_Append(tags, Py_BuildValue("(si)", "Keep",  (int)FinalizeAction::Keep));
    PyList_Append(tags, Py_BuildValue("(si)", "BuildIncremental",  (int)FinalizeAction::BuildIncremental));
    PyList_Append(tags, Py_BuildValue("(si)", "BuildFull",  (int)FinalizeAction::BuildFull));
    PyList_Append(tags, Py_BuildValue("(si)", "BuildNoHistogram",  (int)FinalizeAction::BuildNoHistogram));
    PyObject* args = Py_BuildValue("sO", "FinalizeAction", tags);
    PyObject* kwargs = Py_BuildValue("{s:s}", "module", FULL_MODULE_NAME);
    _enum_FinalizeAction = PyObject_Call(enum_class, args, kwargs);
    ok = ok && (_enum_FinalizeAction != NULL);
    Py_DECREF(kwargs);
    Py_DECREF(args);
    for (Py_ssize_t ii=0; ii<PyList_Size(tags); ++ii)
      Py_DECREF(PyList_GetItem(tags, ii));
    Py_DECREF(tags);
  }

  if (ok) {
    PyObject* tags = PyList_New(0);
    PyList_Append(tags, Py_BuildValue("(si)", "Default", (int)LodMode::Default));
    PyList_Append(tags, Py_BuildValue("(si)", "Early",   (int)LodMode::Early));
    PyList_Append(tags, Py_BuildValue("(si)", "Early1",  (int)LodMode::Early1));
    PyList_Append(tags, Py_BuildValue("(si)", "Never",   (int)LodMode::Never));
    PyList_Append(tags, Py_BuildValue("(si)", "Rebuild", (int)LodMode::Rebuild));
    PyObject* args = Py_BuildValue("sO", "LodMode", tags);
    PyObject* kwargs = Py_BuildValue("{s:s}", "module", FULL_MODULE_NAME);
    _enum_LodMode = PyObject_Call(enum_class, args, kwargs);
    ok = ok && (_enum_LodMode != NULL);
    Py_DECREF(kwargs);
    Py_DECREF(args);
    for (Py_ssize_t ii=0; ii<PyList_Size(tags); ++ii)
      Py_DECREF(PyList_GetItem(tags, ii));
    Py_DECREF(tags);
  }

  Py_DECREF(enum_class);
  Py_DECREF(enum_module);

  return ok;
}

static PyObject*
createNamedTuple(const char *name, const char *fields)
{
  PyObject* collections_module = PyImport_ImportModule("collections");
  if (!collections_module)
    return nullptr;

  PyObject* namedtuple_class = PyObject_GetAttrString(collections_module, "namedtuple");
  if (!namedtuple_class) {
    Py_DECREF(collections_module);
    return nullptr;
  }

  PyObject* args = Py_BuildValue("ss", name, fields);
#if PY_VERSION_HEX >= 0x03060000
  // The "module" assignment is to help unpickling; I hope.
  // It is only supported for Python >= 3.6
  PyObject* kwargs = Py_BuildValue("{s:s}", "module", FULL_MODULE_NAME);
#else
  PyObject* kwargs = nullptr;
#endif
  PyObject* result = PyObject_Call(namedtuple_class, args, kwargs);

  Py_XDECREF(kwargs);
  Py_DECREF(args);
  return result;
}

PyMODINIT_FUNC
PyInit_wrapper(void)
{
  // logger has not been initialized yet.
  // Too bad, because this output might have been useful.
  //self->pimpl_->logger_(1, "Module has been loaded\n");

  if (ZgyReaderClassType.tp_new == NULL)
    ZgyReaderClassType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&ZgyReaderClassType) < 0)
    return NULL;

  if (ZgyWriterClassType.tp_new == NULL)
    ZgyWriterClassType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&ZgyWriterClassType) < 0)
    return NULL;

  if (ZgyUtilsClassType.tp_new == NULL)
    ZgyUtilsClassType.tp_new = PyType_GenericNew;
  if (PyType_Ready(&ZgyUtilsClassType) < 0)
    return NULL;

  _zgy_error = PyErr_NewException(FULL_MODULE_NAME ".ZgyError", NULL, NULL);
  if (_zgy_error == NULL)
    return NULL;
  _zgy_error_FormatError = PyErr_NewException(FULL_MODULE_NAME ".ZgyFormatError", _zgy_error, NULL);
  if (_zgy_error_FormatError == NULL)
    return NULL;
  _zgy_error_NeedOldLibrary = PyErr_NewException(FULL_MODULE_NAME ".ZgyNeedOldLibrary", _zgy_error_FormatError, NULL);
  if (_zgy_error_FormatError == NULL)
    return NULL;
  _zgy_error_UpdateRules = PyErr_NewException(FULL_MODULE_NAME ".ZgyUpdateRules", _zgy_error_FormatError, NULL);
  if (_zgy_error_FormatError == NULL)
    return NULL;
  _zgy_error_CorruptedFile = PyErr_NewException(FULL_MODULE_NAME ".ZgyCorruptedFile", _zgy_error, NULL);
  if (_zgy_error_CorruptedFile == NULL)
    return NULL;
  _zgy_error_UserError = PyErr_NewException(FULL_MODULE_NAME ".ZgyUserError", _zgy_error, NULL);
  if (_zgy_error_UserError == NULL)
    return NULL;
  _zgy_error_InternalError = PyErr_NewException(FULL_MODULE_NAME ".ZgyInternalError", _zgy_error, NULL);
  if (_zgy_error_InternalError == NULL)
    return NULL;
  _zgy_error_EndOfFile = PyErr_NewException(FULL_MODULE_NAME ".ZgyEndOfFile", _zgy_error, NULL);
  if (_zgy_error_EndOfFile == NULL)
    return NULL;
  _zgy_error_SegmentIsClosed = PyErr_NewException(FULL_MODULE_NAME ".ZgySegmentIsClosed", _zgy_error, NULL);
  if (_zgy_error_SegmentIsClosed == NULL)
    return NULL;
  _zgy_error_Aborted = PyErr_NewException(FULL_MODULE_NAME ".ZgyAborted", _zgy_error, NULL);
  if (_zgy_error_Aborted == NULL)
    return NULL;
  _zgy_error_MissingFeature = PyErr_NewException(FULL_MODULE_NAME ".ZgyMissingFeature", _zgy_error, NULL);
  if (_zgy_error_MissingFeature == NULL)
    return NULL;
  _zgy_error_IoError = PyErr_NewException(FULL_MODULE_NAME ".ZgyIoError", _zgy_error, NULL);
  if (_zgy_error_IoError == NULL)
    return NULL;

  if (!createEnums())
    return NULL;

  _namedtuple_Statistics = createNamedTuple("SampleStatistics", "cnt sum ssq min max");
  if (_namedtuple_Statistics == nullptr)
    return nullptr;

  _namedtuple_Histogram = createNamedTuple("SampleHistogram", "cnt min max bin");
  if (_namedtuple_Histogram == nullptr)
    return nullptr;

  _namedtuple_FileStats = createNamedTuple("FileStats", "fileVersion fileSize headerSize dataStart alphaNormalCount alphaNormalSizePerEntry alphaCompressedCount alphaCcompressedSize alphaMissingCount alphaConstantCount brickNormalCount brickNormalSizePerEntry brickCompressedCount brickCompressedSize brickMissingCount brickConstantCount usedSize usedIfUncompressed compressionFactor isCompressed segmentSizes");
  if (_namedtuple_FileStats == nullptr)
    return nullptr;

  PyObject* m = PyModule_Create(&zgymodule);
  if (m == NULL)
    return NULL;

  // Add classes to the interpreter now that we know all went well.
  // AddObject will steal the reference we send, only on success.
  // I should be checking for errors (ret<0) but I won't.
  Py_INCREF(&ZgyReaderClassType);
  Py_INCREF(&ZgyWriterClassType);
  Py_INCREF(&ZgyUtilsClassType);
  Py_INCREF(_zgy_error); // Intentionally being leaked
  Py_INCREF(_zgy_error_FormatError);
  Py_INCREF(_zgy_error_CorruptedFile);
  Py_INCREF(_zgy_error_UserError);
  Py_INCREF(_zgy_error_InternalError);
  Py_INCREF(_zgy_error_EndOfFile);
  Py_INCREF(_zgy_error_SegmentIsClosed);
  Py_INCREF(_zgy_error_Aborted);
  Py_INCREF(_zgy_error_MissingFeature);
  Py_INCREF(_zgy_error_IoError);

  PyModule_AddObject(m, "ZgyReader", (PyObject *)&ZgyReaderClassType);
  PyModule_AddObject(m, "ZgyWriter", (PyObject *)&ZgyWriterClassType);
  PyModule_AddObject(m, "ZgyUtils",  (PyObject *)&ZgyUtilsClassType);
  PyModule_AddObject(m, "ZgyError",           _zgy_error);
  PyModule_AddObject(m, "ZgyFormatError",     _zgy_error_FormatError);
  PyModule_AddObject(m, "ZgyCorruptedFile",   _zgy_error_CorruptedFile);
  PyModule_AddObject(m, "ZgyUserError",       _zgy_error_UserError);
  PyModule_AddObject(m, "ZgyInternalError",   _zgy_error_InternalError);
  PyModule_AddObject(m, "ZgyEndOfFile",       _zgy_error_EndOfFile);
  PyModule_AddObject(m, "ZgySegmentIsClosed", _zgy_error_SegmentIsClosed);
  PyModule_AddObject(m, "ZgyAborted",         _zgy_error_Aborted);
  PyModule_AddObject(m, "ZgyMissingFeature",  _zgy_error_MissingFeature);
  PyModule_AddObject(m, "ZgyIoError",         _zgy_error_IoError);
#if 1
  // Exported just for testing. Normally the app only sees
  // objects created with this type, not the type itself.
  Py_INCREF(_namedtuple_Statistics);
  Py_INCREF(_namedtuple_Histogram);
  Py_INCREF(_namedtuple_FileStats);
  PyModule_AddObject(m, "SampleStatistics",     _namedtuple_Statistics);
  PyModule_AddObject(m, "SampleHistogram",      _namedtuple_Histogram);
  PyModule_AddObject(m, "FileStats",            _namedtuple_FileStats);
#endif

  Py_INCREF(_enum_SampleDataType);
  Py_INCREF(_enum_UnitDimension);
  Py_INCREF(_enum_DecimationType);
  Py_INCREF(_enum_FinalizeAction);
  Py_INCREF(_enum_LodMode);
  PyModule_AddObject(m, "SampleDataType", _enum_SampleDataType);
  PyModule_AddObject(m, "UnitDimension",  _enum_UnitDimension);
  PyModule_AddObject(m, "DecimationType", _enum_DecimationType);
  PyModule_AddObject(m, "FinalizeAction", _enum_FinalizeAction);
  PyModule_AddObject(m, "LodMode",        _enum_LodMode);

  return m;
}
