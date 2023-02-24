// Copyright 2017-2020, Schlumberger
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

#pragma once

#include "declspec.h"
#include <stdexcept>

/** \file exception.h
 *
 * \brief Defines exceptions that may be raised by %OpenZGY.
 *
 * These classes are both visible to the OpenZGY public API and referenced
 * directly from the implementation classes. I apologize for the broken
 * encapsulation. Re-mapping the exceptions in the API layer didn't seem
 * worth the trouble.
 */

namespace OpenZGY { namespace Errors {
#if 0
}}
#endif

/**
 * \defgroup exceptions List of exceptions
 * @{
 */

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4275) // std::runtine_error not dll-exported
#endif

/**
 * \brief Base class for all exceptions thrown by %OpenZGY.
 * \details Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyError: public std::runtime_error
{
protected:
  /** \copybrief OpenZGY::Errors::ZgyError */
  ZgyError(const std::string& arg);
};

#ifdef _MSC_VER
#pragma warning(pop)
#endif

/**
 * \brief Corrupted or unsupported ZGY file.
 *
 * In some cases a corrupted file might lead to a ZgyInternalError
 * or ZgyEndOfFile being thrown instead of this one. Because it isn't
 * always easy to figure out the root cause.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyFormatError: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyFormatError */
  ZgyFormatError(const std::string& arg);
};

/**
 * \brief Please use the old ZGY access library.
 *
 * This file uses features that are only supported in the old ZGY
 * library. E.g. it might be an old style compressed file.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyNeedOldLibrary: public ZgyFormatError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyNeedOldLibrary */
  ZgyNeedOldLibrary(const std::string& arg);
};

/**
 * \brief Cannot be opened for update.
 *
 * The file is valid but some condition prevents it from being opened
 * for update. The application may need to delete and re-create the file.
 * Possible reasons are:
 *
 * - The file is not the latest version.
 * - The file is on the cloud but was created by the old accessor.
 * - The file is on the cloud but was copied there by sdutil.
 * - Possibly others.
 *
 * The reason that OpenZGY can even know that a cloud file is not a
 * file it created itself has to do with alignment and how the file
 * is split into segments.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyUpdateRules: public ZgyFormatError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyUpdateRules */
  ZgyUpdateRules(const std::string& arg);
};

/**
 * \brief The ZGY file became corrupted while writing to it.
 *
 * No further writes are allowed on this file because a previous write
 * raised an exception and we don't know the file's state. Subsequent
 * writes will also throw this exception.

 * The safe approach is to assume that the error caused the file to
 * become corrupted. It is recommended that the application closes and
 * deletes the file.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyCorruptedFile: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyCorruptedFile */
  ZgyCorruptedFile(const std::string& arg);
};

/**
 * \brief Exception that might be caused by the calling application.
 *
 * Determining whether a problem is the fault of the calling application
 * or the %OpenZGY library itself can be guesswork. Application code
 * might choose to treat ZgyUserError and ZgyInternalError the same way.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyUserError: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyUserError */
  ZgyUserError(const std::string& arg);
};

/**
 * \brief Exception that might be caused by a bug in %OpenZGY.
 *
 * Determining whether a problem is the fault of the calling application
 * or the %OpenZGY library itself can be guesswork. Application code
 * might choose to treat ZgyUserError and ZgyInternalError the same way.
 *
 * A corrupt file might also be reported as ZgyInternalError instead of
 * the more appropriate ZgyFormatError.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyInternalError: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyInternalError */
  ZgyInternalError(const std::string& arg);
};

/**
 * \brief Trying to read past EOF.
 *
 * This is always considered an error, and is often due to a corrupted
 * ZGY file. So this error should probably be treated as a ZgyFormatError.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyEndOfFile: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyEndOfFile */
  ZgyEndOfFile(const std::string& arg);
};

/**
 * \brief Exception used internally to request a retry.
 *
 * A write to the cloud failed because the region that was attempted
 * written had already been flushed. And the cloud back-end does not
 * allow writing it again. The calling code, still inside the OpenZGY
 * library, should be able to catch and recover from this problem.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgySegmentIsClosed: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgySegmentIsClosed */
  ZgySegmentIsClosed(const std::string& arg);
};

/**
 * \brief User aborted the computation.
 *
 * If the user supplied a progress callback and this callback returned
 * false then the operation in progress will and by throwing this
 * exception. Which means that this is not an error; it is a consequence
 * of the abort.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyAborted: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyAborted */
  ZgyAborted(const std::string& arg);
};

/**
 * \brief Missing feature.
 *
 * Raised if some optional plug-in (e.g. some cloud back end or a
 * compressor) was loaded or explicitly requested, so we know about
 * it, but the plug-in is not operational for some reason.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyMissingFeature: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyMissingFeature */
  ZgyMissingFeature(const std::string& arg);
};

/**
 * \brief Exception from the I/O layer.
 *
 * Some error was received from a linux syscall acting on a file.
 * For cloud access the actual exception thrown by the back end
 * might be reported instead of wrapping it into a ZgyIoError.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyIoError: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyIoError */
  ZgyIoError(const std::string& filename, int system_errno);
};

/**
 * \brief Exception from the Windows I/O layer.
 *
 * Some error was received from a Windows syscall acting on a file.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyWindowsError : public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyWindowsError */
  ZgyWindowsError(const std::string& filename, unsigned long windows_errno);
};

/**
 * \brief Data must be read-only.
 *
 * Some operations such as alturl require that the data is flagged
 * read-only. This was not the case here.
 *
 * Thread safety: Exceptions defined in this module are safe.
 */
class OPENZGY_API ZgyNotReadOnlyError: public ZgyError
{
public:
  /** \copybrief OpenZGY::Errors::ZgyNotReadOnlyError */
  ZgyNotReadOnlyError(const std::string& filename);
};

/**
 * @}
 */

}} // namespace
