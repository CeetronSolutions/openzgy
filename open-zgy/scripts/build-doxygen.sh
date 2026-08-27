#!/bin/bash -e

set -x

# Build both the internal and the public Doxygen repoprts.
# Build both single page pdf and html with hyperlinks.
# The latter is skipped if latex is not installed.

BUILDROOT=$(dirname $(realpath ${0:?}))/..
DOX_DIR=${BUILDROOT}/build/deploy/native
DOX_TMP=${BUILDROOT}/build/temp/native/doxygen
SOURCE=${BUILDROOT}/native/src
LATEX=/usr/bin/latex
PDFTK=/usr/bin/pdftk
PDFTK=
HAS_LATEX=NO
if [ -n "${LATEX}" -a -x "${LATEX}" ]
then
    HAS_LATEX=YES
fi

cd ${BUILDROOT}/native/src

# Public documentation

if [ "$1" = "" -o "$1" = "apidoc" ]
then
    /bin/rm -rf ${DOX_TMP}/apidoc
    /bin/mkdir -p ${DOX_TMP}/apidoc ${DOX_DIR}
    (cat Doxyfile; \
     echo GENERATE_LATEX=${HAS_LATEX}; \
     echo ENABLED_SECTIONS=SSTORE; \
     echo 'PROJECT_NAME = "OpenZGY/C++ Public API"'; \
     echo OUTPUT_DIRECTORY=${DOX_TMP}/apidoc) | doxygen -
    /bin/tar -C ${DOX_TMP} -zcf ${DOX_DIR}/apidoc.tgz apidoc
    test -r ${DOX_TMP}/apidoc/html/index.html
    ls -l ${DOX_DIR}
    if [ -n "${LATEX}" -a -x "${LATEX}" ]
    then
        make -C ${DOX_TMP}/apidoc/latex > ${DOX_TMP}/apidoc/latexlog.txt 2>&1 < /dev/null
	      /bin/mv ${DOX_TMP}/apidoc/latex/refman.pdf ${DOX_DIR}/apidoc.pdf
        if [ -n "${PDFTK}" -a -x "${PDFTK}" ]
        then
	          ${PDFTK} ${DOX_DIR}/apidoc.pdf \
                     background ../../doc/images/draft-bg.pdf \
                     output ${DOX_DIR}/apidoc.pdf.tmp && \
                /bin/mv ${DOX_DIR}/apidoc.pdf.tmp ${DOX_DIR}/apidoc.pdf
        fi
    fi
    echo "Main output: ${DOX_TMP}/apidoc/html/index.html"
    # Some installers may want the unzipped version.
	  #/bin/rm -rf ${DOX_TMP}/apidoc
fi

# Internal documentation

if [ "$1" = "" -o "$1" = "intdoc" ]
then
    /bin/rm -rf ${DOX_TMP}/intdoc
    /bin/mkdir -p ${DOX_TMP}/intdoc ${DOX_DIR}
    (cat Doxyfile; \
     echo GENERATE_LATEX=${HAS_LATEX}; \
     echo ENABLED_SECTIONS=IMPL SSTORE; \
     echo WARNINGS=NO; \
     echo WARN_IF_UNDOCUMENTED=NO; \
     echo EXTRACT_ALL=YES; \
     echo EXTRACT_PRIVATE=YES; \
     echo EXTRACT_STATIC=YES; \
     echo "# For caller and callee graphs"; \
     echo EXTRACT_ALL            = YES; \
     echo EXTRACT_PRIVATE        = YES; \
     echo EXTRACT_PACKAGE        = YES; \
     echo EXTRACT_STATIC         = YES; \
     echo EXTRACT_ANON_NSPACES   = YES; \
     echo HAVE_DOT               = YES; \
     echo DOT_IMAGE_FORMAT       = svg; \
     echo DOT_CLEANUP            = YES; \
     echo CALL_GRAPH             = NO; \
     echo CALLER_GRAPH           = NO; \
     echo INPUT = . impl ../../doc; \
     echo 'PROJECT_NAME = "OpenZGY/C++ API and Internals"'; \
     echo OUTPUT_DIRECTORY=${DOX_TMP}/intdoc) | doxygen -
    /bin/tar -C ${DOX_TMP} -zcf ${DOX_DIR}/intdoc.tgz intdoc
    test -r ${DOX_TMP}/intdoc/html/index.html
    if [ -n "${LATEX}" -a -x "${LATEX}" ]
    then
        make -C ${DOX_TMP}/intdoc/latex > ${DOX_TMP}/intdoc/latexlog.txt 2>&1 < /dev/null
	      /bin/mv ${DOX_TMP}/intdoc/latex/refman.pdf ${DOX_DIR}/intdoc.pdf
        if [ -n "${PDFTK}" -a -x "${PDFTK}" ]
        then
	          ${PDFTK} ${DOX_DIR}/intdoc.pdf \
                     background ../../doc/images/draft-bg.pdf \
                     output ${DOX_DIR}/intdoc.pdf.tmp && \
                /bin/mv ${DOX_DIR}/intdoc.pdf.tmp ${DOX_DIR}/intdoc.pdf
        fi
    fi
    echo "Main output: ${DOX_TMP}/intdoc/html/index.html"
    # Some installers may want the unzipped version.
	  #/bin/rm -rf ${DOX_TMP}/intdoc
fi
