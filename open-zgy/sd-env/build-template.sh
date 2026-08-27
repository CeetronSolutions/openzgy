#!/bin/bash -e

# Note: To build the SAuth aware version of SDAPI the osdu source code
# needs to be patched with 7 files (4 new, 3 modified) from the old
# seismic-store-client-api-cpp repository:
#
#    src/src/core/SDManager.cc
#    src/src/core/SDManager.h
#    src/src/core/SDManagerImpl.h
#    src/src/lib/auth/slb_auth_provider.cc
#    src/src/lib/auth/slb_auth_provider.h
#    src/src/lib/auth/slb_service_auth_provider.cc
#    src/src/lib/auth/slb_service_auth_provider.h
#
# Note: Submodules are required as of March 2019.
# remember to do on the host:
#   git pull
#   git submodule update --recursive --remote
# Also, cmake now needs to be version 3.1 or higher.
# Also, on debian / ubuntu there is a bug in the rules
# causing crc32c and sdapi to disagreee about where
# the libraries should be installed.

CMAKE=cmake3; type 2>/dev/null $CMAKE || CMAKE=cmake

SRC=../sdapi/src

cd ${HOME}
rm -rf build sdapi_linux64_osdu.tar.gz realname.txt; mkdir build; cd build

AZUREMSG="(with Google support only)"
AZUREDEF=
if [ -r /opt/vcpkg/scripts/buildsystems/vcpkg.cmake ]
then
    AZUREMSG="(with Google and Azure support)"
    AZUREDEF="-DCMAKE_TOOLCHAIN_FILE=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake -DOPTIONAL_STORAGE_PROVIDERS_ENABLED:STRING=azure-curl"
elif [ -r /usr/local/lib/libazurestorage.so ]
then
    AZUREMSG="(with Google and (dynamic?) Azure support)"
    AZUREDEF="-DOPTIONAL_STORAGE_PROVIDERS_ENABLED:STRING=azure-curl"
elif [ -r /usr/local/lib/libazurestorage.a ]
then
    AZUREMSG="(with Google and (static?) Azure support)"
    AZUREDEF="-DOPTIONAL_STORAGE_PROVIDERS_ENABLED:STRING=azure-curl"
fi
echo "Starting build of ${LINUXDISTRO:-unknown} ${1:-Release} configuration ${AZUREMSG}"

#Version on commandline?
#VERSION=$(echo "${2}" | tr . ' ')
#VERSIONDEF=$(bash -c "set -- $VERSION;"' echo -DLIB_VERSION_MAJOR=${1:0} -DLIB_VERSION_MINOR=${2:0} -DLIB_VERSION_PATCH=${3:0}')

#Version from source? Patch will likely be zero.
#Because patch is expected to be set on the commandline.
fv=${SRC}/src/core/SDVersion.h
v_major=$(sed -n -e 's;//.*;;' -e 's;^.*Major *= *\([0-9][0-9]*\).*$;\1;p' $fv)
v_minor=$(sed -n -e 's;//.*;;' -e 's;^.*Minor *= *\([0-9][0-9]*\).*$;\1;p' $fv)
v_patch=$(sed -n -e 's;//.*;;' -e 's;^.*Patch *= *\([0-9][0-9]*\).*$;\1;p' $fv)
VERSIONDEF="-DLIB_VERSION_MAJOR=${v_major:-0} -DLIB_VERSION_MINOR=${v_minor:-0} -DLIB_VERSION_PATCH=${v_patch:-0}"
v_version="${v_major:-0}.${v_minor:-0}.${v_patch:-0}"
echo "Version: $v_version"

# Make sure we link with the default curl i.e. whatever is installed
# in the default location. If you need to build curl from sources
# (CentOS 7) then that custom build should simply be installed
# overwriting the original one. When running inside a docker
# container this is both simplar and safer. Note:
# It isn't enough to leave CURL_ENV unset; must be explicit.
export CURL_USE_SYSTEM_CURL=YES
$CMAKE -DCMAKE_BUILD_TYPE=${1:-Release} ${VERSIONDEF} ${AZUREDEF} ${SRC}

# This kludge is hopefully not needed any longer. On some platforms
# there was confusion about where the binaries ended up.
#make || true
#test -x crc32c/lib64 || ln -s lib crc32c/lib64
# This is new... centos7 builds, should't hurt others I think.
test -e /usr/local/include/crc32c || sudo ln -s  /home/me/sdapi/src/third-party/crc32c/include/crc32c /usr/local/include/

# This kludge is needed as of June 2021 but has to be applied outside
# the build container because inside the source is readonly.
# After some experimenting I found that this approach is the minimal
# change that makes all 6 config build cleanly when using vcpkg.
# I have not tested the same configs without vcpkg.
#sed -ie '/find_package.*[Ll]ib[Xx]ml2.*/s/\(find.*\)/set( LibXml2_FOUND true ) # \1/' src/CMakeLists.txt

# This kludge is for integrating the slb-only extensions to sdapi.
# As above, it needs to be run outside the container. The assumption
# is that the private repository folder is the parent of the osdu one.
# This REALLY ought to be done as a plug-in.
#cp -t src/src/core ../src/src/core/*
#cp -t src/src/lib/auth ../src/src/lib/auth/*

make
rm -rf tmp
mkdir -p tmp/include tmp/lib/linux64

function gethead {
    gitdir="${1:-.}/.git"
    if [ -r "${gitdir}/HEAD" ]
    then
	read hash ref < "${gitdir}/HEAD"
	if [ "x$hash" = "xref:" -a "x$ref" != "x" -a -r "${gitdir}/$ref" ]
	then
	    read hash < "${gitdir}/$ref"
	fi
    else
	hash=no-hash-available
    fi
    echo $hash
}

# Capture the hash of the OSDU part of SDAPI, and the version of some
# packages installed by vcpkg.
echo "LINUX: ${LINUXDISTRO:-unknown}" > tmp/version.txt
echo "SDAPI: $(gethead /home/me/sdapi)" >> tmp/version.txt
if [ -x /opt/vcpkg/vcpkg ]
then
    echo "VCPKG: $(gethead /opt/vcpkg)" >> tmp/version.txt
    echo "" >> tmp/version.txt
    echo "AZURE:" >> tmp/version.txt
    sudo date > /dev/null 2>&1 # Skip the first-time lecture.
    sudo /opt/vcpkg/vcpkg list >> tmp/version.txt || true
fi
echo "" >> tmp/version.txt
echo "DYNAMIC libsdapi.so" >> tmp/version.txt
ldd libsdapi.so >> tmp/version.txt || true

cp -a -t tmp/include ${SRC}/src/core/*.h
cp -a -t tmp/include ${SRC}/src/core/*.hh 2>/dev/null || true
cp -a -t tmp/lib/linux64 libsdapi*
if [ -r /opt/vcpkg/installed/x64-linux/lib ]
then
    (cd /opt/vcpkg/installed/x64-linux/lib; \
	   LIBS=$(find . -name \*.so\* ! -name libboost\* -print); \
     test -z "${LIBS}" || \
         /bin/cp -a -t ${HOME}/build/tmp/lib/linux64 ${LIBS})
    # Note that ${LIBS} might be empty if building all static.
    # The above copies both the actual shared objecs and the symlimks.
    # Not all symlinks are needed but it is tricky to figure out which.
    # As long as they remain symlinks it is not a big deal, but if
    # packaged somewhere that don't support symlinks (Python wheels
    # and most Microsoft stuff) this can lean to massive bloat.
    # All real files i.e. not symlinks need to be kept.
    # Files that will be dynamically loaded are need and tricky to detect.
    # To get the actual list of required files, except dynamically loaded:
    # readelf -d *.so | sed -ne '/SONAME/s;^.*[[]\(.*\)[]].*$;\1;p'
    # Or consider this heuristic: symlinks nameed *.so can in general be
    # removed because if there is any version handling at all it would
    # be pointless to have the library's SONAME set to the unadorned .so.
    # For files with multiple links this only saves the space of one of them.
fi
for file in $(find tmp/lib/linux64 -type f -name \*.so\* -print)
do
    old_rpath=$(patchelf --print-rpath "$file")
    patchelf --set-rpath '$ORIGIN' "$file"
    new_rpath=$(patchelf --print-rpath "$file")
    echo "Changed rpath of $file from ${old_rpath} to ${new_rpath}"
    readelf -d "$file" | grep PATH || true
done

tar zcf ../sdapi_linux64.tar.gz -C tmp version.txt include lib

# Embed Linux distro and compiler version in the file name.
# (Python version gets embedded also, but only for the wheels).
#
# TODO-Low fix this; Using the Salmon scheme which initially
# assumed there would be just one Linux distro.
#
#   1) Use just major.minor, or maybe just major, compiler version.
#      Bumping the patch number is very unlikely to break anything
#      and we normally end up ignoring it anyway.
#
#   2) So far I have been lucky that no two linux distros have chosen
#      the exact same gcc major.minor.patch. So the distro is sort of
#      implied by the compiler and not the other way around. Sorry.
#      Need to include the distro name like OpenZGY does.

set +e
GCCVERSION=$((g++ -dumpfullversion 2>/dev/null || g++ -dumpversion) | sed -e 's;\.[0-9][0-9][0-9]*;.9;g' -e 's;\.;;g')
NEWGCCVERS=$(echo ${GCCVERSION} | sed -e 's;[0-9]$;;')
PLATFORM=Lin64_gcc$(echo "$GCCVERSION" | sed -e 's;\.[0-9][0-9][0-9]*;.9;g' -e 's;\.;;g')
test ${PLATFORM} != Lin64_gcc492 || PLATFORM=Lin64_gcc485
NEWPLATFORM=${LINUXDISTRO:-unknown}-gcc${NEWGCCVERS}

tarball=sdapi_linux64_local.tar.gz
if [ ! -r ${HOME}/sdapi/src/src/lib/auth/slb_auth_provider.h ]
then
    # This looks like the open source SDAPI without any SLB extensions.
    # Nope... SAuth is now in the open source part.
    #tarball=sdapi_linux64_osdu.tar.gz
    echo "Info: No patching of SLB extensions is done, or needed."
else
    echo "Info: slb_auth_provider.h is present."
fi

# If the SDK is to be used e.g. in ZGY-Cloud it may need to be renamed and
# copied out of the container. The following might help.
echo ${PLATFORM}/${tarball} > ${HOME}/realname.txt
echo ${NEWPLATFORM}/${tarball} > ${HOME}/realname2.txt

echo "Successful build of ${LINUXDISTRO:-unknown} ${1:-Release} configuration ${AZUREMSG}"
echo "The following only applies if you used 'make manual'."
echo "Now run this command outside the docker container:"
echo
echo "docker cp ${HOSTNAME}:/home/me/sdapi_linux64.tar.gz ${PLATFORM}/${tarball}"
