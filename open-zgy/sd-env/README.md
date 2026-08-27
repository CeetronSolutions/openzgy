## Contents of sd-env

Make targets for unattended builds. This is the normal usage.

See [Building and Testing](../README.md#building-on-linux)
in the parent [README.md](../README.md]).

Output is stored as ../sdapi-pkgs/DIST-sdapi_linux64.tar.gz
It needs to be copied to ../pkg/sdapi_linux64.tar.gz before
building OpenZGY for DIST. The files in ../pkg are not named
after the distro, so you can only have one set in that folder.

    make sdapi      -- Build SDAPI only
    make sdapi-DIST -- Build SDAPI only, for DIST only

Make targets for manual builds, typically run in this order:

    make clobber    -- Removes docker instances but not images
    make manual     -- Rebuilds images used for manual builds
                       Source mapped from $(WORKAREA)
    make run        -- Create instances but don't start them
    make xterm      -- Start instances in individual windows
    ./build.sh      -- to be executed inside each xterm

To manually build one distro only, use e.g. manual-focal and
run-focal. Then start with "docker start -ai sd-focal"
instead of "make xterm".

Note that manual builds expect to find the SDAPI source at
../seismic-store-cpp-lib
in the host. The source will be mapped as a volume into the
container. Any patching needs to be done outside. "make sdapi"
will instead copy the source and will handle patching itself.

<!--
NOTE: The help text in Makefile should be kept in sync with README.md
-->
