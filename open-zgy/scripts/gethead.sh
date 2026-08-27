#!/bin/bash

# Crude replacement for "git rev-parse HEAD" to be used if git is not
# installed or if most of the .git folder has been excluded by .dockerignore.
# Beware that this might give the wrong result in some circumstances.

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

gethead "$1"
