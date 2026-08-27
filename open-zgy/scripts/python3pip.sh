#!/bin/bash

# Run pip3 with up to 3 retries, with longer and longer timeouts.
# This is an attempt to fix problems in the gitlab build. Note that
# all failures are retried, not just timeouts. Because I am lazy. Also
# time each execution. Note that there are two versions of this
# script, using "pip3" and "python3 -m pip". Technically they should
# be interchangeable (with the second one preferred) but there are
# subtle differences. If it works, don't fix it. those. So, when
# changing to use this wrapper, choose the version that used to work

TIMEFORMAT='python3pip.sh elapsed: %R seconds - user %U - system %S - load %P'
time python3 -m pip "$@" ||
    (echo "RETRY#1"; time python3 -m pip --timeout=60 --default-timeout=60 "$@") ||
    (echo "RETRY#2"; time python3 -m pip --timeout=60 --default-timeout=60 "$@") ||
    (echo "RETRY#3"; time python3 -m pip --timeout=300 --default-timeout=300 "$@") ||
    (echo python3pip.sh FAILED "$@"; exit 1)
