#! /usr/bin/env bash

set -e
set -o pipefail

# $1 - file
# $2 - expected error code
test_script () {
    file="$1"
    experr="$2"
    echo "Testing $1"
    tmp=$(mktemp)
    ((./$file 2>&1)> "$tmp" || (code=$? && exit $(($code - $experr))) ) || ( echo "...failed ($code)" && cat "$tmp"  && exit 1)
}

DIR=`dirname $0`
export -f test_script

# argparse should return exit status 2
ls "$DIR"/*.py | parallel test_script {} 2
ls "$DIR"/*.R | parallel test_script {} 1
