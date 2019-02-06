#! /usr/bin/env bash

set -e
set -o pipefail
set -u

tmp="$(mktemp)"

cat > "$tmp"

head -n 1 "$tmp"  && tail -n +2 "$tmp" | sort $@


