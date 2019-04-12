#! /usr/bin/env bash

set -e
#set -o pipefail
set -u

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(dirname $0)
readonly -a ARGS=("$@")
readonly NARGS="$#"

if [[ $NARGS -ne 1 ]]; then
	>&2 echo "Summarize RASE prediction statistics (extract last lines)."
	>&2 echo "usage: $PROGNAME directory"
	exit 1
fi

ls "$1"/*.predict.tsv > /dev/null 2> /dev/null \
	|| {
	echo "No RASE prediction files found (*.predict.tsv)."
	exit 1
}

dbs=$(ls "$1"/*.predict.tsv | perl -pe 's/.*__(.*)\.predict.tsv.*/\1/g' | sort | uniq)

for db in $dbs; do
	{
		cd "$1"
		head -n1 $(ls *__$db.predict.tsv | head -n1)
		tail -n1 *__$db.predict.tsv
	} \
		| perl -pe 's/( <==\n|__)/\t/g' \
		| perl -pe 's/(\.predict\.tsv|==> )//g' \
		| perl -pe 's/datetime/experiment\tdb\tdatetime/g' \
		| egrep -v ^$ \
		> $db.summary.tsv
done

