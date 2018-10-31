#!/bin/bash

if [ $# -ne 2 ]; then
  echo 1>&2 "Usage: $0 <db-path> <fasta-file>"
  exit 3
fi

set -e

if [ ! -d $1 ]
  then
	echo "database directory does not exist!"
	exit 1
  else
	db="$(dirname $(readlink -e $1))/$(basename $1)"
fi

if [ ! -f $2 ]
  then
	echo "fasta file does not exist or is not a file!"
	exit 1
  else
	data="$(dirname $(readlink -e $2))"
	genome="$(basename $2)"
fi

echo "db: $db"
echo "data: $data"
echo "genome: $genome"


sudo docker run \
	--rm \
	-v $db:/db \
	-v $data:/data \
	oschwengers/referenceseeker:latest \
        -v \
	$genome

