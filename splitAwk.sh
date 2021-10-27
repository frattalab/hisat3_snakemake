#!/bin/bash

input=$2
prefix=$3
echo "i'm input"
echo $input
echo "i'm prefix"
echo $prefix
awk -v dirname=${prefix} 'NR>1 {filename=dirname$3".txt"; print $0 >> filename }' $input