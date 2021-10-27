#!/bin/bash

input=$2
prefix=$3

awk -v dirname=${prefix} 'NR>1 {filename=dirname$3".txt"; print $0 >> filename }' $input