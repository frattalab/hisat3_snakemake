#!/bin/bash

input=$1
output=$2

awk  -F "\t" '{print $1 "\t" $2 "\t" $2 + 1 "\t" $5}' $input > $output