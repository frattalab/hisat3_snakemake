#!/bin/bash

input=$2
output=$3

awk  -F "\t" '{print $1 "\t" $2 "\t" $2 + 1 "\t" $5}' $input $output