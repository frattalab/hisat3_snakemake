#!/bin/bash

input=$1
output=$2

awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' $input > $output   
