#!/bin/bash

input=$1
output=$2
awk  -F "\t" '{if($7+$5 ==0) $4 = 0; else $4 = ($5)/($7+$5)} {print $1 "\t" $2 "\t" $2 + 1 "\t" $4}' $input > $output