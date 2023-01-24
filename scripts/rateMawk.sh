#!/bin/bash

input=$1
output=$2
mawk -F '\t' '{sum6+=$6; sum8+=$8}END{print sum6,sum8}' $input > $output