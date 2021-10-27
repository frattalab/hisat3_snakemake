#!/bin/bash

input=$2
prefix=$3

awk -F "\t" "{print > $prefix $3 ".txt"}" $input