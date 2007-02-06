#!/bin/bash

path=`echo $0 | sed 's/testperfs.sh//'`
executab=`echo -e $path'crlibm_testperf'`

if [ "$1" = "" ]; then
    iter=100000
else 
    iter=$1
fi

$executab log RN $iter
$executab exp RN $iter
$executab sin RN $iter
$executab cos RN $iter
$executab tan RN $iter
$executab asin RN $iter
$executab acos RN $iter
$executab atan RN $iter
$executab log10 RN $iter
$executab log2 RN $iter
$executab sinpi RN $iter
$executab cospi RN $iter
$executab tanpi RN $iter

