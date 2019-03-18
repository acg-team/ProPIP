#!/usr/bin/env bash

var1=$(cat ./out_align_pip_likelihood.log);
var2=-101.389573720629116;
#status=$(echo $var1'=='$var2 | bc -l);
status=$(awk "function abs(value){return (value<0?-value:value);}; BEGIN {print (abs($var1)-abs($var2))}");

if [ "$status" -eq "0" ]; then
   exit 0;
else
   exit 1;
fi
