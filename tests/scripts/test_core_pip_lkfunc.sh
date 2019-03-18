#!/usr/bin/env bash

var1=$(cat ./out_core_pip_lkfunc.estimates.json | grep LogLikelihood | awk -F":" '{ print $2}' | sed 's/\"//g' | sed 's/ //g');
var2=-21.358154531926616;
#status=$(echo $var1'=='$var2 | bc -l);
status=$(awk "function abs(value){return (value<0?-value:value);}; BEGIN {print (abs($var1)-abs($var2))}");

if [ "$status" -eq "0" ]; then
   exit 0;
else
   exit 1;
fi
