#!/bin/bash
a='5'
anguloB="`calc "${a}*atan(1.0)/45.0"`"
angulo="`echo ${anguloB} | sed 's/^~//g'`"
echo ${anguloB}
echo ${angulo}
