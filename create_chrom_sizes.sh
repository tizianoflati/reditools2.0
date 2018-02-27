#!/bin/bash

#cat $1 | grep "^@SQ"| cut -f 2-3 | sed 's/SN://;s/LN://g'
cat $1 | cut -f 1-2
