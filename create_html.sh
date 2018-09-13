#!/bin/bash

TEMPDIR=$1
#cat template.html | sed "s@EVENTS_DATA@$(cat "$1"/times.txt)@g" | sed "s@GROUPS_DATA@$(cat "$1"/groups.txt)@g" > reditools.html
cat template.html | sed "/EVENTS_DATA/{s/EVENTS_DATA//g
r "$1"/times.txt
}" | sed "/GROUPS_DATA/{s/GROUPS_DATA//g
r "$1"/groups.txt
}" > reditools.html
