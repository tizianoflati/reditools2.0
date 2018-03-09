#!/bin/bash

TEMPDIR=$1
cat template.html | sed "s@EVENTS_DATA@$(cat "$1"/times.txt)@g" | sed "s@GROUPS_DATA@$(cat "$1"/groups.txt)@g" > reditools.html
