#!/bin/bash
cd 0
for name in base*
do
    newname="$(echo "$name" | cut -c6-)"
    cp "$name" "$newname"
done
