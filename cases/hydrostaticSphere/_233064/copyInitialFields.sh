#!/bin/bash
cd 0
gunzip *
rm *0
for field in p U T rho; do
  cp "${field}" "${field}0"
  sed -i -e "s/object.*;/object      ${field}0;/g" "${field}0"
done
gzip *
cd ..
echo [0-9]* | xargs -n 1 cp 0/*0.gz
