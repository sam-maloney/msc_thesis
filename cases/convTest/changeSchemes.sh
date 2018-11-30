#!/bin/bash
for file in ./*x/system/fvSchemes; do
   sed -i -e "s/limiter.*;/limiter        $1;/g" $file 
   sed -i -e "s/balanced.*;/balanced       $2;/g" $file 
done
