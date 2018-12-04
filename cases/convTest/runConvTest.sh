#!/bin/bash
echo
echo Removing previous results...
rm -r "$1$2$3"
rm -r *x/1/ *x/0.5/ *x/sets/
echo Creating new folder...
mkdir "$1$2$3"
cp "ref_rho$3" "$1$2$3/ref_rho"
cp makeConvergencePlot.m "$1$2$3"
./changeSchemes.sh $1 $2
echo Beginning simulations:
for dir in *x/; do
    cd $dir
    N="${dir//x\/}"
    echo "    Running N = $N case..."
    ./computeFields 1 $N $3
    dbnsPotentialFoam > log
    sample > /dev/null
    cp sets/1/lineY_rho.xy "../$1$2$3/${N}_rho"
    cd ..
done
./fixRhoSamples.sh "$1$2$3"
echo Convergence Testing Finished!
echo
