#!/bin/bash

echo "Started!"

cd ./EFz
echo "I am in $(pwd)"
ncl EFz_Ave_VerticalProfile.ncl
matlab -nodisplay -nojvm < plot_timeseries.m
cd ..

cd ./HFz
o "I am in $(pwd)"
ncl HFz_Ave_VerticalProfile.ncl
matlab -nodisplay -nojvm < plot_timeseries.m
cd ..

cd ./rainfall
echo "I am in $(pwd)"
grads -cblx step1.rainfall.run
matlab -nodisplay -nojvm < plot_timeseries.m
cd ..


echo "Done!"
