#!/bin/bash

echo "Started!"

cd ./EFx
echo "I am in $(pwd)"
ncl EFx_Ave_VerticalProfile.ncl
matlab -nodisplay -nojvm < plot_mean.m
cd ..

cd ./EFz
o "I am in $(pwd)"
ncl EFz_Ave_VerticalProfile.ncl
matlab -nodisplay -nojvm < plot_mean.m
cd ..

cd ./MFx
echo "I am in $(pwd)"
ncl MFx_Ave_VerticalProfile.ncl
matlab -nodisplay -nojvm < plot_mean.m
cd ..

cd ./U
echo "I am in $(pwd)"
ncl U_Ave_VerticalProfile.ncl
matlab -nodisplay -nojvm < plot_mean.m
cd ..

echo "Done!"
