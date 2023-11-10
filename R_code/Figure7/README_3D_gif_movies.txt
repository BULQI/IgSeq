# ReadMe for generate 3D gif movies as in analysis_v1.6_MFA_3Dgifmovie_v1.0.R

## Dependent R packages
The R script depends on two packages among others, magick and rgl.
It is the rgl package need magick to render 3d gif generation. When you install the magick R package under conda environment, you might see troubles saying there is not library magick++ to use, even after you have installed libmagick++-dev (for debian linux or ubuntu). The trouble comes from that the conda env can not correct identify system-wise installed libmagick++. In this case, we need to maually provide the include and library directory to R by

    R CMD INSTALL --configure-vars='INCLUDE_DIR="PATH" LIB_DIR="PATH"'

You need figure out the path to include and library directory. To do that you need to use pkg-config

    pkg-config --print-variables Magick++  
    pkg-config --variable pcfiledir Magick++ 
    pkg-config --variable includedir Magick++ 

For detailed commands please check man pkg-config. 

Eventually the following command works for me

    R CMD INSTALL --configure-vars='INCLUDE_DIR="/usr/include/x86_64-linux-gnu//ImageMagick-6 -I/usr/include/ImageMagick-6" LIB_DIR="/usr/lib/x86_64-linux-gnu"' ~/Downloads/magick_2.8.0.tar.gz

Of course, you need to download the magick package 2.8.0 from R CRAN at https://cran.r-project.org/src/contrib/Archive/magick/ .


