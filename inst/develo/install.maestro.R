# module load libiconv
# Currently Loaded Modulefiles:
#   1) gcc/9.2.0(default)   4) freetype/2.10.1      7) fribidi/1.0.12  10) java/1.8.0  13) geos/3.8.0     16) libtiff/4.1.0                 19) glib/2.66.7    22) libffi/3.3     25) ImageMagick/6.9.10-91  28) udunits/2.2.26       
# 2) libxml2/2.9.10       5) fontconfig/2.13.91   8) curl/7.68.0     11) gdal/3.0.4  14) R/4.3.0        17) cairo/1.17.4                  20) pcre/8.43      23) xorg-libs/1.1  26) gsl/2.6                29) libjpeg-turbo/2.0.6  
# 3) libiconv/1.16        6) harfbuzz/2.6.7       9) cmake/3.19.4    12) proj/7.0.0  15) libpng/1.6.37  18) gobject-introspection/1.66.1  21) pixman/0.38.4  24) expat/2.4.1    27) java/13.0.2            
# 
# 
# Currently Loaded Modulefiles:
#   1) R/4.1.0              4) libiconv/1.16        7) harfbuzz/2.6.7  10) cmake/3.19.4  13) proj/7.0.0      16) gsl/2.6        19) libtiff/4.1.0  22) libffi/3.3     25) pcre/8.43                     28) ImageMagick/6.9.10-91  
# 2) gcc/9.2.0(default)   5) freetype/2.10.1      8) fribidi/1.0.12  11) java/1.8.0    14) geos/3.8.0      17) java/13.0.2    20) cairo/1.17.4   23) xorg-libs/1.1  26) pixman/0.38.4                 29) libjpeg-turbo/2.0.6    
# 3) libxml2/2.9.10       6) fontconfig/2.13.91   9) curl/7.68.0     12) gdal/3.0.4    15) udunits/2.2.26  18) libpng/1.6.37  21) glib/2.66.7    24) expat/2.4.1    27) gobject-introspection/1.66.1  


# R_install_packages --git-hub C3BI-pasteur-fr/UTechSCB-SCHNAPPs
# R_install_packages shinythemes
# R_install_packages ggalluvial 
# R_install_packages debugme
# R_install_packages kableExtra
# 
# R -e "withr::with_makevars(c(PKG_LIBS = '-liconv'), install.packages('haven',repos = 'http://cran.us.r-project.org'), assignment = '+=')"
# R -e "withr::with_makevars(c(PKG_LIBS = '-liconv'), install.packages('readxl',repos = 'http://cran.us.r-project.org'), assignment = '+=')"
# R -e "withr::with_makevars(c(PKG_LIBS = '-liconv'), install.packages('tidyverse',repos = 'http://cran.us.r-project.org'), assignment = '+=')"
# 
# R_install_packages --bioclite InteractiveComplexHeatmap


