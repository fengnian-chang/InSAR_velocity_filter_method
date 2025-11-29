source_ve=ve.grd
res=0.01
width=$1
ve_filtered=ve_mf${width}.grd
R=-R73/112E/18N/46N
inc=0.01

T=0   # 0-1
vn=vn_T${T}

gmt grdfilter $source_ve -G$ve_filtered -Fm${width} -D1 -I${res} -V

#awk '{print $1,$2,$4}' geotiffs/velfit.dat | gmt surface -G${vn}.grd -I$inc -T$T $R -V
#gdal_translate -of GTiff ${vn}.grd ${vn}.tif
#gdalwarp2match.py ${vn}.tif decomposed/ve.tif Vn_matched.tif
#gdal_translate -of NetCDF Vn_matched.tif ${vn}.grd

gmt grdmath -M $ve_filtered DDX 1000000 MUL = dvedx${width}.grd
gmt grdmath -M $ve_filtered DDY 1000000 MUL = dvedy${width}.grd
gmt grdmath -M ${vn}.grd DDX 1000000 MUL = dvndx${width}.grd
gmt grdmath -M ${vn}.grd DDY 1000000 MUL = dvndy${width}.grd
gmt grdmath dvedy${width}.grd dvndx${width}.grd ADD 2 DIV  = dxy${width}.grd
gmt grdmath dxy${width}.grd 2 POW 2 MUL dvndy${width}.grd 2 POW ADD dvedx${width}.grd 2 POW ADD SQRT = EII${width}_T${T}.grd
gmt grdmath dxy${width}.grd 2 POW dvedx${width}.grd dvndy${width}.grd SUB 2 POW 4 DIV ADD SQRT = Eshear${width}_T${T}.grd
gmt grdmath dvndy${width}.grd dvedx${width}.grd ADD = Edil${width}_T${T}.grd
gmt grdmath dvndx${width}.grd dvedy${width}.grd SUB = Evort${width}_T${T}.grd
