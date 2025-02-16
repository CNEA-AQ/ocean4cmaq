#!/bin/bash
#
# Create OPEN and SURF files for CMAQ from 
# a land & coastline file
#
#---------------------------------------
#INPUTs:
start_date="2024-06-11"
LAND=data/land_10m_global.gpkg
#COAST=original/coastline.gpkg
DMS_file=data/DMS_monthly_global_2011.nc
CHL_file=data/CHLO_monthly_global_2019.nc
GRIDDESC_file="GRIDDESC"
GRIDNAME="US_SD"
#--------------------------------------
crs_latlon="EPSG:4326"
YYYY=${start_date:0:4}
  MM=${start_date:5:2}
  DD=${start_date:9:2}

#Read GRIDDESC
read projName xorig yorig dx dy nx ny nz <<< $(sed -n "/${GRIDNAME}/{n;p;}" $GRIDDESC_file)
read projId alp bet gam clon clat <<< $(sed -n "/${projName}/{n;p;}" $GRIDDESC_file)

if [[ $projId == "2" ]]  # (lambert conformal conic)
then
   crsOut="+proj=lcc +lat_1=${alp} +lat_2=${bet} +lon_0=${gam} +lat_0=${clat} +a=6370000.0 +b=6370000.0 +units=m" 
   read xc yc ellipsoidh <<<$( gdaltransform -s_srs "${crs_latlon}" -t_srs "${crsOut}" <<< $( echo "${clon} ${clat}" ) )
   xmin=$(bc -l <<<"$xorig")
   xmax=$(bc -l <<<"($xorig+ $dx*$nx)")
   ymin=$(bc -l <<<"$yorig")
   ymax=$(bc -l <<<"($yorig+ $dy*$ny)")

   read lon[1] lat[1] ell <<<$( gdaltransform -t_srs "${crs_latlon}" -s_srs "${crsOut}" <<< $( echo "${xmin} ${ymin}" ) )
   read lon[2] lat[2] ell <<<$( gdaltransform -t_srs "${crs_latlon}" -s_srs "${crsOut}" <<< $( echo "${xmax} ${ymax}" ) )
   read lon[3] lat[3] ell <<<$( gdaltransform -t_srs "${crs_latlon}" -s_srs "${crsOut}" <<< $( echo "${xmin} ${ymax}" ) )
   read lon[4] lat[4] ell <<<$( gdaltransform -t_srs "${crs_latlon}" -s_srs "${crsOut}" <<< $( echo "${xmax} ${ymin}" ) )
   
   #Find boundaries in latlon coordinates:
   IFS=$'\n'
   minlon=$(echo "${lon[*]}" | sort -nr | tail -n1)
   maxlon=$(echo "${lon[*]}" | sort -nr | head -n1)
   minlat=$(echo "${lat[*]}" | sort -nr | tail -n1)
   maxlat=$(echo "${lat[*]}" | sort -nr | head -n1)

else
   echo "Projection: $projId, not suported."
fi

echo "----"
echo "CHECKS:"
echo "proj: $crsOut"
echo "xc, yc, xmin ymin xmax ymax: $xc, $yc, $xmin,$ymin, $xmax, $ymax"
echo "nx ny: $nx, $ny"
echo "center min max (lon-lat): $clon $clat. $minlon $minlat. $maxlon $maxlat"
echo "----"
################################################################################
##---
## Remap the CHLO and DMS fields to local grid:
echo "Remap DMS.."
gdal_translate -of GTiff -b 1 NETCDF:"$DMS_file":"dms$MM" tmp.tiff
gdalwarp -overwrite -srcnodata 99999 -s_srs "$crs_latlon" -te ${minlon} ${minlat} ${maxlon} ${maxlat} tmp.tiff tmp1.tiff
gdalwarp -overwrite -s_srs "$crs_latlon" -t_srs "$crsOut" tmp1.tiff tmp2.tiff
gdalwarp -overwrite -srcnodata -999.99 -t_srs "$crsOut" -tr $dx $dy -te $xmin $ymin $xmax $ymax -r bilinear -ot Float32 -of Gtiff tmp2.tiff tmp3.tiff
gdal_translate -of NETCDF tmp3.tiff DMS.nc
rm tmp*

echo "Remap CHLO.."
gdal_translate -of GTiff -b 1 NETCDF:"$CHL_file":"chlo$MM" tmp.tiff
gdalwarp -overwrite -srcnodata 99999 -s_srs "$crs_latlon" -te ${minlon} ${minlat} ${maxlon} ${maxlat} tmp.tiff tmp1.tiff
gdalwarp -overwrite -s_srs "$crs_latlon" -t_srs "$crsOut" tmp1.tiff tmp2.tiff
gdalwarp -overwrite -srcnodata 99999 -t_srs "$crsOut" -tr $dx $dy -te $xmin $ymin $xmax $ymax -r bilinear -ot Float32 -of Gtiff tmp2.tiff tmp3.tiff
gdal_translate -of NETCDF tmp3.tiff CHLO.nc
rm tmp*

echo "PREP."
echo "   Clip land shapefile"
ogr2ogr -skipfailures -f GPKG land1.gpkg $LAND -clipsrc ${minlon} ${minlat} ${maxlon} ${maxlat}
echo "   Reproject cropped shapefile"
ogr2ogr -skipfailures -f GPKG land2.gpkg land1.gpkg -t_srs "$crsOut" -clipsrc ${xmin} ${ymin} ${xmax} ${ymax}

echo "Create OPEN." #(pixels at open ocean)
gdal_rasterize -burn 0 -a_nodata -999 -init 1 -ot Int32 -tr $dx $dy -te $xmin $ymin $xmax $ymax -of NetCDF land2.gpkg OPEN.nc #

echo "Create SURF." #(fraction of surface at less than 50m of the coastline)
ogr2ogr -f GPKG coast.gpkg land2.gpkg  -nln coastline -dialect sqlite -sql "SELECT ST_Boundary(geom) AS geom,* FROM land_global"
ogr2ogr -f GPKG buffer.gpkg coast.gpkg -nln buffer    -dialect sqlite -sql "SELECT ST_Buffer(geom, 50) AS geom, * FROM coastline"

#base grid. burn buffer with value == 1.0 on "fine" resolution.
gdal_rasterize -a_nodata -999 -init 0.0 -burn 1.0 -tr 20 20 -te $xmin $ymin $xmax $ymax -ot Float32 -of GTiff buffer.gpkg buffer_raster.tif
gdalwarp -overwrite -tr $dx $dy -te $xmin $ymin $xmax $ymax -r average -ot Float32 -of GTiff buffer_raster.tif surf.tif 
gdal_translate -of netCDF surf.tif SURF.nc

echo "cleaning intermediate files.."
rm buffer_raster.tif surf.tif land1.gpkg land2.gpkg coast.gpkg buffer.gpkg

########################################################################################################################################
# Merge netcdfs. and change dim,var and attr names to ioapi/cmaq format.
echo "Change var.  names.."
ncrename -v Band1,OPEN OPEN.nc
ncrename -v Band1,SURF SURF.nc
ncrename -v chlo$MM,CHLO CHLO.nc
ncrename -v dms$MM,DMS DMS.nc
#

if [ -f ./ocean2cmaq.exe ]
then

cat << EOF >> namelist.ocean
&control
  start_date="$start_date",
  griddesc_file="$GRIDDESC_file",
  gridname="$GRIDNAME",
  open_file="OPEN.nc",
  surf_file="SURF.nc",
  chlo_file="CHLO.nc",
  dms_file="DMS.nc"
/
EOF
	./ocean2cmaq.exe < namelist.ocean
else
	echo "ERROR: ocean2cmaq.exe not found."
	echo "       In order to create it do:"
	echo "        > cd src                "
	echo "        > make                  "
fi;
