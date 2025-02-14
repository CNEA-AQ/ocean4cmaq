#!/bin/bash
#
# Create OPEN and SURF files for CMAQ from 
# a land & coastline file
#
#---------------------------------------
#INPUTs:
start_date="2024-06-11"
LAND=original/land_10m_global.gpkg
DMS_file=original/DMS_monthly_global_2011.nc
CHL_file=original/CHLO_monthly_global_2019.nc
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
   read xc yc ellipsoidh <<<$( gdaltransform -s_srs "epsg:4326" -t_srs "${crsOut}" <<< $( echo "${clon} ${clat}" ) )
   xmin=$(bc -l <<<"$xorig")
   xmax=$(bc -l <<<"($xorig+ $dx*$nx)")
   ymin=$(bc -l <<<"$yorig")
   ymax=$(bc -l <<<"($yorig+ $dy*$ny)")

   read lon[1] lat[1] ell <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${crsOut}" <<< $( echo "${xmin} ${ymin}" ) )
   read lon[2] lat[2] ell <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${crsOut}" <<< $( echo "${xmax} ${ymax}" ) )
   read lon[3] lat[3] ell <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${crsOut}" <<< $( echo "${xmin} ${ymax}" ) )
   read lon[4] lat[4] ell <<<$( gdaltransform -t_srs "epsg:4326" -s_srs "${crsOut}" <<< $( echo "${xmax} ${ymin}" ) )
   
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
#clip
ogr2ogr -skipfailures -f GPKG land1.gpkg $LAND -clipsrc ${minlon} ${minlat} ${maxlon} ${maxlat}
#reproject
ogr2ogr -skipfailures -f GPKG land2.gpkg land1.gpkg -t_srs "$crsOut" -clipsrc ${xmin} ${ymin} ${xmax} ${ymax}

echo "Create OPEN." #(pixels at open ocean)
gdal_rasterize -burn 0 -a_nodata -999 -init 1 -ot Int32 -tr $dx $dy -te $xmin $ymin $xmax $ymax -of NetCDF land2.gpkg OPEN.nc #

echo "Create SURF." #(fraction of surface at less than 50m of the coastline)
ogr2ogr -f GPKG coast.gpkg land2.gpkg  -nln coastline -dialect sqlite -sql "SELECT ST_Boundary(geom) AS geom,* FROM land_global"
ogr2ogr -f GPKG buffer.gpkg coast.gpkg -nln buffer    -dialect sqlite -sql "SELECT ST_Buffer(geom, 50) AS geom, * FROM coastline"

#base grid. burn buffer with value == 1.0 on "fine" resolution.
gdal_rasterize -a_nodata -999 -init 0.0 -burn 1.0 -tr 25 25 -te $xmin $ymin $xmax $ymax -ot Float32 -of GTiff buffer.gpkg buffer_raster.tif
gdalwarp -overwrite -tr $dx $dy -te $xmin $ymin $xmax $ymax -r average -ot Float32 -of GTiff buffer_raster.tif surf.tif 
gdal_translate -of netCDF surf.tif SURF.nc

echo "cleaning intermediate files.."
rm buffer_raster.tif surf.tif land1.gpkg land2.gpkg coast.gpkg buffer.gpkg

# Merge netcdfs. and change dim,var and attr names to ioapi/cmaq format.
echo "Change var.  names.."
ncrename -v Band1,OPEN OPEN.nc
ncrename -v Band1,SURF SURF.nc
ncrename -v chlo$MM,CHLO CHLO.nc
ncrename -v dms$MM,DMS DMS.nc
echo "Merge [ OPEN.nc, SURF.nc, CHLO.nc and DMS.nc] files."
ncks -A SURF.nc -o OPEN.nc
ncks -A DMS.nc -o OPEN.nc
ncks -A CHLO.nc -o OPEN.nc
echo "Change dim.  names.."
ncrename -d x,NCOLS OPEN.nc
ncrename -d y,NROWS OPEN.nc
echo "Change attr. names.."
ncatted -O -a IOAPI_VERSION,global,a,c,"ioapi-3.2: \$Id: init3" OPEN.nc
ncatted -O -a EXEC_ID,global,a,c,"???????????????? " OPEN.nc
ncatted -O -a FTYPE,global,a,l,"1" OPEN.nc
ncatted -O -a SDATE,global,a,l,"0" OPEN.nc
ncatted -O -a STIME,global,a,l,"0" OPEN.nc
ncatted -O -a WDATE,global,a,l,$YYYY$MM$DD OPEN.nc
ncatted -O -a WTIME,global,a,l,0 OPEN.nc
ncatted -O -a CDATE,global,a,l,$YYYY$MM$DD OPEN.nc
ncatted -O -a CTIME,global,a,l,0 OPEN.nc
ncatted -O -a TSTEP,global,a,l,10000 OPEN.nc
ncatted -O -a NTHIK,global,a,l,1 OPEN.nc
ncatted -O -a NCOLS,global,a,l,$nx OPEN.nc
ncatted -O -a NROWS,global,a,l,$ny OPEN.nc
ncatted -O -a NLAYS,global,a,l,1 OPEN.nc
ncatted -O -a NVARS,global,a,l,2 OPEN.nc
ncatted -O -a P_ALP,global,a,f,$alp OPEN.nc
ncatted -O -a P_BET,global,a,f,$bet OPEN.nc
ncatted -O -a P_GAM,global,a,f,$gam OPEN.nc
ncatted -O -a XCENT,global,a,f,$xc OPEN.nc
ncatted -O -a YCENT,global,a,f,$yc OPEN.nc
ncatted -O -a XORIG,global,a,f,$xmin OPEN.nc
ncatted -O -a YORIG,global,a,f,$ymin OPEN.nc
ncatted -O -a XCELL,global,a,f,$dx OPEN.nc
ncatted -O -a YCELL,global,a,f,$dy OPEN.nc
ncatted -O -a VGTYP,global,a,s,-9999 OPEN.nc
ncatted -O -a VGTOP,global,a,f,5000 OPEN.nc
ncatted -O -a VGLVLS,global,a,f,1 OPEN.nc
ncatted -O -a GDNAM,global,a,c,"$GRIDNAME" OPEN.nc
ncatted -O -a UPNAM,global,a,c,"attribute_value" OPEN.nc
ncatted -O -a VAR-LIST,global,a,c,"OPEN SURF" OPEN.nc
ncatted -O -a FILEDESC,global,a,c,"Ocean file for sea spray emissions." OPEN.nc
ncatted -O -a HISTORY,global,a,c,"" OPEN.nc
ncatted -O -a history,global,a,c,"" OPEN.nc


echo "Cleaning directory.."
mv OPEN.nc ocean_${MM}_${GRIDNAME}.nc 
rm SURF.nc CHLO.nc DMS.nc

echo "end."
