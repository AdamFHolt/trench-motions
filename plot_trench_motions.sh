#!/bin/bash
predicted=${1-"rms_samodel.lith10e20.0_asthen10e20.75"}
observed=${2-"tnew.sa.dat"}
matches=${3-""}
rms_or_signs=${4-""} # 0 = signs, 1 = rms

ofilepng="$predicted".png
rm $ofile $ofilepng 2> /dev/null

datasets=/home/aholt/Dropbox/Izu_3D/data_stuff/datasets
age_grd=$datasets/age/age.3.6.NaN.grd
vt=data/vt/$observed

gmtset FONT_LABEL 13p,1,black FONT_ANNOT_PRIMARY 11.5p,1,black
gmtset MAP_FRAME_PEN thin,black MAP_TICK_PEN_PRIMARY thin,black MAP_TICK_PEN_SECONDARY thin,black
gmtset MAP_TICK_LENGTH_PRIMARY 4p MAP_TICK_LENGTH_SECONDARY 1.5p PS_CHAR_ENCODING Standard+

reg=-R0/360/-70/70
proj=-Jm0.018i

makecpt -T0/200/10 -I -Cgray > age.cpt
makecpt -T-60/60/5 -Cred2green -D > vt.cpt

psbasemap $reg $proj  -Ba60f30/a30f15:."":Wesn  -K -Y+9.6 > tmp.ps
grdview -Qc $age_grd   -Cage.cpt $reg $proj  -K -P -O >> tmp.ps 
pscoast $reg $proj -Di -K  -O -Gdarkgrey -P >> tmp.ps
psxy $datasets/plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt $reg $proj -O  -K   -W0.5,darkmagenta  >> tmp.ps
psscale --FONT_LABEL=10p,Helvetica,black  --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba50:"age  [Ma]": -O -D0.75i/3.4i/2.5c/0.25ch -E -Cage.cpt  -K >> tmp.ps 
psscale --FONT_LABEL=10p,Helvetica,black  --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba25:"v@-T@-  [mm/yr]": -O -D2i/3.4i/2.5c/0.25ch -E -K -Cvt.cpt  >> tmp.ps 
scale=90 # trench motion, in mm/yr
## lon, lat, color, azimuth, length
gawk -v s=$scale  '{print($1,$2,$6,$3,-1.0*($6/s))}' $vt | psxy $reg $proj  -O  -SV0.14i+e+a55 -W1.3p,darkslategray  -Cvt.cpt -K  >> tmp.ps

psbasemap $reg $proj  -Ba60f30/a30f15:."":WeSn -O -K -Y-9.6 >> tmp.ps
grdview -Qc $age_grd   -Cage.cpt $reg $proj  -K -P -O >> tmp.ps
pscoast $reg $proj -Di -K  -O -Gdarkgrey -P >> tmp.ps
psxy $datasets/plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt $reg $proj -O  -K   -W0.5,darkmagenta  >> tmp.ps
gawk -v s=$scale  '{print($2,$1,$3*10.0,$4,-1.0*(($3*10.0)/s))}' "$predicted".txt | psxy $reg $proj  -O  -SV0.14i+e+a55 -W1.3p,darkslategray  -Cvt.cpt  -K >> tmp.ps
if [ $rms_or_signs -eq 1 ];
then
   makecpt -T0/10/1 -Cseis -D -I > rms.cpt
   gawk '{print($2,$1,$3)}' $matches.txt | psxy $reg $proj -Crms.cpt -Sc0.1c -Wthinnest,grey  -O -K  >> tmp.ps
   psscale --FONT_LABEL=10p,Helvetica,black  --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba2:"misfit  [cm/yr]": -O -D0.75i/3.4i/2.5c/0.25ch -E -Crms.cpt  >> tmp.ps
   rm rms.cpt $matches.txt 2> /dev/null
else
   makecpt -T0/1/0.5 -Cpolar -D -I > signs.cpt
   gawk '{print($2,$1,$3)}' $matches.txt | psxy $reg $proj -Csigns.cpt -Sc0.1c -Wthinnest,grey  -O -K  >> tmp.ps
   psscale --FONT_LABEL=10p,Helvetica,black  --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba1:"signs match?": -O -D0.75i/3.4i/2.5c/0.25ch -E -Csigns.cpt  >> tmp.ps
   rm signs.cpt $matches.txt 2> /dev/null
fi



eps2eps tmp.ps tmpn.ps
mv tmpn.ps tmp.ps
convert  -density 600 -flatten -rotate 90 tmp.ps $ofilepng
echo plot output in $ofilepng

rm "$predicted".txt 2> /dev/null
rm tmp.ps 2> /dev/null
rm age.cpt vt.cpt 2> /dev/null
