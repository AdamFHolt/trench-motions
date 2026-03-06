#!/bin/bash
set -euo pipefail

predicted=${1-"rms_samodel.lith10e20.0_asthen10e20.75"}
observed=${2-"tnew.sa.dat"}
matches=${3-""}
rms_or_signs=${4-""} # 0 = signs, 1 = rms

ofilepng="${predicted}.png"
predicted_txt="${predicted}.txt"
matches_txt="${matches}.txt"
rm -f "$ofilepng" 2>/dev/null || true

gmt_cmd() {
  if command -v "$1" >/dev/null 2>&1; then
    "$@"
  else
    gmt "$@"
  fi
}

gmt_set() {
  if command -v gmtset >/dev/null 2>&1; then
    gmtset "$@"
  else
    gmt set "$@"
  fi
}

if ! command -v gmt >/dev/null 2>&1; then
  echo "Missing GMT executable ('gmt'). Install GMT to generate maps." >&2
  exit 1
fi
if ! command -v gawk >/dev/null 2>&1; then
  echo "Missing gawk executable. Install gawk to generate maps." >&2
  exit 1
fi

datasets="${DATASETS_DIR:-data/gmt_datasets}"
age_grd="$datasets/age/age.3.6.NaN.grd"
pb_file="$datasets/plate_boundaries/bird_PB2002/PB2002_tdiddy.gmt"
vt="data/vt/$observed"

# Support flat dataset layout directly under data/.
if [ ! -f "$age_grd" ] && [ -f "data/age.3.6.NaN.grd" ]; then
  age_grd="data/age.3.6.NaN.grd"
fi
if [ ! -f "$pb_file" ] && [ -f "data/PB2002_tdiddy.gmt" ]; then
  pb_file="data/PB2002_tdiddy.gmt"
fi

if [ ! -f "$vt" ]; then
  echo "Missing observed trench velocity file: $vt" >&2
  exit 1
fi
if [ ! -f "$predicted_txt" ]; then
  echo "Missing predicted file: $predicted_txt" >&2
  exit 1
fi
if [ ! -f "$matches_txt" ]; then
  echo "Missing matches file: $matches_txt" >&2
  exit 1
fi

has_age=0
has_pb=0
if [ -f "$age_grd" ]; then
  has_age=1
else
  echo "warning: age grid not found at $age_grd; plotting without age raster background" >&2
fi
if [ -f "$pb_file" ]; then
  has_pb=1
else
  echo "warning: plate boundary file not found at $pb_file; plotting without PB overlay" >&2
fi

gmt_set FONT_LABEL 13p,1,black FONT_ANNOT_PRIMARY 11.5p,1,black
gmt_set MAP_FRAME_PEN thin,black MAP_TICK_PEN_PRIMARY thin,black MAP_TICK_PEN_SECONDARY thin,black
gmt_set MAP_TICK_LENGTH_PRIMARY 4p MAP_TICK_LENGTH_SECONDARY 1.5p PS_CHAR_ENCODING Standard+

reg=-R0/360/-70/70
proj=-Jm0.018i

if [ "$has_age" -eq 1 ]; then
  gmt_cmd makecpt -T0/200/10 -I -Cgray > age.cpt
fi
gmt_cmd makecpt -T-60/60/5 -Cred2green -D > vt.cpt

gmt_cmd psbasemap $reg $proj -Ba60f30/a30f15:."":Wesn -K -Y+9.6 > tmp.ps
if [ "$has_age" -eq 1 ]; then
  gmt_cmd grdimage "$age_grd" -Cage.cpt $reg $proj -K -P -O >> tmp.ps
fi
gmt_cmd pscoast $reg $proj -Di -K -O -Gdarkgrey -P >> tmp.ps
if [ "$has_pb" -eq 1 ]; then
  gmt_cmd psxy "$pb_file" $reg $proj -O -K -W0.5,darkmagenta >> tmp.ps
fi
if [ "$has_age" -eq 1 ]; then
  gmt_cmd psscale --FONT_LABEL=10p,Helvetica,black --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba50:"age  [Ma]": -O -D0.75i/3.4i/2.5c/0.25ch -E -Cage.cpt -K >> tmp.ps
fi
gmt_cmd psscale --FONT_LABEL=10p,Helvetica,black --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba25:"v@-T@-  [mm/yr]": -O -D2i/3.4i/2.5c/0.25ch -E -K -Cvt.cpt >> tmp.ps

scale=90 # trench motion, in mm/yr
gawk -v s="$scale" '{print($1,$2,$6,$3,-1.0*($6/s))}' "$vt" | gmt_cmd psxy $reg $proj -O -SV0.14i+e+a55 -W1.3p,darkslategray -Cvt.cpt -K >> tmp.ps

gmt_cmd psbasemap $reg $proj -Ba60f30/a30f15:."":WeSn -O -K -Y-9.6 >> tmp.ps
if [ "$has_age" -eq 1 ]; then
  gmt_cmd grdimage "$age_grd" -Cage.cpt $reg $proj -K -P -O >> tmp.ps
fi
gmt_cmd pscoast $reg $proj -Di -K -O -Gdarkgrey -P >> tmp.ps
if [ "$has_pb" -eq 1 ]; then
  gmt_cmd psxy "$pb_file" $reg $proj -O -K -W0.5,darkmagenta >> tmp.ps
fi
gawk -v s="$scale" '{print($2,$1,$3*10.0,$4,-1.0*(($3*10.0)/s))}' "$predicted_txt" | gmt_cmd psxy $reg $proj -O -SV0.14i+e+a55 -W1.3p,darkslategray -Cvt.cpt -K >> tmp.ps

if [ "$rms_or_signs" -eq 1 ]; then
  gmt_cmd makecpt -T0/10/1 -Cseis -D -I > rms.cpt
  gawk '{print($2,$1,$3)}' "$matches_txt" | gmt_cmd psxy $reg $proj -Crms.cpt -Sc0.1c -Wthinnest,grey -O -K >> tmp.ps
  gmt_cmd psscale --FONT_LABEL=10p,Helvetica,black --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba2:"misfit  [cm/yr]": -O -D0.75i/3.4i/2.5c/0.25ch -E -Crms.cpt >> tmp.ps
  rm -f rms.cpt 2>/dev/null || true
else
  gmt_cmd makecpt -T0/1/0.5 -Cpolar -D -I > signs.cpt
  gawk '{print($2,$1,$3)}' "$matches_txt" | gmt_cmd psxy $reg $proj -Csigns.cpt -Sc0.1c -Wthinnest,grey -O -K >> tmp.ps
  gmt_cmd psscale --FONT_LABEL=10p,Helvetica,black --FONT_ANNOT_PRIMARY=8p,Helvetica,black -Ba1:"signs match?": -O -D0.75i/3.4i/2.5c/0.25ch -E -Csigns.cpt >> tmp.ps
  rm -f signs.cpt 2>/dev/null || true
fi

if command -v eps2eps >/dev/null 2>&1 && command -v convert >/dev/null 2>&1; then
  eps2eps tmp.ps tmpn.ps
  mv tmpn.ps tmp.ps
  convert -density 600 -flatten -rotate 90 tmp.ps "$ofilepng"
else
  gmt psconvert tmp.ps -A -Tg -E300 -F"$predicted"
fi

echo "plot output in $ofilepng"

rm -f "$matches_txt" 2>/dev/null || true
rm -f tmp.ps tmpn.ps 2>/dev/null || true
rm -f age.cpt vt.cpt 2>/dev/null || true
