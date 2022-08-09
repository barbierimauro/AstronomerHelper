#!/bin/bash
# 2011, mauro barbieri, maurobarbieri.science@gmail.com
#
# modify the following parameters according to your
# installation of TEMPO2 and your observatory

T2BIN=/astrosoft/tempo2/bin/tempo2 # path to tempo2 executable
CODE_SITE=ae                       # observatory code site
EPHEM=DE405                        # solar system ephemeris to use
dt=0.01                            # time step to use, in days
JD_START=50000.321421412           # starting MJD
JD_END=50005.21421                 # ending MJD


# DO NOT MODIFY THE FOLLOWING PART OF THE SCRIPT !!!!


TOA=toa.$CODE_SITE.tim
TOA_SMALL=toa_small.$CODE_SITE.tim

CLOCK_CORR_GPS=$CODE_SITE"2gps.clk gps2utc.clk"

IN_T1=tempo1.$EPHEM.tdb.$CODE_SITE.par
OUT_T1=check.ephem.$EPHEM.tempo1.$CODE_SITE.tdb

IN_T2=tempo2.$EPHEM.tdb.$CODE_SITE.par
OUT_T2=check.ephem.$EPHEM.tempo2.$CODE_SITE.tdb

IN_T3=tempo2.$EPHEM.tcb.$CODE_SITE.par
OUT_T3=check.ephem.$EPHEM.tempo2.$CODE_SITE.tcb

OUT_CFR_TDB=check.ephem.$EPHEM.tempo1vs2.$CODE_SITE.tdb




leng=$(echo $JD_END "-" $JD_START | bc)
nstep=$(echo $leng "/" $dt | bc)
echo "creating TIM files with " $nstep " points"


echo "FORMAT 1" > $TOA
i=0
j=0
JD=0
while [ $i -lt $nstep ]
do
   let i=i+1
   let j=i-1
   time_aux1=$(echo $dt "*" $j | bc )
   JD=$(echo $time_aux1 "+" $JD_START | bc )
   echo  "null_string 1440.0 " $JD " 0.1 "$CODE_SITE >> $TOA
done


head $TOA > $TOA_SMALL

echo "making PAR files"


cat << EOD1> $IN_T1
PSRJ           J05314+2200
RAJ             05 34 31.95
DECJ           +22 00 52.1
PEPOCH         51544.0
POSEPOCH       51544.0
DM             0
DMEPOCH        51544.0
TEMPO1
EPHEM          $EPHEM
CLK_CORR_CHAIN $CLOCK_CORR_GPS utc2tai.clk tai2tt_tai.clk
CORRECT_TROPOSPHERE N
EOD1

cat <<EOD2> $IN_T2
PSRJ           J05314+2200
RAJ             05 34 31.95
DECJ           +22 00 52.1
PEPOCH         51544.0
POSEPOCH       51544.0
DM             0
DMEPOCH        51544.0
UNITS          TDB
EPHEM          $EPHEM
CLK_CORR_CHAIN $CLOCK_CORR_GPS utc2tai.clk tai2tt_tai.clk
CORRECT_TROPOSPHERE N
EOD2

cat <<EOD3> $IN_T3
PSRJ           J05314+2200
RAJ             05 34 31.95
DECJ           +22 00 52.1
PEPOCH         51544.0
POSEPOCH       51544.0
DM             0
DMEPOCH        51544.0
UNITS          TCB
EPHEM          $EPHEM
CLK_CORR_CHAIN $CLOCK_CORR_GPS utc2tai.clk tai2tt_tai.clk
CORRECT_TROPOSPHERE N
EOD3

echo "barycentering in TEMPO1(TDB)"
$T2BIN -output general2 -nobs 410000 -npsr 1 -f $IN_T1 $TOA -s "{sat} {bat} {roemer} {shapiro} {tt} {tt2tb} {clock} {earth_ssb1} {earth_ssb2} {earth_ssb3} {sun_earth1} {sun_earth2} {sun_earth3}\n"> $OUT_T1

echo "barycentering in TEMPO2(TDB)"
$T2BIN -output general2 -nobs 410000 -npsr 1 -f $IN_T2 $TOA -s "{sat} {bat} {roemer} {shapiro} {tt} {tt2tb} {clock} {earth_ssb1} {earth_ssb2} {earth_ssb3} {sun_earth1} {sun_earth2} {sun_earth3}\n"> $OUT_T2

echo "barycentering in TEMPO2(TCB)"
$T2BIN -output general2 -nobs 410000 -npsr 1 -f $IN_T3 $TOA -s "{sat} {bat} {roemer} {shapiro} {tt} {tt2tb} {clock} {earth_ssb1} {earth_ssb2} {earth_ssb3} {sun_earth1} {sun_earth2} {sun_earth3}\n"> $OUT_T3

echo "#sat1 bat1 roemer1 shapiro1 tt1 tt2tb1 clock1 earth_ssb11 earth_ssb21 earth_ssb31 sun_earth11 sun_earth21 sun_earth31 sat2 bat2 roemer2 shapiro2 tt2 tt2tb2 clock2 earth_ssb12 earth_ssb22 earth_ssb32 sun_earth12 sun_earth22 sun_earth32" > $OUT_CFR_TDB
grep -v 'i' $OUT_T1 | grep -v 'ret' | grep . > tmp_file  && mv tmp_file $OUT_T1
grep -v 'i' $OUT_T2 | grep -v 'ret' | grep . > tmp_file  && mv tmp_file $OUT_T2
paste $OUT_T1 $OUT_T2 | awk '{printf "%40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f %40.24f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26}'>> $OUT_CFR_TDB



# 1 14 sat
# 2 15 bat
# 3 16 roemer
# 4 17 shapiro
# 5 18 tt
# 6 19 tt2tb
# 7 20 clock
# 8 21 earth_ssb1
# 9 22 earth_ssb2
#10 23 earth_ssb3
#11 24 sun_earth1
#12 25 sun_earth2
#13 26 sun_earth3


GCOM=cfr.T1vsT2.tdb.$EPHEM.$CODE_SITE
GNUPLOT_SCRIPT=$GCOM.plt


echo "making graphics PS and PNG"
echo "set terminal png notransparent size 1000,750 enhanced" > $GNUPLOT_SCRIPT
echo "unset k" >> $GNUPLOT_SCRIPT
echo "set xlabel 'MJD'">> $GNUPLOT_SCRIPT
echo "set ylabel 's'">> $GNUPLOT_SCRIPT
echo "set xrange ["$time_start":"$time_end"]" >> $GNUPLOT_SCRIPT
# roemer
echo "set out '"$GCOM".roemer.png'" >> $GNUPLOT_SCRIPT
echo "plot '"$OUT_CFR_TDB"' u 1:(\$3-\$16) w p pt 7">>  $GNUPLOT_SCRIPT
# shapiro
echo "set out '"$GCOM".shapiro.png'" >> $GNUPLOT_SCRIPT
echo "plot '"$OUT_CFR_TDB"' u 1:(\$4-\$17) w p pt 7" >>  $GNUPLOT_SCRIPT
echo "set terminal postscript enhanced color" >> $GNUPLOT_SCRIPT
echo "unset k" >> $GNUPLOT_SCRIPT
echo "set xlabel 'MJD'">> $GNUPLOT_SCRIPT
echo "set ylabel 's'">> $GNUPLOT_SCRIPT
echo "set xrange ["$time_start":"$time_end"]" >> $GNUPLOT_SCRIPT
# roemer
echo "set out '"$GCOM".roemer.ps'" >> $GNUPLOT_SCRIPT
echo "plot '"$OUT_CFR_TDB"' u 1:(\$3-\$16) w p pt 7">>  $GNUPLOT_SCRIPT
# shapiro
echo "set out '"$GCOM".shapiro.ps'" >> $GNUPLOT_SCRIPT
echo "plot '"$OUT_CFR_TDB"' u 1:(\$4-\$17) w p pt 7" >>  $GNUPLOT_SCRIPT

gnuplot $GNUPLOT_SCRIPT
display $GCOM.roemer.png &
display $GCOM.shapiro.png &
