awk '{print $1*1e10}' *??_det:1_throughput_atmosphere.ascii > wl.dat
awk '{print $2}' *??_det:1_throughput_detector.ascii > det.dat
awk '{print $2}' *??_det:1_throughput_fibinj.ascii > fib.dat
awk '{print $2}' *??_det:1_throughput_instrument.ascii > ins.dat
awk '{print $2}' *??_det:1_throughput_telescope.ascii > tel.dat
awk '{print $2}' *??_det:1_throughput_blaze.ascii > bla.dat
awk '{print $2}' *??_det:1_dispersion_dispersion.ascii > dis.dat
awk '{print $2}' *??_det:1_throughput_atmosphere.ascii > atm.dat
awk '{print $2}' *??_det:1_sed_sky.ascii > sky.dat

echo 'wl,trasm,det,fib,ins,tel,bla,dis,atm,sky'> harps.csv
paste wl.dat det.dat fib.dat ins.dat tel.dat bla.dat dis.dat atm.dat sky.dat | sort -nk1 > harps.dat
awk '{printf"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",$1,$2*$3*$4*$5,$2,$3,$4,$5,$6,$7,$8,$9}' harps.dat >> harps.csv

wc -l harps.csv harps.dat
