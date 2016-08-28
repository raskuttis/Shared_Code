#!/bin/bash

res="256"
rtime="71:55:00"
masses=("5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4" "5.0e4")
names=("a" "b" "c" "d" "e" "f" "g" "h" "i" "j" "k" "l" "i")
radii=("15.0" "15.0" "15.0" "15.0" "15.0" "15.0" "15.0" "15.0" "15.0" "15.0" "15.0" "15.0" "15.0")
alphas=("2.0" "2.0" "2.0" "2.0" "2.0" "2.0" "2.0" "2.0" "2.0" "2.0" "2.0" "2.0" "2.0")
betas=("0.05" "0.1" "0.2" "0.5" "1.0" "2.0" "5.0" "10.0" "20.0" "50.0" "100.0" "200.0" "50000.0")
suffix="LD_SDs"
fr="0"
sr="1"
psi="2000.0"

mnorm="5.0e4"
rnorm="15.0"
anorm="2.0"
xlimnorm=$(echo "2*$rnorm" | bc -l)
tffnorm="1.76e1"
thstnorm="1.76e-3"
rhofdtnorm="1.76e-2"

for i in {0..12}
do
	fname="${names[$i]}"
	mass="${masses[$i]}"
	alpha="${alphas[$i]}"
	beta="${betas[$i]}"
	radius="${radii[$i]}"
	xlim=$(echo "2*$radius" | bc -l)
	tffnormns=$(echo "$tffnorm" | sed 's/e/*10^/g;s/ /*/' | bc)
	massns=$(echo "$mass" | sed 's/e/*10^/g;s/ /*/' | bc)
	mnormns=$(echo "$mnorm" | sed 's/e/*10^/g;s/ /*/' | bc)
	tffns=$(echo "$tffnormns/sqrt(($massns/$mnormns)/($radius/$rnorm)^3)" | bc -l)
	rhofdt=$(echo "$tffns/1000.0" | bc -l)
	thst=$(echo "$tffns/10000.0" | bc -l)

	echo $i $mass $radius

	newf="sbrun.radpg.run$fname"
	cp sbrun.radpg.based $newf

	sed -i "s/fileextout=UV_M5.0e4_R15.0_N128/fileextout=UV_M5.0e4_R15.0_N$res/g" $newf
	sed -i "s/fileextout=UV_M5.0e4_R15.0/fileextout=UV_M5.0e4_R$radius/g" $newf
	sed -i "s/fileextout=UV_M$mnorm/fileextout=UV_M$mass/g" $newf
	
	sed -i "s/M_GMC=$mnorm/M_GMC=$mass/g" $newf	
	sed -i "s/rcloud=$rnorm/rcloud=$radius/g" $newf	
	sed -i "s/alpha_vir=$anorm/alpha_vir=$alpha/g" $newf	
	
	sed -i "s/fluxrad_out=0/fluxrad_out=$fr/g" $newf	
	sed -i "s/surfd_out=1/surfd_out=$sr/g" $newf	
	sed -i "s/Psi=0.0/Psi=$psi/g" $newf	
	sed -i "s/beta=0.05/beta=$beta/g" $newf	
	sed -i "s/Tf4_B0.05_SDs/Tf4_B0.05_$suffix/g" $newf	
	sed -i "s/Tf4_B0.05/Tf4_B$beta/g" $newf	
	
	sed -i "s/Nx1=128/Nx1=$res/g" $newf	
	sed -i "s/Nx2=128/Nx2=$res/g" $newf	
	sed -i "s/Nx3=128/Nx3=$res/g" $newf	
	sed -i "s/07:55:00/$rtime/g" $newf	

	sed -i "s/tlim=$tffnorm/tlim=$tffns/g" $newf	
	sed -i "s/dt=$thstnorm/dt=$thst/g" $newf	
	sed -i "s/rho_fdt=$rhofdtnorm/rho_fdt=$rhofdt/g" $newf	
	
	sed -i "s/x1min=-$xlimnorm/x1min=-$xlim/g" $newf	
	sed -i "s/x1max=$xlimnorm/x1max=$xlim/g" $newf	
	sed -i "s/x2min=-$xlimnorm/x2min=-$xlim/g" $newf	
	sed -i "s/x2max=$xlimnorm/x2max=$xlim/g" $newf	
	sed -i "s/x3min=-$xlimnorm/x3min=-$xlim/g" $newf	
	sed -i "s/x3max=$xlimnorm/x3max=$xlim/g" $newf	

	sbatch $newf
done
