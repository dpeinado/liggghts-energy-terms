#!/bin/bash
#model="gran/hertz/incremental/energy  1 0 "
model="gran/hooke/history/energy  1 "
#model="gran/hertz/history/energy  1 "
#modelName="SubStep_hertz_incremental_energy"
modelName="SubStep_hooke_history_energy"
#modelName="SubStep_hertz_integral_energy"
#for option in 0 2 4
#do
option=0
cofI=0
cofD=5
enI=1
enD=0
poI=0
poD=3
Vmod=1
RelM=1000
rootName=${modelName}-${option}_En${enI}_${enD}_COF${cofI}_${cofD}_PO${poI}_${poD}
echo pair_style ${model} ${option}
#for a in 0 #1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
for a in 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85
do
#	angulo="`calc "atan(${a}/20.0)"`"
	anguloB="`calc "${a}*atan(1.0)/45.0"`"
	angulo="`echo ${anguloB} | sed 's/^~//g'`"
	echo ${angulo}, "angulo"
	Vx1B="`calc "${Vmod}*cos(${angulo})"`"
	echo ${Vx1B}, "Vx1B"
	Vy1B="`calc "${Vmod}*sin(${angulo})"`"
	echo ${Vy1B}, "Vy1B"
	Vx1="`echo ${Vx1B} | sed 's/^~//g'`"
	Vy1="`echo ${Vy1B} | sed 's/^~//g'`"
	echo ${Vx1}, "Vx1"
	echo ${Vy1}, "Vy1"
	Vx2B="`calc "(-${Vx1}*${RelM})"`"
	echo ${Vx2B}, "Vx2B"
	Vy2B="`calc "(-${Vy1}*${RelM})"`"
	Vx2="`echo ${Vx2B} | sed 's/^~//g'`"
	Vy2="`echo ${Vy2B} | sed 's/^~//g'`"
	echo ${angulo}, ${Vx1},${Vy1}, ${Vx2}, ${Vy2}
#	read -n 1 -s "Press any key ..."
	./lmp_serial_debug << EOF
	echo both
	dimension 3
	atom_style granular
	boundary f f f
	newton off
	region dominio block 0 100 0 100 0 100 units box
	create_box 2 dominio
	create_atoms 1 single 49.0 50.0 50 units box
	create_atoms 2 single 54.5 50.0 50 units box
	set atom 1 vx ${Vx1}
	set atom 1 vy ${Vy1}
	set atom 1 diameter 10
	set atom 1 density  1.90985485103132161955691
	set atom 1 omegay   0
	set	atom 2 vx  ${Vx2}
	set	atom 2 vy ${Vy2}
	set atom 2 diameter 1
	set atom 2 density  1.90985485103132161955691
	neighbor		1 bin
	neigh_modify	delay 0
	fix 			m1 all property/global youngsModulus peratomtype  763.94194041252864782277e7  763.94194041252864782277e7
	fix 			m2 all property/global poissonsRatio peratomtype ${poI}.${poD} ${poI}.${poD}
	fix 			m3 all property/global coefficientRestitution peratomtypepair 2 ${enI}.${enD} ${enI}.${enD} ${enI}.${enD} ${enI}.${enD}
	fix 			m4 all property/global coefficientFriction peratomtypepair 2 ${cofI}.${cofD} ${cofI}.${cofD} ${cofI}.${cofD} ${cofI}.${cofD}
	fix 		 	m5 all property/global characteristicVelocity scalar 1000
	pair_style		${model} ${option}
	pair_coeff		* *
	communicate		single vel yes
	fix    			1 all nve/sphere
	timestep		0.000000001
	fix 			ts_check all check/timestep/gran 10 0.1 0.1
	run 0
	compute			rot_e all erotate/sphere
	compute		    	epotN all reduce sum f_CPEn
	compute			edisN all reduce sum f_CDEn
	compute			edisTV all reduce sum f_CDEVt
	compute			edisTF all reduce sum f_CDEFt
	compute			workT all reduce sum f_CTFW
	compute			edisH all reduce sum f_DEH
	variable		eCon equal "ke + c_rot_e + c_epotN"
	variable		eTot equal "v_eCon + c_edisN + c_workT + c_edisTF+c_edisTV+c_edisH"
	variable		eKin equal "ke"
	thermo_style		custom step atoms ke c_rot_e c_epotN c_edisN c_edisTV c_edisTF c_workT v_eCon v_eTot c_edisH
	thermo			100000
	thermo_modify		lost ignore norm no
	compute_modify		thermo_temp dynamic yes
	dump			mydmp all custom 1000 files/dump-${rootName}.${a} id type mass x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius f_CPEn f_CDEn f_CDEVt f_CDEFt f_CTFW f_DEH
	run			1000000
EOF
python ~/liggghts-energy-terms/pyPost/pyGRAPH.py files/dump-${rootName}.${a} files/${rootName}_${a}_
sh plot_collision_energy.sh "${rootName} ANG = ${a}" ${rootName}_${a}_
done
python ~/liggghts-energy-terms/pyPost/pyCB.py "files/dump-${rootName}.*" files/plot_${rootName}.dat ${cofI}"."${cofD}
#done
echo
echo

#	xsecondB="`calc "49.0+1.0*cos(${angulo})"`"
#	xsecond="`echo ${xsecondB} | sed 's/^~//g'`"
#	ysecondB="`calc "50.0-1.0*sin(${angulo})"`"
#	ysecond="`echo ${ysecondB} | sed 's/^~//g'`"
