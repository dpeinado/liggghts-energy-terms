#!/bin/bash
Vangle=0
Vmod=10
#model="gran/hertz/incremental/energy  1 0 "
model="gran/hooke/history/energy  1 "
#model="gran/hertz/history/energy  1 "
#modelName="SubStep_hertz_incremental_energy"
modelName="SubStep_hooke_history_energy"
#modelName="SubStep_hertz_integral_energy"
for option in 0 # 0 2 4
do
cofI=1
cofD=0
enI=1
enD=0
poI=0
poD=3
rootName=${modelName}-${option}_En${enI}_${enD}_COF${cofI}_${cofD}_PO${poI}_${poD}
echo pair_style ${model} ${option}
for a in 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85
do
	Vangle=0
	Vmod=10
	Vx="`calc "${Vmod}*cos(3.141592/180.0*(${a}+180))"`"
	Vy="`calc "${Vmod}*sin(3.141592/180.0*(${a}+180))"`"
	./lmp_serial_debug << EOF
	echo both
	dimension 3
	atom_style granular
	boundary f f f
	newton off
	region dominio block 0 100 0 100 0 100 units box
	create_box 1 dominio
	create_atoms 1 single 50 50.0 0.5001 units box
	set atom 1 vz ${Vx}
	set atom 1 vy ${Vy}
	set atom 1 diameter 1
	set atom 1 density  1.90985485103132161955691
	set atom 1 omegay   0
	neighbor		1 bin
	neigh_modify	delay 0
	fix 			m1 all property/global youngsModulus peratomtype 763.94194041252864782277e7
	fix 			m2 all property/global poissonsRatio peratomtype ${poI}.${poD}
	fix 			m3 all property/global coefficientRestitution peratomtypepair 1 ${enI}.${enD}
	fix 			m4 all property/global coefficientFriction peratomtypepair 1 ${cofI}.${cofD}
	fix 		 	m5 all property/global characteristicVelocity scalar 10
	pair_style		${model} ${option}
	pair_coeff		* *
	fix			zwall all wall/${model} 0 zplane 0 NULL 1
	communicate		single vel yes
	fix    			1 all nve/sphere
	timestep		0.0000001
	run 0
	compute			rot_e all erotate/sphere
	compute		    	epotN all reduce sum f_CPEn
	compute			edisN all reduce sum f_CDEn
	compute			epotT all reduce sum f_CPEt
	compute			edisTV all reduce sum f_CDEVt
	compute			edisTF all reduce sum f_CDEFt
	compute			workT all reduce sum f_CTFW
	compute			edisH all reduce sum f_DEH
	variable		eCon equal "ke + c_rot_e + c_epotN"
	variable		eTot equal "v_eCon + c_edisN + c_workT + c_edisTF+c_edisTV+c_edisH"
	variable		eKin equal "ke"
	thermo_style		custom step atoms ke c_rot_e c_epotN c_epotT c_edisN c_edisTV c_edisTF c_workT v_eCon v_eTot c_edisH
	thermo			1000
	thermo_modify	lost ignore norm no
	compute_modify	thermo_temp dynamic yes

	dump			mydmp all custom 10 files/dump-${rootName}.${a} id type mass x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius f_CPEn f_CDEn f_CPEt f_CDEVt f_CDEFt f_CTFW f_DEH
	run			    5500
EOF
python ~/liggghts-energy-terms/pyPost/pyGRAPH.py files/dump-${rootName}.${a} files/${rootName}_${a}_
sh plot_collision_energy_wall.sh "${rootName} ANG = ${a}" ${rootName}_${a}_
done
python ~/liggghts-energy-terms/pyPost/pyCW.py "files/dump-${rootName}.*" files/plot_${rootName}.dat 0 ${cofI}"."${cofD}
done
echo
echo
