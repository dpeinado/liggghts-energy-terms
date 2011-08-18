#!/bin/bash
Vangle=0
Vmod=10
#model="gran/hertz/incremental/energy  1 0 "
model="gran/hooke/history/energy  1 "
#modelName="hertz_incremental_energy"
modelName="NonSubStep2_hooke_history_energy"
#for option in 0 1 2
for option in 0 # 0 2 4
do
#option=1
cofI=1
cofD=0
enI=0
enD=7
poI=0
poD=3
rootName=${modelName}-${option}_En${enI}_${enD}_COF${cofI}_${cofD}_PO${poI}_${poD}
echo pair_style ${model} ${option}
for a in 0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85
do
	xsecond="`calc "49.5+cos(${a}*3.141592/180.)"`"
	ysecond="`calc "50.0-sin(${a}*3.141592/180.)"`"
	Vangle=0
	Vmod=10
	Vx="`calc "${Vmod}*cos(3.141592/180.0*${Vangle})"`"
	Vy="`calc "${Vmod}*sin(3.141592/180.0*${Vangle})"`"
	./lmp_non_substep2 << EOF
	echo both
	dimension 3
	atom_style granular
	boundary f f f
	newton off
	region dominio block 0 100 0 100 0 100 units box
	create_box 1 dominio
	create_atoms 1 single 49.5 50.0 50 units box
	create_atoms 1 single ${xsecond} ${ysecond} 50 units box
	set atom 1 vx ${Vx}
	set atom 1 vy ${Vy}
	set atom 1 diameter 1
	set atom 1 density  1.90985485103132161955691
	set atom 1 omegay   0
	set	atom 2 vx 0
	set	atom 2 vy 0
	set atom 2 diameter 1
	set atom 2 density  1.90985485103132161955691
	neighbor		1 bin
	neigh_modify	delay 0
	fix 			m1 all property/global youngsModulus peratomtype 763.94194041252864782277e7
	fix 			m2 all property/global poissonsRatio peratomtype ${poI}.${poD}
	fix 			m3 all property/global coefficientRestitution peratomtypepair 1 ${enI}.${enD}
	fix 			m4 all property/global coefficientFriction peratomtypepair 1 ${cofI}.${cofD}
	fix 		 	m5 all property/global characteristicVelocity scalar 10
	pair_style		${model} ${option}
	pair_coeff		* *
	communicate		single vel yes
	fix    			1 all nve/sphere
	timestep		0.0000001
	run 1
	compute			rot_e all erotate/sphere
	compute		    epotN all reduce sum f_CPEn
	compute			edisN all reduce sum f_CDEn
	compute			epotT all reduce sum f_CPEt
	compute			edisTV all reduce sum f_CDEVt
	compute			edisTF all reduce sum f_CDEFt
	compute			workT all reduce sum f_CTFW
	compute			edisH all reduce sum f_DEH
	variable		eCon equal "ke + c_rot_e + c_epotN"
	variable		eTot equal "v_eCon + c_edisN + c_workT + c_edisTF+c_edisTV+c_edisH"
	variable		eKin equal "ke"
	thermo_style	custom step atoms ke c_rot_e c_epotN c_epotT c_edisN c_edisTV c_edisTF c_workT v_eCon v_eTot c_edisH
	thermo			1000
	thermo_modify	lost ignore norm no
	compute_modify	thermo_temp dynamic yes
	dump			mydmp all custom 10 files/dump-${rootName}.${a} id type mass x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius f_CPEn f_CDEn f_CPEt f_CDEVt f_CDEFt f_CTFW f_DEH
	run			    5500
EOF
python /home/lagar/LIGGGHTS/pyPost/pyGRAPH.py files/dump-${rootName}.${a} files/${rootName}_${a}_
sh plot_collision_energy.sh "${rootName} ANG = ${a}" ${rootName}_${a}_
done
python /home/lagar/LIGGGHTS/pyPost/pyCB.py "files/dump-${rootName}.*" files/plot_${rootName}.dat
done
echo
echo
