#!/bin/bash

mkdir -p min
gmx grompp -f mdp/min.mdp -c test.gro -p ghost.top -o min/min.tpr
cd min
gmx mdrun -deffnm min

cd ..
mkdir -p npt
gmx grompp -f mdp/npt.mdp -c min/min.gro -p ghost.top -o npt/npt.tpr
cd npt 
gmx mdrun -deffnm npt

cd ..
mkdir -p FE
for i in `seq 0 19`
do 
cp mdp/FE.mdp FE/"$i".mdp
sed -i.bak "s|XXX|$i|g" FE/"$i".mdp
done 
rm -rf FE/*.bak

#cp topol.top ghost.top
# Edit ghost.top to add ghost molecule type and change (NA 1 and CL 1) to (ghost 1) at the end

cd FE
for i in `seq 0 19`
do 
mkdir -p topol"$i"
done 
cp ../npt/npt.gro topol0/start.gro 


for i in `seq 0 19`
do 
	gmx grompp -f "$i".mdp -c topol"$i"/start.gro -p ../ghost.top -o topol"$i"/topol.tpr -maxwarn 1
	cd topol"$i"
	gmx mdrun -deffnm topol -nsteps 100000 -v
	if [ $i -ne 19 ]; then
		cp topol.gro ../topol"$((i+1))"/start.gro
	fi
	cd ..
done

mkdir -p dhdl

for i in `seq 0 19`
do 
	cp topol"$i"/topol.xvg dhdl/topol."$i".dhdl.xvg
done
