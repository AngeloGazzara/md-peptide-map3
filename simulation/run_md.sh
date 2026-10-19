#!/bin/bash
# Equilibrium MD of MAP_3 in water, 2 ns at 300 K (GROMACS), followed by the trajectory analysis.
# Force field AMBER99SB-ILDN, TIP3P water, 0.15 M NaCl. The .xvg outputs are the ones in ../data/md.
set -e
GMX=${GMX:-gmx}

# --- 1. Topology, dodecahedral box, water and ions ---
$GMX pdb2gmx -f map_3.pdb -i map3_posre.itp -p map3.top -o map3.gro \
            -ff amber99sb-ildn -ignh -heavyh -water tip3p
$GMX editconf -f map3.gro -o dodecahedron.gro -c -d 0.8 -bt dodecahedron   # fewer waters than a cubic box
$GMX solvate -cp dodecahedron.gro -o dodecahedron_solvated.gro -p map3.top
touch dummy.mdp
$GMX grompp -f dummy.mdp -c dodecahedron_solvated.gro -p map3.top -o ions.tpr   # dummy.mdp: empty file
echo "SOL" | $GMX genion -s ions.tpr -p map3.top -o box_ions.gro -conc 0.15 -neutral

# --- 2. Energy minimization ---
$GMX grompp -f mdp/minimization.mdp -p map3.top -c box_ions.gro -o em.tpr -maxwarn 1
$GMX mdrun -s em.tpr -deffnm em -v
echo -e "Potential \n 0" | $GMX energy -f em.edr -o potential.xvg

# --- 3. Production MD: 2 ns at 300 K (v-rescale thermostat) ---
$GMX grompp -f mdp/md_300K.mdp -c em.gro -p map3.top -o prod.tpr -maxwarn 1
$GMX mdrun -deffnm prod -s prod.tpr -v

# --- 4. Analysis (outputs in ../data/md) ---
echo -e "Temperature \n Pressure \n 0" | $GMX energy -f prod.edr -o Temp_Press.xvg
echo -e "0 \n" | $GMX trjconv -s prod.tpr -f prod.trr -pbc mol -o noPBC.xtc -ur compact
echo -e "Backbone\nBackbone" | $GMX rms -s map_3.pdb -f noPBC.xtc -o rmsd.xvg -tu ns
echo -e "System" | $GMX trjconv -s prod.tpr -f noPBC.xtc -dump 1.5 -o Struct15.pdb -tu ns
echo -e "Backbone\nBackbone" | $GMX rms -s Struct15.pdb -f noPBC.xtc -o rmsd15.xvg -tu ns
echo "C-alpha" | $GMX rmsf -f noPBC.xtc -s prod.tpr -o rmsf.xvg -res
$GMX rama -f noPBC.xtc -s prod.tpr -o rama.xvg
echo -e "r 1\nname 10 E1\nr 18\nname 11 E2\nq" | $GMX make_ndx -f map_3.pdb -o index_NEW.ndx
$GMX distance -f noPBC.xtc -s prod.tpr -select 'com of group "E1" plus com of group "E2"' \
             -n index_NEW.ndx -oav distance.xvg
echo -e "Protein \n Water \n " | $GMX hbond -f noPBC.xtc -s prod.tpr -num hb_water_protein.xvg
echo -e "Protein \n Protein \n " | $GMX hbond -f noPBC.xtc -s prod.tpr -num hb_protein_protein.xvg
