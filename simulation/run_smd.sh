#!/bin/bash
# Steered MD of MAP_3 in water at 300 K and 10 K (GROMACS).
# Force field AMBER99SB-ILDN, TIP3P water, 0.15 M NaCl, pulling along x at 0.01 nm/ps.
set -e
GMX=${GMX:-gmx}

# --- 1. Topology and simulation box (elongated along x, the pulling direction) ---
$GMX pdb2gmx -f map_3.pdb -i map_3_posre.itp -p map_3.top -o map_3.gro \
             -ff amber99sb-ildn -ignh -heavyh -water tip3p
$GMX editconf -f map_3.gro -o newbox.gro -center 3.280 2.181 2.4775 -box 18 6.560 4.362

# --- 2. Solvation and ions ---
$GMX solvate -cp newbox.gro -o solvated.gro -p map_3.top
touch dummy.mdp
$GMX grompp -f dummy.mdp -c solvated.gro -p map_3.top -o ions.tpr
echo "SOL" | $GMX genion -s ions.tpr -p map_3.top -o box_ions.gro -conc 0.15 -neutral

# --- 3. Pull groups: residue 1 is pulled, residue 18 is the reference (restrained) ---
echo q | $GMX make_ndx -f box_ions.gro -o index.ndx
N=$(grep -c '^\[' index.ndx)           # number of default groups -> the new ones get indices N and N+1
echo -e "r 1 \n name $N pull \n r 18 \n name $((N+1)) ref \n q \n" | $GMX make_ndx -f box_ions.gro -n index.ndx -o index.ndx
echo "ref" | $GMX genrestr -f box_ions.gro -n index.ndx -o posre_ref.itp
POSRES_LINES='#ifdef POSRES_ref\n#include "posre_ref.itp"\n#endif\n'
# insert the restraint right after the protein position restraints (inside the protein moleculetype)
sed -i "/^#include \"map_3_posre.itp\"/,/^#endif/{/^#endif/a\\
\\
#ifdef POSRES_ref\\
#include \"posre_ref.itp\"\\
#endif
}" map_3.top

# --- 4. Energy minimization ---
$GMX grompp -f mdp/minimization.mdp -p map_3.top -c box_ions.gro -o em.tpr -maxwarn 1
$GMX mdrun -s em.tpr -deffnm em -v -ntmpi 1

# --- 5. Equilibration and pulling, at each temperature ---
for T in 300K 10K; do
    mkdir -p $T
    $GMX grompp -f mdp/equilibration_$T.mdp -c em.gro -p map_3.top -o $T/eq.tpr -maxwarn 1
    $GMX mdrun -deffnm $T/eq -s $T/eq.tpr -v
    echo -e "Temperature \n 0" | $GMX energy -f $T/eq.edr -o $T/temperature.xvg

    $GMX grompp -f mdp/pulling_$T.mdp -c $T/eq.gro -r $T/eq.gro -t $T/eq.cpt \
                -p map_3.top -n index.ndx -o $T/pull.tpr -maxwarn 2
    $GMX mdrun -deffnm $T/pull -s $T/pull.tpr -v
done
# Outputs used in the analysis: $T/pull_pullf.xvg (force) and $T/pull_pullx.xvg (COM distance)
