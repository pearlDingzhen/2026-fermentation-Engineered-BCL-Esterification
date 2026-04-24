#!/bin/bash

for i in di_yc_ywc_ph7
do
  pdb4amber -i ${i}.pdb -o ${i}_amber.pdb -y --dry

  cat << EOF >> ${i}.tleap
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff2
loadoff YC.lib
loadamberparams YC.frcmod
loadAmberPrep SEO.prepin
loadAmberParams frcmod2.SEO
loadAmberParams frcmod1.SEO
${i} = loadpdb ${i}_amber.pdb
saveAmberParm ${i} ${i}.nosol.top ${i}.nosol.inpcrd
solvateBox ${i} TIP3PBOX 8
addIonsRand ${i} Na+ 32
addIonsRand ${i} Cl- 0
charge ${i}
saveAmberParm ${i} ${i}.top ${i}.inpcrd
savepdb ${i} ${i}.check.pdb
quit
EOF

  # Execute tleap and subsequent commands with the updated variable names
  tleap -f ${i}.tleap
  parmed ${i}.top << EOFPAR
    HMassRepartition
    parmout ${i}.hmr.top
    go
EOFPAR

  python a_mod_top_RSFF2C_strict.py ${i}
done