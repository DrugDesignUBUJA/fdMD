#!/bin/csh -v
# ......................
# fdMD_MMGBSA_Send.cmd .
# ......................

# top file without repulsion
set nom_top  = 2ohk_onelig_norep.top

#defining the Protein 
set ini_res = 1 
### ....  Without the ligand.......
set ifi_res = 387 
set num_res_prot = @ifi_res

# Number of Ligands in the simulation
@ num_lig      =  26

### ....  Do PB calculations  ?
set pb = 'NO'

set prefix   = '2ohk'
set lig_parm = 'lig.parm'
set lig_prep = 'lig.prep'

@ startf = 1     
@ endf   = 200   
@ nfreq  = 1

set indi  = 1.0
set exdi = 80.0

set num_proc  = 8
  
set AMBERHOME = '/aplic/amber/amber16_ompi'

@ last_lig  = $num_res_prot + $num_lig
@ count_lig = $num_res_prot + 1

@ pdb_exist = 0

echo $count_lig $last_lig

while ( $count_lig <= $last_lig )
  set name_pdb =  'lig_'$count_lig'_LAST.pdb'
  if ( -e $name_pdb ) then
    set initial_pdb = $name_pdb
    @ pdb_exist = 1
  endif
  @ count_lig++
end
if ( $pdb_exist == 0 ) then
  echo "NO initial PDB = Prot + Lig "
  echo "STOP "
  exit
endif

rm send_mmpbsa_all

set nom_top      = $prefix'_onelig_norep.top'
set nom_rst      = $prefix'_onelig_norep.rst'
set nom_pdb_leap = $prefix'_onelig_norep_leap.pdb'

########
#LEAP
########

mkdir leap
cd    leap

cp ../$lig_parm     .
cp ../$lig_prep     .
cp ../$initial_pdb  .

#############
#file tleap
#############

echo 'source leaprc.protein.ff14SB'                                     >>tleap.script
echo 'source leaprc.water.tip3p'                                        >>tleap.script
echo 'source leaprc.gaff2'                                              >>tleap.script
echo 'set default pbradii mbondi2'                                      >>tleap.script
echo 'loadAmberParams frcmod.ionsjc_tip3p'                              >>tleap.script
echo "NEW_parm =  loadamberparams $lig_parm"                            >>tleap.script
echo "loadamberprep $lig_prep   "                                       >>tleap.script
echo "pd98 = loadpdb  $initial_pdb "                                    >>tleap.script
echo "saveamberparm pd98 $nom_top $nom_rst    "                         >>tleap.script
echo "savePdb pd98 $nom_pdb_leap             "                          >>tleap.script
echo 'quit'                                                             >>tleap.script

################
#file Runtleap
################

echo  "#"!"/bin/csh                              "                     >>Run_tleap
echo ' source $AMBERHOME/amber.sh                '                     >>Run_tleap
echo ' $AMBERHOME/AmberTools/bin/tleap -f   tleap.script '             >>Run_tleap

chmod u+x Run_tleap
./Run_tleap

cp $nom_top ../
cd ..

rm fftop
mkdir fftop

cd ./fftop

#############
#anteMMPBSA
#############

echo '   '$AMBERHOME'/bin/ante-MMPBSA.py     \'  >> run_anteMMPBSA
echo "   -p ../$nom_top                      \"  >> run_anteMMPBSA
echo "   -c ./complex_com.top                \"  >> run_anteMMPBSA
echo "   -r ./complex_rec.top                \"  >> run_anteMMPBSA
echo "   -l ./complex_lig.top                \"  >> run_anteMMPBSA
echo '   -m  '"'":$ini_res-$ifi_res"'" "     \"  >> run_anteMMPBSA
echo '   -s  ':WAT,Na+,Cl-' \'                   >> run_anteMMPBSA
echo "   --radii=mbondi2    >> ante-MMPBSA.log    "  >> run_anteMMPBSA

   chmod u+x run_anteMMPBSA
 
  ./run_anteMMPBSA

cd ../

foreach file (*.nc)
  echo $file
  set pose = ${file:r}
  echo $pose
  rm -f -r $pose
  mkdir $pose
  cd $pose
  mkdir ffdir
  mkdir fftop
  cp -r ../fftop ./
  cd ./ffdir
  
#................
# Inputs MMPBSA .
#................

  rm  mmpbsa_py.inp

  echo '&general                                                  ' >> mmpbsa_py.inp
  echo " startframe=$startf, endframe=$endf, interval=$nfreq      " >> mmpbsa_py.inp
  echo ' keep_files=0,                                            ' >> mmpbsa_py.inp
  echo '/                                                         ' >> mmpbsa_py.inp
  if $pb == 'YES' then
    echo '&pb                                                ' >> mmpbsa_py.inp
    echo ' cavity_surften = 0.005420, cavity_offset = 0.9200,' >> mmpbsa_py.inp
    echo " indi = $indi, exdi = $exdi,                       " >> mmpbsa_py.inp
    echo ' istrng=0.0, radiopt = 0, inp = 1,                 ' >> mmpbsa_py.inp
    echo '/                                                  ' >> mmpbsa_py.inp
  endif
  echo '&gb                                                ' >> mmpbsa_py.inp
  echo ' igb=2,                                            ' >> mmpbsa_py.inp
  echo ' surften = 0.0072,                                 ' >> mmpbsa_py.inp
  echo ' saltcon=0.0                                       ' >> mmpbsa_py.inp
  echo '/                                                  ' >> mmpbsa_py.inp

#.............
# Run MMGBSA . 
#.............
# This is a file to run the MMGBSA calculation : 
# MUST be changed for each Computational Center
#
  
  rm  send_mmpbsa

  echo "#"!"/bin/bash                              "                >> send_mmpbsa           
  echo "########################################## "                >> send_mmpbsa
  echo "# Opcions i parametres de l'SGE            "                >> send_mmpbsa
  echo "########################################## "                >> send_mmpbsa
  echo "# (1) Name of the Work ( Indentification ) "                >> send_mmpbsa
  echo "#$ -N $pose                                "                >> send_mmpbsa
  echo "# (2) Computational Resources (not needed) "                >> send_mmpbsa
  echo "##$ -l h_rt                                "                >> send_mmpbsa
  echo "##$ -l mem_free                            "                >> send_mmpbsa
  echo "#$ -pe smp $num_proc                       "                >> send_mmpbsa
  echo "##$ -l exclusive=true                      "                >> send_mmpbsa
  echo "# (2) Output information                   "                >> send_mmpbsa
  echo "#$ -cwd                                    "                >> send_mmpbsa
  echo "#$ -o $pose.out          		   "                >> send_mmpbsa
  echo "#$ -e $pose.err 		           "                >> send_mmpbsa
  echo "# (4) Send e-mail at the end of the work   "                >> send_mmpbsa
  echo "##$ -m e                                   "                >> send_mmpbsa
  echo "##$ -M  Name_User@ub.edu                   "                >> send_mmpbsa
  echo "########################################## "                >> send_mmpbsa
  echo "# Software resources                       "                >> send_mmpbsa
  echo "########################################## "                >> send_mmpbsa
  echo "# Modules                                  "                >> send_mmpbsa
  echo "module load numpy/1.6.1                    "                >> send_mmpbsa
  echo "module load python/2.7.10                  "                >> send_mmpbsa
  echo "module load ambertools/16_ompi             "                >> send_mmpbsa
  echo "########################################## "                >> send_mmpbsa
  echo 'echo " SGE_O_WORKDIR : $SGE_O_WORKDIR"     '                >> send_mmpbsa
  echo 'echo "nslots  : $NSLOTS "                  '                >> send_mmpbsa
  echo 'echo "TMP DIR   : $TMPDIR"                 '                >> send_mmpbsa
  echo "########################################## "                >> send_mmpbsa
  echo "# calcul                                   "                >> send_mmpbsa
  echo "########################################## "                >> send_mmpbsa
  echo "set RUN=/home/g6jaime/jaime/run            "                >> send_mmpbsa
  echo 'echo "RUN DIR   : $RUN "                   '                >> send_mmpbsa
  echo 'echo " AMBERHOME   : $AMBERHOME "          '                >> send_mmpbsa
  echo 'source "$AMBERHOME/amber.csh"              '                >> send_mmpbsa
  echo "#                                          "                >> send_mmpbsa
  echo 'cd ../fftop                                '                >> send_mmpbsa
  echo 'set TOP_DIR=`pwd`                          '                >> send_mmpbsa
  echo 'cd $SGE_O_WORKDIR                          '                >> send_mmpbsa
  echo 'cd ../../                                  '                >> send_mmpbsa
  echo 'set DIN_DIR=`pwd`                          '                >> send_mmpbsa
  echo 'cd $SGE_O_WORKDIR                          '                >> send_mmpbsa
  echo 'cp $SGE_O_WORKDIR/mmpbsa_py.inp $TMPDIR    '                >> send_mmpbsa
  echo 'cd $TMPDIR                                 '                >> send_mmpbsa

  echo " mpirun  -np $num_proc  MMPBSA.py.MPI -O   \"               >> send_mmpbsa
  echo '              -i  ./mmpbsa_py.inp \'                        >> send_mmpbsa
  echo '              -o  ./mmpbsa_py.out \'                        >> send_mmpbsa
  echo '              -eo ./DeltaG_BySnaps.csv \'                   >> send_mmpbsa
  echo '              -sp  $DIN_DIR/'$nom_top"     \"               >> send_mmpbsa
  echo '              -cp  $TOP_DIR/complex_com.top   \'            >> send_mmpbsa
  echo '              -rp  $TOP_DIR/complex_rec.top   \'            >> send_mmpbsa
  echo '              -lp  $TOP_DIR/complex_lig.top   \'            >> send_mmpbsa
  echo '              -y   $DIN_DIR/'$file'      \'                 >> send_mmpbsa
  echo "              >>  ./progres.log               "             >> send_mmpbsa
  echo '#'                                                          >> send_mmpbsa
  echo 'rm _M* refer*        '                                      >> send_mmpbsa
  echo 'cp * $SGE_O_WORKDIR  '                                      >> send_mmpbsa



  chmod u+x send_mmpbsa

  echo "cd $SGE_O_WORKDIR  "                             >> ../../send_mmpbsa_all
  echo "cd  $pose/ffdir    "                             >> ../../send_mmpbsa_all
  echo "qsub -q $nom_computer"'.q send_mmpbsa '          >> ../../send_mmpbsa_all


  cd ../../

end

chmod u+x send_mmpbsa_all
rm -r fftop
 
