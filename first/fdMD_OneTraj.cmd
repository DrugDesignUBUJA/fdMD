#!/bin/csh -v 
###
# -----------------------
# -     fdMD_OneTraj    -
# -----------------------
#  
# AMBER version and directory
set amber_version = 18
set amber_dir = /home/prog/amber$amber_version
set cpptraj   = $amber_dir/bin/cpptraj

# Number of Protein residues
# Number of Ligands in the simulation

@ num_res_prot = 387
@ num_lig      =  26

# Residues to Superpose

@ ini_res_sup  =   1
@ ifi_res_sup  = 387

# Refers to the Protein PLUS ONE ligand: Without wat,Na+,Cl-.. but both prepared !!! 
# xray can be = yes/no

set xray = 'no' 

if ( $xray == 'yes' ) then
  set ref_pdb  = "Name_Of_XRay_PDB.pdb"
else
  set ref_pdb  = "first_lig.pdb"
endif

set dir_traj     = .
set top_name     = "2ohk_Sol_vdw2301_ZeroB.top"

set top_file     = $dir_traj/$top_name

set prefix_traj  = "2ohk"

@ ini_traj = 1  
@ end_traj = 10 

@ ini_read   =    1
set end_read = last
@ inc_read   =    1

set last_pdb = 'yes'
if ( $last_pdb == 'yes' ) then
  @ num_snaps_total  = 200
  echo " Number of Snapshots    : " $num_snaps_total
endif

set DIR="`pwd`"
echo ' Working Directory : ' $DIR

@ tot_traj = $end_traj - $ini_traj
@ tot_traj++
echo " Number of Trajectories : " $tot_traj

set wat         = '[wat]'
set nowat       = '[NoWat]'
set nowatOnelig = '[OneLigNoWat]'
set lig         = '[lig]'

rm  inptrj_one_lig

echo " $cpptraj > gen_one_traj.out << EOF     "  >> inptrj_one_lig
echo "                                        "  >> inptrj_one_lig
echo " parm  $top_file $wat                   "  >> inptrj_one_lig

@ num_traj = 0
while ($num_traj  < $tot_traj)
  @ count_traj = $num_traj + $ini_traj 
  set trajec = $prefix_traj'_'$count_traj'_dyn.nc'
  if ( $amber_version > 16 ) then
    echo " trajin $dir_traj/$trajec  $ini_read $end_read $inc_read parm $wat "  >> inptrj_one_lig
  else
    echo " trajin $dir_traj/$trajec  $ini_read $end_read $inc_read $wat "  >> inptrj_one_lig
  endif
  @ num_traj++
end
echo "                                        "  >> inptrj_one_lig
echo " autoimage                              "  >> inptrj_one_lig
echo " strip :WAT,Na+,Cl-  outprefix NoWat    "  >> inptrj_one_lig
echo " rms first out RMS_first.dat "':'"$ini_res_sup"'-'"$ifi_res_sup"'@CA ' >> inptrj_one_lig
echo " trajout RMSD_FIRST_NoWat_alone.nc      "  >> inptrj_one_lig
echo " run                                    "  >> inptrj_one_lig
echo "                                        "  >> inptrj_one_lig

#
#........
#

@ count_lig = $num_res_prot
@ count_lig++                  

@ anal_lig  = 1

while ($anal_lig  <= $num_lig)
  if ( $xray == 'yes' ) then
   set nc_out = first_lig_$count_lig.nc 
  else
   set nc_out =       lig_$count_lig.nc 
  endif
  echo " clear all                              "  >> inptrj_one_lig
  echo " parm  NoWat.$top_name $nowat           "  >> inptrj_one_lig
  if ( $amber_version > 16 ) then
     echo " trajin RMSD_FIRST_NoWat_alone.nc parm $nowat"  >> inptrj_one_lig
  else
     echo " trajin RMSD_FIRST_NoWat_alone.nc $nowat"  >> inptrj_one_lig
  endif
  if ( $anal_lig == 1 ) then
    echo " strip "'!'":1-$num_res_prot,$count_lig outprefix OneLig "  >> inptrj_one_lig
  else
    echo " strip "'!'":1-$num_res_prot,$count_lig                  "  >> inptrj_one_lig
  endif
  echo " trajout $nc_out  $nowat                "  >> inptrj_one_lig
  echo " run                                    "  >> inptrj_one_lig
  echo "                                        "  >> inptrj_one_lig

  if ( $xray == 'yes' ) then
    echo " clear all                                   "  >> inptrj_one_lig
    echo " parm  OneLig.NoWat.$top_name  $nowatOnelig  "  >> inptrj_one_lig
    if ( $amber_version > 16 ) then
      echo " trajin $nc_out         parm   $nowatOnelig  "  >> inptrj_one_lig
    else
      echo " trajin $nc_out                $nowatOnelig  "  >> inptrj_one_lig
    endif
    echo " reference "'./'"$ref_pdb      $nowatOnelig  "  >> inptrj_one_lig
    echo " rms reference  out RMS_XRay.dat "':'"$ini_res_sup"'-'"$ifi_res_sup"'@CA ' >> inptrj_one_lig
    echo " trajout lig_$count_lig.nc     $nowatOnelig  "  >> inptrj_one_lig
    echo " run                                         "  >> inptrj_one_lig
    echo "                                             "  >> inptrj_one_lig
  endif


  @ count_lig++                  
  @ anal_lig++                  
end

#
#........
#

@ count_lig = $num_res_prot
@ count_lig++                  

@ anal_lig  = 1

while ($anal_lig  <= $num_lig)
  echo " clear all                              "  >> inptrj_one_lig
  echo " parm  OneLig.NoWat.$top_name $lig      "  >> inptrj_one_lig
  if ( $amber_version > 16 ) then
    echo " trajin  lig_$count_lig.nc  parm $lig      "  >> inptrj_one_lig
  else
    echo " trajin  lig_$count_lig.nc       $lig      "  >> inptrj_one_lig
  endif
  echo " strip "'!'"@C99"                          >> inptrj_one_lig
  echo " trajout lig_$count_lig""_c99.pdb  pdb  "  >> inptrj_one_lig
  echo " run                                    "  >> inptrj_one_lig
  echo "                                        "  >> inptrj_one_lig
  if ( $last_pdb == 'yes' ) then
    echo " clear all                              "  >> inptrj_one_lig
    echo " parm  OneLig.NoWat.$top_name $lig      "  >> inptrj_one_lig
    if ( $amber_version > 16 ) then
      echo " trajin  lig_$count_lig.nc $num_snaps_total $num_snaps_total 1 parm $lig "  >> inptrj_one_lig
    else
      echo " trajin  lig_$count_lig.nc $num_snaps_total $num_snaps_total 1      $lig "  >> inptrj_one_lig
    endif
    echo " trajout lig_$count_lig"'_LAST.pdb pdb'     >> inptrj_one_lig
    echo " run                                    "  >> inptrj_one_lig
    echo "                                        "  >> inptrj_one_lig
  endif
  @ count_lig++                  
  @ anal_lig++                  
end


echo "EOF"                                         >> inptrj_one_lig
chmod u+x inptrj_one_lig
