#!/usr/bin/env python3
import numpy as np
import os, sys, stat
import time
from math import sqrt
import  matplotlib.pyplot as plt
from argparse import ArgumentParser

""" __Input Data__  """
"""
1. use_pocket: True/False
        : True : The reference structure is defined by residues in a pocket 
2. file_pocket : One line with the residues of the protein defining the pocket
                 ( in AMBER Mask format ) 
3. XRay : False/ True 
        : False :  The reference structure is not the X-Ray Structure
        : True  :  The reference structure is     the X-Ray Structure

# If use_pocket = False AND XRay= False, the LAST structure of each trajectory is
  used as reference for the plots time/distance. Where distance is defined between 
  the position of C99 in the last point and the position of C99 along the trajectory. 

4. protein_pdb : Must be defined and ALL the trajectories must be superposed to this structure.
                 The protein structure to be used to analyse reactivity.
                 Can include the ligand or not. If XRay=True The Ligand must be present. 

5. dis_lig_prot_min : A trajectory in which all the ligand atoms are farther than 
                      this value with respect to any protein atoms is a NON REACTIVE trajectory.

6. num_plots    : Number of different plots to be done.

7. dis_plots    : Limits for the axis Y .

8. snaps_by_one : Number of snapshots per ns .

9. num_ns_anal  : Number of ns to be analysed in depth.  

10. percent_anal : % of points in the analysed ns ( num_ns_anal ) that are allowed to be
                  greater than " dis_lig_prot_min " to stay the trajectory as REACTIVE.

11. prefix_filpdb : Prefix for the files containing the positions of the atoms with wdW repulsion.

12. name_atomrep  : Name of the atom with wdW repulsion.
        
"""
###########MAIN
def main() :

  prog_info ()
  
  args = cmdlineparse()

  cha_pocket =      args.use_pocket 
  cha_XRay   =      args.XRay        
  file_pocket=      args.file_pocket 
  protein_pdb=      args.protein_pdb 

  XRay = False
  if cha_XRay   == 'True' :
    XRay = True

  use_pocket = False
  if cha_pocket == 'True' :
    use_pocket = True

  dis_lig_prot_min = float ( args.dis_lig_prot_min )
  num_plots = int ( args.num_plots )
  dis_plots =     ( args.dis_plots )
  snaps_by_one = int ( args.snaps_by_one )
  num_ns_anal  = int ( args.num_ns_anal  )
  percent_anal = int ( args.percent_anal )

  prefix_filrep=     ( args.prefix_filrep )
  name_atomrep =     ( args.name_atomrep  )

  print ( " Atom name to analyse            : ",name_atomrep )

  print ( " Use pocket to define distances  : ",use_pocket )
  if use_pocket :
     print ( " Pocket Definition            : ",file_pocket )
  print ( " Use XRay structure as reference : ",XRay )
  if XRay :
     print ( " Using X-Ray structure           : ",protein_pdb )
  else :
     print ( " Reference Protein structure     : ",protein_pdb )
  print ( " Minimum distance Ligand-Prot    : ",dis_lig_prot_min )
  print ( " Num Plots to be done            : ",num_plots )       
  print ( " Dis Plots  [x,y,..]             : ",dis_plots )       

  print ( " Snaps by nanosecond             : ",snaps_by_one )    
  print ( " Num ns to analyse               : ",num_ns_anal  )    
  print ( " PerCent out of the limits       : ",percent_anal )    

  print ( " Prefix added to pdb files       : ",prefix_filrep )

  dis_plots = dis_plots.split(",")
  for ndim in range (len(dis_plots)) :
    dis_plots[ndim] = float (dis_plots[ndim] )

  initial_dir = os.getcwd()
  print ( "  " )
  print ( " Initial Dir = ",initial_dir )
  print ( "  " )

# Directory for pictures

  dir_fig_dist = 'figures_distances'
  exist_dir = os.path.isdir(dir_fig_dist)
  if exist_dir :
    if os.listdir(dir_fig_dist) != [] :
      os.chdir(dir_fig_dist)
      os.system ( "rm * " )
      os.chdir('..')
    os.rmdir(dir_fig_dist)
  os.mkdir(dir_fig_dist)
#

  num_reac = 0

  num_at_prot, num_at_lig, pos_name_atomrep, resi_name_lig = read_prot (protein_pdb,name_atomrep)
  print ( "... There are = ",num_at_prot," PROTEIN Atoms " )
  if num_at_lig != 0 :
    print ( "... There are = ",num_at_lig, " LIGAND  Atoms " )
    print ( "... With Residue name = ",resi_name_lig )
    print ( "... And ",name_atomrep, " at position ",pos_name_atomrep+1 )

#file_delete_nc = ("Pre_Grafic_{}.out".format(nombres[KUAL]))

  file_delete_nom = ("delete_Non_Reactive")
  file_delete = open (file_delete_nom,"w")

  if use_pocket :
    print ( "... Pocket ACTIVE for reference definition " ) 
    if XRay :
      print ( " ... Pocket ACTIVE thus XRay turned OFF " )
      XRay = False
    res_pocket  = read_pocket ( file_pocket ) 
    x_ref,y_ref,z_ref = calc_cent_geom (res_pocket,num_at_prot,x_prot,y_prot,z_prot,at_name,res_num)
    print ( "    ... Pocket REF  {:.3f} {:.3f} {:.3f} ".format(x_ref,y_ref,z_ref)) 

  if num_at_lig != 0 and XRay : 
    x_ref = x_lig[pos_name_atomrep]
    y_ref = y_lig[pos_name_atomrep]
    z_ref = z_lig[pos_name_atomrep]
    print ( " ... XRay REF  {:.3f} {:.3f} {:.3f} ".format(x_ref,y_ref,z_ref)) 

  print_ref = True
  lista = os.listdir( os.getcwd() )
  for pdb in lista:
    if prefix_filrep in pdb:
      print("                     ")
      print("... Pocessing File  :",pdb)

      pdb_split = pdb.split("_")
      lig_num = pdb_split [1]

      x  = []
      y  = []
      z  = []
      num_at = read_pdb_c99(pdb,name_atomrep,x,y,z)

      num_ns = num_at / snaps_by_one
      num_snaps_anal = num_at - snaps_by_one * num_ns_anal

      if num_at_lig == 0 or ( not XRay ) and ( not use_pocket)  : 
        x_ref = x[-1]
        y_ref = y[-1]
        z_ref = z[-1]
        if   print_ref  :
          print ( "    ... SEL  REF  {:.3f} {:.3f} {:.3f} ".format(x_ref,y_ref,z_ref)) 
          print_ref = False

      min_lig_prot, num_near_ref = dis_min_prot (num_at_prot,dis_lig_prot_min,x_ref,y_ref,z_ref)

      dis_at = []
      d_max = 0.0
      for i in range( num_at ) :
        x_tm = x [i]
        y_tm = y [i]
        z_tm = z [i]
        x2 = ( x_tm - x_ref ) * ( x_tm - x_ref )
        y2 = ( y_tm - y_ref ) * ( y_tm - y_ref )
        z2 = ( z_tm - z_ref ) * ( z_tm - z_ref )
        d = sqrt ( x2 + y2 + z2 )
        if d > d_max :
          d_max = d
        dis_at.append(d)

      print ( "    ... D_max ( {}_last --> {}_dyn  )    = {:.4f} ".format(name_atomrep,name_atomrep,d_max) )
      print ( "    ... D_min ( {}_last --> Protein )    = {:.4f} ".format(name_atomrep,min_lig_prot) )
      print ( "    ... Number of PROTEIN atoms near REF = {:6d}  ".format(num_near_ref) )

      if min_lig_prot >= dis_lig_prot_min :

        print("    --> File :",pdb, ' Is a NON reactive trajectory withing ',dis_lig_prot_min, ' A ')

        name_to_delete     = 'lig_' + lig_num + '*'
        name_nc_to_delete  = 'lig_' + lig_num + '.nc' 
        name_pdb_to_delete = 'lig_' + lig_num + prefix_filrep + '.pdb' 
        name_end_to_delete = 'lig_' + lig_num + '_LAST.pdb'
        name_csv_to_delete = 'lig_' + lig_num + '_disat.csv'

        print("       --> Files To Delete :",name_to_delete)

        name_nc_to_delete  = ' rm -f ' + name_nc_to_delete  + '\n'
        name_pdb_to_delete = ' rm -f ' + name_pdb_to_delete + '\n'
        name_end_to_delete = ' rm -f ' + name_end_to_delete + '\n'
        name_csv_to_delete = ' rm -f ' + name_csv_to_delete + '\n'

        file_delete.write( name_nc_to_delete  )
        file_delete.write( name_pdb_to_delete )
        file_delete.write( name_end_to_delete )
        file_delete.write( name_csv_to_delete )

      else :
        num_no_react   = 0
        num_snaps_work =  num_at - num_snaps_anal
        for i in range ( num_snaps_anal,num_at ) :
          x_tm = x [i]
          y_tm = y [i]
          z_tm = z [i]
          min_snap_prot,num_near_atom = dis_min_prot (num_at_prot,dis_lig_prot_min,x_tm,y_tm,z_tm)
          if min_snap_prot >= dis_lig_prot_min :
            num_no_react  += 1
        print ( "    ... There are ",num_no_react," Farther than ",dis_lig_prot_min," in a total of ",num_snaps_work )

        num_no_react   = 0
        for i in range ( num_snaps_anal,num_at ) :
          if dis_at[i] >= dis_lig_prot_min :
            num_no_react  += 1
        print ( "    ... Non_Reactive snaps ",num_no_react," in a total of ",num_snaps_work )
        per_cent = 100 * num_no_react / num_snaps_work

        if per_cent  >= percent_anal :
          print("    --> File :",pdb, ' Is a NON-REACTIVE reactive trajectory with ',per_cent, ' % OUT')

          name_to_delete     = 'lig_' + lig_num + '*'
          name_nc_to_delete  = 'lig_' + lig_num + '.nc' 
          name_pdb_to_delete = 'lig_' + lig_num + prefix_filrep + '.pdb' 
          name_end_to_delete = 'lig_' + lig_num + '_LAST.pdb'
          name_csv_to_delete = 'lig_' + lig_num + '_disat.csv'

          print("       --> Files To Delete :",name_to_delete)

          name_nc_to_delete  = ' rm -f ' + name_nc_to_delete  + '\n'
          name_pdb_to_delete = ' rm -f ' + name_pdb_to_delete + '\n'
          name_end_to_delete = ' rm -f ' + name_end_to_delete + '\n'
          name_csv_to_delete = ' rm -f ' + name_csv_to_delete + '\n'

          file_delete.write( name_nc_to_delete  )
          file_delete.write( name_pdb_to_delete )
          file_delete.write( name_end_to_delete )
          file_delete.write( name_csv_to_delete )
        else :
          print("    --> File :",pdb, ' Is a REACTIVE trajectory withing ',dis_lig_prot_min, ' A ')
          num_reac += 1
         

      x_din = np.linspace( 0, num_ns, num_at )
      y_dis = np.array ( dis_at )
    
      print ( "    ... Analysing ",num_ns_anal," ns with ",snaps_by_one * num_ns_anal," snaps in a total of ",num_at)
      y_ave, sd_y, y_max   = stat_x_ns (num_at,num_snaps_anal,y_dis )
      print ( "        ... y_ave ({:.3f} ),sd( {:.3f} ) y_max( {:.3f} )" .format(y_ave, sd_y, y_max )) 

      file_excel_disat = 'lig_'+lig_num+'_disat.csv'
      file_excel       = open (file_excel_disat,"w")

      for rowx in range ( num_at ) :
        strw = ( "{:15.9} {:15.9}".format(x_din[rowx],y_dis[rowx]) )
        strw = str ( strw ) 
        strw  = strw + '\n'
        file_excel.write(strw)

      file_plot        = 'lig_'+lig_num
    
      for num_plt in range ( num_plots ) :

        file_plot_disat  = 'lig_'+lig_num+'_disat'+str(dis_plots[num_plt])+'.png'
        distan_to_plot = dis_plots[num_plt]

        #"""
        plt.xlim(0,num_ns)
        plt.ylim(0,distan_to_plot)

        plt.title(file_plot)
        plt.xlabel("t(ns)")
        plt.ylabel("Distance (Angstrom)")
        plt.plot (x_din,y_dis,'g', label='Dist_Last')
        plt.legend()
        plt.savefig (file_plot_disat, format='png',dpi=600)
        plt.show()

        move = 'mv ' + file_plot_disat + ' ./' + dir_fig_dist + '/'
        os.system ( move  )
        #""" 

  os.chmod (file_delete_nom,0o0764)

  print ( "..." ) 
  print ( "... There are = ",num_reac," REACTIVE trajectories " ) 
  print ( "..." ) 

  exit ()

""" read_pocket """

def read_pocket (data_inp):
  residues = []
  file_pocket =  open (data_inp)
  num_lines = 0
  for lines in file_pocket :
    num_lines += 1
  if num_lines != 1 :
    print ( " Only ONE line is allowed : STOP ")
    exit ()
  file_pocket =  open (data_inp)
  line = file_pocket.readline ().replace("\n","")
  num_char = len (line)
# print ( " Num_Char = ",num_char )
  double_point = line.find(":")
# print ( "  :       = ",double_point )
  if double_point == -1 :
    print("Please use : ")
    exit ()
  line_nop = line[double_point+1:num_char]
# print( line_nop )
  coma = line.find(",")
  if coma == -1 :
    hairline = line.find("-")
    if hairline == -1 :
      residues.append ( int ( line_nop ) )
    else:
      range_res = line_nop.split("-")
      for i in range(len(range_res)):
        range_res[i] = int ( range_res[i] )
#     print ( range_res )
      for i in range(range_res[0],range_res[1]+1):
        residues.append ( i )
    print ( residues )
  else:
    all_range_res = line_nop.split(",")
#   print ( all_range_res )
    for i in range ( len(all_range_res) ) :
      if len (all_range_res[i]) == 0 :
        print ( " Error in Residue definition ")
        exit ()
    for i in range ( len(all_range_res) ) :
      range_res = all_range_res[i].split("-")
      for i in range(len(range_res)):
        range_res[i] = int ( range_res[i] )
      if len (range_res) == 2 :
        for j in range(range_res[0],range_res[1]+1):
          residues.append ( j )
      else:
          residues.append ( range_res[0] )

# print ( residues )
  return residues

""" Geometry Center for the pocket """

def calc_cent_geom ( res_pocket,num_at_prot,x_prot,y_prot,z_prot,at_name,res_num) :
  num_res_pocket = len(res_pocket)
# print ( "**",num_res_pocket)
  x_ref = 0.0
  y_ref = 0.0
  z_ref = 0.0
  for i in range(num_res_pocket):
    for n_at in range (num_at_prot):
      if  res_num[n_at] == res_pocket[i] and at_name[n_at] == 'CA' :
#       print (  res_num[n_at],res_pocket[i],at_name[n_at],x_prot[n_at] )
        x_ref += x_prot[n_at]
        y_ref += y_prot[n_at]
        z_ref += z_prot[n_at]
  x_ref /= num_res_pocket
  y_ref /= num_res_pocket
  z_ref /= num_res_pocket
  return x_ref,y_ref,z_ref

""" stat_x_ns  """     

def stat_x_ns (num_at,num_snaps_anal,y_dis ) :
  sy    =  0.0
  sd_y  =  0.0
  y_max =  0.0
  num_snaps_work =  num_at - num_snaps_anal
  for ns2anal in range ( num_snaps_anal,num_at ) :
    y   = y_dis [ns2anal] 
    sy  = y + sy
#   print ( ns2anal, y) 
    if y > y_max :
      y_max = y
  y_ave = sy / num_snaps_work
  for ns2anal in range ( num_snaps_anal,num_at ) :
    y  = y_dis[ns2anal] 
    sd_y = sd_y + ( y - y_ave ) * ( y - y_ave )
  sd_y = sqrt ( sd_y / num_snaps_work )
# print ( " y_ave ({:.3f} ),sd( {:.3f} ) y_max( {:.3f} )" .format(y_ave, sd_y, y_max )) 
  return y_ave, sd_y, y_max

""" dis_min_prot """

def dis_min_prot (num_at_prot,dis_lig_prot_min,x_cal,y_cal,z_cal):
  d_min = 100 * dis_lig_prot_min  
  num_nearest = 0
  for i in range( num_at_prot ) :
    x_tm = x_prot [i]
    y_tm = y_prot [i]
    z_tm = z_prot [i]
    x2 = ( x_tm - x_cal ) * ( x_tm - x_cal )
    y2 = ( y_tm - y_cal ) * ( y_tm - y_cal )
    z2 = ( z_tm - z_cal ) * ( z_tm - z_cal )
    d = sqrt ( x2 + y2 + z2 )
    if d < d_min :
      d_min = d
    if d <= dis_lig_prot_min :
      num_nearest += 1
# print ( "    ... D_min = {:.4f} ".format(d_min) )
  return d_min,num_nearest

""" read_prot """

def read_prot (pdb,name_atomrep):
  file_pdb  = open (pdb)
  num_lines = 0
  for lines in file_pdb :
    num_lines += 1
# print ( "    ... There are = ",num_lines," Lines " )
  file_pdb  = open (pdb)
  num_atoms     = 0
  num_atoms_lig = 0
  protein_on = True 
  resi_name_lig   = 'UNK'

  for i in range ( num_lines )  :

    line = file_pdb.readline ().replace("\n","") 
    if 'ATOM' in line and protein_on :
      values=line.split()
      at_name.append (values [2])
      res_num.append (int (values [4]))
      x_prot.append (float(values [5]))
      y_prot.append (float(values [6]))
      z_prot.append (float(values [7]))
      num_atoms += 1
    elif  num_atoms != 0 :
      protein_on = False

    if 'ATOM' in line and (not protein_on) :
      values=line.split()
#     print ( values )
      atom_name_lig.append (values [2])
      resi_name_lig   = values [3]
      x_lig.append (float(values [5]))
      y_lig.append (float(values [6]))
      z_lig.append (float(values [7]))
      num_atoms_lig += 1

  if  num_atoms_lig != 0 :
    for j in range ( num_atoms_lig )  :
      if atom_name_lig[j] == name_atomrep :      
        pos_name_atomrep = j
  else :
    pos_name_atomrep = 0

  return num_atoms,num_atoms_lig,pos_name_atomrep,resi_name_lig
"""
  print ( "... ... There are = ",num_atoms," Atoms Protein " )
  for i in range ( num_atoms )  :
    print ( " ATOM {:.3f} {:.3f} {:.3f} ".format(x_prot[i],y_prot[i],z_prot[i])) 
  if num_atoms_lig != 0 :
    print ( "... ... There are = ",num_atoms_lig," Atoms Ligand " )
    for i in range ( num_atoms_lig )  :
      print ( " ATOM {:.3f} {:.3f} {:.3f} ".format(x_lig[i],y_lig[i],z_lig[i])) 
# exit ()
"""

""" read_pdb_c99 """

def read_pdb_c99 (pdb,name_atomrep,x,y,z):
  file_pdb  = open (pdb)
  num_lines = 0
  for lines in file_pdb :
    num_lines += 1
# print ( "... ... There are = ",num_lines," Lines " )
  file_pdb  = open (pdb)
  num_atoms = 0
  for i in range ( num_lines )  :
    line = file_pdb.readline ().replace("\n","") 
    if name_atomrep in line :
      values=line.split()
      x.append (float(values [5]))
      y.append (float(values [6]))
      z.append (float(values [7]))
      num_atoms += 1
  print ( "    ... There are = ",num_atoms," Atoms(Structures) " )
  return num_atoms

""" END read_pdb_c99 """

# ... READ datafiles

def cmdlineparse():

    parser = ArgumentParser(description="command line arguments")

    parser.add_argument("-pocket",     dest="use_pocket"   , required=True,   help=" Use pocket to define distances (False) ")
    parser.add_argument("-fil_pocket", dest="file_pocket"  , required=False,  help=" File containing the aa of the pocket ")

    parser.add_argument("-xray"    , dest="XRay"         , required=True,  help=" Use XRay structure as reference (False) ")
    parser.add_argument("-fil_prot", dest="protein_pdb"  , required=True,  help=" Reference Protein structure ")

    parser.add_argument("-dis_min" , dest="dis_lig_prot_min" , required=True, help=" Minimum distance Ligand-Prot ")

    parser.add_argument("-prefix_filpdb" , dest="prefix_filrep"  , required=True, help=" Posfix added to pdb files ")
    parser.add_argument("-name_atom"     , dest="name_atomrep"  , required=True, help=" Atom to analyse (C99) ") 

    parser.add_argument("-num_plots" , dest="num_plots"  , required=True, help=" Num Plots to be done (2)  ")
    parser.add_argument("-dis_plots" , dest="dis_plots"  , required=True, help=" Distances to be Plot  [5,25] ")

    parser.add_argument("-snaps_byone"  , dest="snaps_by_one"   , required=True, help=" Snaps by nanosecond ")
    parser.add_argument("-num_ns_anal"  , dest="num_ns_anal"    , required=True, help=" Num ns to analyse   ")
    parser.add_argument("-percent_anal" , dest="percent_anal"   , required=True, help=" PerCent out of the limits ")

    args=parser.parse_args()
    return args

def prog_info () :
  print ( "... ... ... ... ... ... ... " )
  print ( "...  fdMD_ReactiveTraj  ... " )
  print ( "...     ( 2019 )        ... " )
  print ( "... ... ... ... ... ... ... " )
  print ( " " )

if __name__ == '__main__':

  x_prot   = []
  y_prot   = []
  z_prot   = []
  at_name  = []
  res_name = []
  res_num  = []
  x_lig    = []
  y_lig    = []
  z_lig    = []
  atom_name_lig = []
  resi_name_lig = ' '

  main()

