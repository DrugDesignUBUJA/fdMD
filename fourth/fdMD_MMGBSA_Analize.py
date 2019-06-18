#!/usr/bin/env python3

"""
 .......................
  fdMD_MMGBSA_Analize.py
 .......................

The structure of directories MUST be :
1) Dir where you run this file
- prefix_Dir_1
- prefix_Dir_2
- .....
Into these directories you MUST have
---- ffdir
---- fftop
"""

import numpy as np
import os
import time
from math import sqrt
import  matplotlib.pyplot as plt

snaps_by_one = 1      # Number of snapshot by nanosecond 
ns_ave       = 10     # Number of ns for the average
prefix       ='lig_'  # Prefix 

gb_mmpbsa = "GENERALIZED BORN:"
pb_mmpbsa = "POISSON BOLTZMANN:"
delta     = "DELTA Energy Terms"

def ave_data (energy) :
  is_div = num_ns%ns_ave
  if is_div == 0 :
    num_interval = int ( num_ns // ns_ave )
    print ( "num_interval",num_interval)
  else :
    num_interval = 1

  ave_interval = []
  num_ns_interval = int ( num_ns / num_interval )
  data_interval    = num_ns_interval * snaps_by_one
  ind = 0
  for interval in range ( num_interval ) :
#   print ( " >>> interval = ",interval )
    ave = 0.0
    for data in range ( data_interval ) :
#     print ( " >>> >>> data = ",data )
      ave = ave + energy [ind]
#     print ( " >>> >>> >>> ind = ",ind )
      ind += 1
    ave /= data_interval
    ave_interval.append ( ave ) 
  print ( ave_interval )

  ave_graph = []
  for interval in range ( num_interval ) :
    for data in range ( data_interval ) :
      ave_graph.append ( ave_interval [interval] )

  num_data_new = len ( ave_graph )
  print ( "num_data_new = ",num_data_new )

  return ave_interval,ave_graph


initial_dir = os.getcwd()
print ( " Working in the Directory : ",initial_dir )

exist_dir = os.path.isdir('pictures')
if exist_dir : 
  if os.listdir('pictures') != [] :
    os.chdir('pictures')
    os.system ( "rm * " )
    os.chdir('..')       
  os.rmdir('pictures')
os.mkdir('pictures')

dir_all = os.listdir(os.getcwd())
num_dir = 0
for dir in dir_all:
  if os.path.isdir(dir):
    if prefix in dir :
      num_dir += 1
      print(" {} {}".format(num_dir,dir))

print ( " There are {} Directories : ".format(num_dir)) 

file_ave  = open ( "ave_gbpb.res","w" )

for dir in dir_all:
  if os.path.isdir(dir):
    if prefix in dir :
      name_dir = dir.split("_")
      molec =  name_dir  [0]
      pose  =  name_dir  [1]
      line_to_write = str (molec) + '_' + str (pose) 
      print ( molec,pose )
      os.chdir(dir)
      print ( " >>> Working in the Directory : ",dir )
#     dir_mmpbsa = ("mmpbsa_{}".format(dir))
#     print ( " >>> >>> mmpbsa_ Directory : ",dir_mmpbsa )
#     os.chdir(dir_mmpbsa)
      dir_ffdir_all = os.listdir(os.getcwd())
      for dir_ffdir in dir_ffdir_all:
        if os.path.isdir(dir_ffdir):
          if 'ffdir' in dir_ffdir :
            os.chdir(dir_ffdir)
            print ( " >>> >>> >>> ffdir_ Directory : ",dir_ffdir )
       
            os.system("cp DeltaG_BySnaps.csv ./tempfile.csv")
            os.system("ls -l")
#
# 1. GENERALIZED BORN
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            gb_data  = False
            gb_delta = False
            lin_gb = 0
            for line in file_mmpbsa :
              lin_gb += 1
              if gb_mmpbsa in line :
                gb_data = True
                gp_pos  = lin_gb
                print ( " >>> >>> >>> >>> GENERALIZED BORN " )
              if gb_data and not gb_delta:
                if delta in line :
                  gb_delta = True
                  gb_delta_pos  = lin_gb
                  print ( " >>> >>> >>> >>> >>> GB Delta ",gb_delta_pos )
                  break
#
# 1. POISSON BOLTZMANN 
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            pb_data  = False
            pb_delta = False
            lin_pb = 0
            for line in file_mmpbsa :
              lin_pb += 1
              if pb_mmpbsa in line :
                pb_data = True
                pp_pos  = lin_pb
                print ( " >>> >>> >>> >>> POISSON BOLTZMANN " )
              if pb_data and not pb_delta:
                if delta in line :
                  pb_delta = True
                  pb_delta_pos  = lin_pb
                  print ( " >>> >>> >>> >>> >>> GB Delta ",pb_delta_pos )
                  break
#
# 2. GENERALIZED BORN
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            if gb_delta :
              Total_Energy_GB = []
              for lin_null in range ( 1,gb_delta_pos+1 ) :
                line = file_mmpbsa.readline ().replace("\n","")
#               print ( " line ",line )
              line_info = file_mmpbsa.readline ().replace("\n","")
              names_ener =  line_info.split(",")
              num_term = len (names_ener)
              for num_t in range(num_term):
#               print ( names_ener[num_t])
                if  names_ener[num_t] == "DELTA TOTAL" :
                  ind_ETot_GB = num_t
#             print ( " Total Energy = ", names_ener[ind_ETot_GB])
              num_ener = 1
              line_ener = file_mmpbsa.readline ().replace("\n","")
              while line_ener != "" :
                energies = line_ener.split(",")
                Total_Energy_GB.append( float (energies [ind_ETot_GB] ))
                line_ener = file_mmpbsa.readline ().replace("\n","")
                num_ener += 1
              num_ener -= 1
              print ( " There are {} GB snapshots  ".format(num_ener) )
              num_ns = num_ener / snaps_by_one
              y_gb = np.array ( Total_Energy_GB )
              ave_gb, y_ave_gb = ave_data (y_gb)
              y_ave_gb = np.array ( y_ave_gb )
              line_to_write = line_to_write + ' GB '
              for i in range ( len(ave_gb) ) :
                line_to_write = line_to_write + str(ave_gb[i]) + ' '
#
# 2. POISSON BOLTZMANN
#
            file_mmpbsa = open ( "tempfile.csv","r+" )
            if pb_delta :
              Total_Energy_PB = []
              for lin_null in range ( 1,pb_delta_pos+1 ) :
                line = file_mmpbsa.readline ().replace("\n","")
#               print ( " line ",line )
              line_info = file_mmpbsa.readline ().replace("\n","")
              names_ener =  line_info.split(",")
              num_term = len (names_ener)
              for num_t in range(num_term):
#               print ( names_ener[num_t])
                if  names_ener[num_t] == "DELTA TOTAL" :
                  ind_ETot_PB = num_t
#             print ( " Total Energy = ", names_ener[ind_ETot_PB])
              num_ener = 1
              line_ener = file_mmpbsa.readline ().replace("\n","")
              while line_ener != "" :
                energies = line_ener.split(",")
                Total_Energy_PB.append( float (energies [ind_ETot_PB] ))
                line_ener = file_mmpbsa.readline ().replace("\n","")
                num_ener += 1
              num_ener -= 1
              print ( " There are {} PB snapshots  ".format(num_ener) )
              num_ns = num_ener / snaps_by_one
              y_pb = np.array ( Total_Energy_PB )
              ave_pb, y_ave_pb = ave_data (y_pb)
              y_ave_pb = np.array ( y_ave_pb )
              line_to_write = line_to_write + ' PB '
              for i in range ( len(ave_gb) ) :
                line_to_write = line_to_write + str(ave_pb[i]) + ' '

            line_to_write = line_to_write + '\n'
            file_ave.write(line_to_write)

            title  = str(molec) + '_' + str(pose) 
            file_plot = '../../pictures/' + dir + '.png'

            label_gb  = 'GB'
            label_pb  = 'PB'

            style_gb = '--'
            style_pb = ':'

            x_din = np.linspace( 0, num_ns, num_ener )

            plt.ylim()
            plt.xlim(0,num_ns )
            plt.xlabel("ns")
            plt.ylabel(r'$\Delta$G (kcal/mol)')

            plt.title(title)

            if  gb_delta and pb_delta :
              plt.plot (x_din,y_gb,color='b', label=label_gb)
              plt.plot (x_din,y_pb,color='r', label=label_pb)
              plt.plot (x_din,y_ave_gb,color='y',linestyle=style_gb,label=label_gb)
              plt.plot (x_din,y_ave_pb,color='y',linestyle=style_pb,label=label_pb)
            elif     gb_delta and not pb_delta :
              plt.plot (x_din,y_gb,'b', label=label_gb)
              plt.plot (x_din,y_ave_gb,color='y',linestyle=style_gb,label=label_gb)
            elif not gb_delta and     pb_delta :
              plt.plot (x_din,y_pb,'r', label=label_pb)
              plt.plot (x_din,y_ave_pb,color='y',linestyle=style_pb,label=label_pb)

            plt.legend()
            plt.savefig (file_plot, format='png',dpi=300)
            plt.show()

            os.system("rm ./tempfile.csv")
            os.system("ls -l")
            os.chdir("..")

#     os.chdir("..")
      os.chdir("..")
 
exit ()

