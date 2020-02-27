import numpy as np
import sys
import MDAnalysis as mda
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d
import scipy.optimize as optimization

'''
methods of selection : 
range of residues : e.g. selecting 100 to 120 residues => 'protein and resid 100:120'
range of atoms : e.g. selecting 2102 to 2485 residues => 'protein and bynum 2102:2485'

'''
system=mda.Universe('step7_production.gro','step7_production.xtc')
output='output.xvg'

#selction of helicies  
helix1=system.select_atoms('protein and resid 237:265 and name CA')

#change the variable coordinates2 ='z' for measuring the angle with respect to Z axis.
coordinates2='x'
if coordinates2=='Z' or coordinates2 =='z' :
   helix2='z'
else :
   helix2=system.select_atoms('protein')

def Range(data):
   x=[]
   y=[]
   for i in range(0, len(data.T)):
       x.append(np.mean(data.T[i])- np.std(data.T[i]))
       y.append(np.mean(data.T[i]) + np.std(data.T[i]))
   normx,normy=np.linalg.norm(x),np.linalg.norm(y)
   minimum=normx-(normx+normy)/2
   maximum=normy-(normx+normy)/2
   return(minimum,maximum)

# Function to perform SVD and extract the eigenvector
def average_axis(data):
   minimum,maximum=Range(data)
   datamean = data.mean(axis=0)
   uu, dd, vv = np.linalg.svd(data - datamean)
   linepts = vv[0] * np.mgrid[minimum:maximum:2j][:, np.newaxis]
   linepts += datamean
   return(linepts)

# function that measure angle using the dot product between two vectors
def calc_angle(first_vector, second_vector,first1,last1, first2, last2):
         angle=np.degrees(np.arccos(abs(np.dot((first_vector[0] - first_vector[1]) , (second_vector[0]- second_vector[1])))/(np.linalg.norm(first_vector[0] - first_vector[1])*(np.linalg.norm(second_vector[0] - second_vector[1])))))
         #angle=180-angle
         return(angle)

Z_angle=[]

#  Loop to iterate through trajectory
for ts in system.trajectory[::]: ## [0:10:] first 10 frames will be considered, Pythonthic in nature. 
   sys.stdout.write("Complete \r%d%% " %((ts.frame*100)/len(system.trajectory)))
   sys.stdout.flush()
   data1=helix1.positions
   print(ts)
   if ts.frame==1:
      timestep=ts.time
   if len(helix2)==1:
      data2=np.array([[0,0,-1],[0,0,-10],[0,0,-20], [0,0,-30]])
      linepts2=np.array([[0,0,10],[0,0,20]])
   else:
      data2=helix2.positions
      linepts2=average_axis(data2)
   linepts1=average_axis(data1)
   data1_first=data1[1][2]                                                                                                                                                                                               
   data1_last=data1[-1][2]                                                                                                                                                                                               
   data2_first=data2[1][2]                                                                                                                                                                                               
   data2_last=data2[-1][2]                                                                                                                                                                                               
   Z_angle.append(calc_angle(linepts1, linepts2, data1_first, data1_last,  data2_first,data2_last))                                                                                                                      
                                                                                                                                                                                                                         
Z_angle=np.asarray(Z_angle)                                                                                                                                                                                              
print("\nAngles are saved %s:\n"%output)                                                                                                                                                                                 
print(Z_angle)                                                                                                                                                                                                           

# writing output in the output xvg file.
with open ("%s"%output,"w") as output:                                                                                                                                                                                   
    output.write("#\n#\n#\n@\n@\n@\n")                                                                                                                                                                                   
    for j in range(0, len(Z_angle)):                                                                                                                                                                                     
       output.write('\t{}\t{}'.format(j*100, Z_angle[j]) + "\n" )                                                                                                                                                        
output.close()                                                                                                                                                                                                           
                                                                                                                                                                                                                         
# Verify that everything looks right using 3D plot. Takes the last frame for plotting.                                                                                                                                                                                 
'''                                                                                                                                                                                                                      
import matplotlib.pyplot as plt                                                                                                                                                                                          
import mpl_toolkits.mplot3d as m3d                                                                                                                                                                                       
                                                                                                                                                                                                                         
ax = m3d.Axes3D(plt.figure())                                                                                                                                                                                            
ax.scatter3D(*data1.T)                                                                                                                                                                                                   
ax.plot3D(*linepts1.T, color='r', lw=4, markersize=10)                                                                                                                                                                   
ax.scatter3D(*data2.T)                                                                                                                                                                                                   
ax.plot3D(*linepts2.T, color='b', lw=4, markersize=10)                                                                                                                                                                   
plt.show()                                                                                                                                                                                                               
'''                                                                                                                                                                                                                      
                                                                                                                                                                                                                         91,1          Bot
                                                                                                                                                                                                                                          1,1           Top
