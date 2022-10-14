#exec(open("C:\\Users\\Berlin\\Documents\\GitHub\\PyFracTim\\examples\\GetFractureSpeed_PyFrac.py").read())

#To retrieve results:
#Sim 1.
#PyFrac_volscl_2.5375_vol_0.225_visc_0.005_c_crit_15.8573_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_49050
#Sim 2.
#PyFrac_volscl_2.5375_vol_1.950_visc_0.015_c_crit_49.0718_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_9810
#Sim 3.
#PyFrac_volscl_2.5375_vol_0.375_visc_0.005_c_crit_27.9923_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_19620
#Sim 4.
#PyFrac_volscl_2.5375_vol_0.150_visc_0.001_c_crit_20.4334_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_29430

###TO EXTRACT!
#10* appendix sim 2:
#PyFrac_volscl_2.5375_vol_7.500_visc_0.05_c_crit_47.4653_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_19620
#10* appendix sim 3
#PyFrac_volscl_2.5375_vol_3.750_visc_0.005_c_crit_47.4653_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_19620
#10* appendix sim 4
#PyFrac_volscl_2.5375_vol_3.750_visc_0.05_c_crit_47.4653_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_19620
#10* appendix sim 5

#Sim1*10


from visualization import *
import numpy
fracturestring='Data/PyFrac_volscl_2.461903399796411_vol_1.950_visc_0.05_c_crit_38.6392_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_9810_Explicit'
#fracturestring='Data/PyFrac_volscl_2.461903399796411_vol_1.950_visc_0.05_c_crit_38.6392_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_9810_Implicit'
Fr_list, properties = load_fractures(address=f'./{fracturestring}')

time_srs = get_fracture_variable(Fr_list,
                                 variable='time')


#print(properties.FluidProperties)

print(f"this is what is inside the tuple 'properties': %s" % (properties,))
print(len(properties))
# this lines UNPACKS values
# of variable a
(MatProp, FluidProp, InjProp, SimProp) = properties #extract from tuple
print(vars(InjProp))
SrcCoords=InjProp.sourceCoordinates
print(SrcCoords)

SrcX=SrcCoords[0]
SrcDepth=SrcCoords[1]
SrcXrnd=numpy.round_(SrcX, decimals=6)
print(SrcXrnd)
if SrcXrnd!=0:
  print(SrcXrnd)
  raise NameError('the injection location isnt at x=0')


#ReyNo=get_fracture_variable(Fr_list, 'Re', return_time=False)
# plot footprint
'''
Fig_FP = plot_fracture_list(Fr_list,
                            variable='mesh',
                            projection='2D')

Fig_FP = plot_fracture_list(Fr_list,
                            variable='footprint',
                            projection='2D',
                            fig=Fig_FP)

# plot slice
plot_prop = PlotProperties(line_style='.-')
Fig_WS = plot_fracture_list_slice(Fr_list,
                                  variable='w',
                                  projection='2D',
                                  plot_prop=plot_prop,
                                  plot_cell_center=True,
                                  orientation='vertical')
plt.show() #block=True
'''

#Starting to calc velocity - in postprocess_fracture.py
#get_fracture_variable(Fr_list, 'd_min', edge=3, return_time=False)


#Check volume is stable
get_fracture_variable(Fr_list, 'volume', edge=3, return_time=False)

data=fracturestring.split('_')
print(data)
print(data[2])
volume=float(data[4])
c_crit=float(data[9])
G=float(data[11])
nu=float(data[13])
K_Ic=2e6#float(data[14])
#deltagamma=float(data[-1])
print(volume)
print(c_crit)
print(G)
print(nu)
print(K_Ic)
#print(deltagamma)



youngs_mod = (2*G)*(1+nu)           # Young's modulus       [pa] 
c_crit2=(-((3*youngs_mod*volume)/(8*(math.pi)**0.5*K_Ic*(nu**2)-8*(math.pi)**0.5*K_Ic)))**(2/5)

#-------------------------------------------------------------------------------------
#Calculating velocity - all fractures
ylength=2
ylim=(round(c_crit2)*ylength)+1 #c_crit2!!!!!!!!!!
ylim=ylim*6

try:
  ReyNo=get_fracture_variable(Fr_list, 'Re', return_time=False) #front velocity
  print(ReyNo[1])
except:
  print("Not sure you saved this property (ReynoldsNo)")


#edge (int) – – the edge of the cell that will be plotted. 
#This is for variables that are evaluated on the cell edges instead of cell center. 
#It can have a value from 0 to 4 (0->left, 1->right, 2->bottom, 3->top, 4->average).

#Velocity of upper edge
Velocity=get_fracture_variable(Fr_list, 'v', edge=3, return_time=False) #front velocity
uppertipheight=np.zeros(len(Velocity))
lowertipheight=np.zeros(len(Velocity))
edgetip=np.zeros(len(Velocity))
TimeSeconds=np.zeros(len(Velocity))
Velocities=np.zeros(len(Velocity))

ylim=ylim*100
#print(ylim)
#error('herer')
ind=0;
ind2=0;
for j in range(0, len(Velocity)):
  ind=ind+1;
  ind2=ind2+1;
  #Get the widths and mesh of the fracture - in 'Pyfrac/src/postprocess_fracture.py'
  widths=get_fracture_variable(Fr_list, 'w', edge=4, return_time=False)
  meshes=get_fracture_variable(Fr_list, 'mesh', edge=4, return_time=False)
  frk1msh=meshes[-ind] #print(dir(frk1msh)) - to show contained fields
  frk1wdth=widths[-ind] 
  frk1velo=Velocity[-ind] 
  #Get xy coords of each element
  xy=frk1msh.CenterCoor
  #finding upper tip - if width is over 1 and largest z
  x=xy[:,0] #python indexing
  z=xy[:,1] #python indexing
  mnx=min(np.abs(x))
  #print(mnx)
  mnxrnd=numpy.round_(mnx, decimals=6)
  if mnxrnd!=0:
    #We assume:
    #The injection point is x=0
    #The central cell straddles x=0
    raise NameError('The centre cell should be at zero...')
  #Filter on x=0
  idxgoodx=np.nonzero(np.abs(x) == mnx)  
  zX=z[idxgoodx]
  frk1veloX=frk1velo[idxgoodx]
  frk1wdthX=frk1wdth[idxgoodx]
  #Filter on w>0
  idxgoodw=np.nonzero(frk1wdthX > 0)
  zXW=zX[idxgoodw]
  frk1veloXW=frk1veloX[idxgoodw]
  try:
    #Find highest z - this is your velocity!
    idxHighestCell=np.nanargmax(zXW)
    Velocities[ind2]=frk1veloXW[idxHighestCell]
    uppertipheight[ind2]=zXW[idxHighestCell]
    #Find lowest z - this is your velocity!
    lowertipheight[ind2]=np.nanmin(zXW)  
    #Get the times of each fracture [last 10]
    TimeSeconds[ind2]=time_srs[-ind];
    edgetip[ind2]=0
  except:
    continue

#Flip so corresponds to increasing time 
edgetip=np.flipud(edgetip)
uppertipheight=np.flipud(uppertipheight)
lowertipheight=np.flipud(lowertipheight)
TimeSeconds=np.flipud(TimeSeconds)
Velocities=np.flipud(Velocities)


#relative to src (0 is src)
heights=uppertipheight-lowertipheight
heightsFromSRC=uppertipheight-SrcDepth

print(lowertipheight)
print(SrcDepth)

print(heights)
print(TimeSeconds)
print(Velocities)

all=np.dstack((heights, TimeSeconds, Velocities, heightsFromSRC,edgetip)).squeeze()
#all=numpy.asarray([ np.transpose(x), np.transpose(y), np.transpose(widths) ])
numpy.savetxt("SpeedData.csv", all, delimiter=",")

#mean of values
#Velocitylst10Propstps=Velocitylst10Propstps[Velocitylst10Propstps!=0] #clip values of zero (and last)




#Calcuate velocity of each step
dz_t=np.zeros(len(Velocity)-1)
for j in range(0, len(dz_t)):
  #
  #print(j)
  #print("srt")
  if j==len(uppertipheight)-1:
    continue
  current_time=TimeSeconds[j]
  next_time=TimeSeconds[j+1]
  diff_time=next_time-current_time
  #
  current_tipheight=uppertipheight[j]
  next_tipheight=uppertipheight[j+1]
  diff_z=next_tipheight-current_tipheight
  #print(diff_z)
  #print(diff_time)
  dz_t[j]=diff_z/diff_time
  #print(dz_t[j])
  #catch some errors
  if uppertipheight[j]==uppertipheight[j+1]:
    dz_t[j]=np.NaN
  if diff_z==-uppertipheight[j] or diff_z==-uppertipheight[j+1] or dz_t[j]==0.0:
    dz_t[j]=np.NaN



print(dz_t)
'''
print(dz_t) 
#mean of my calculated speed values
print("mean of my calculated speed values")
print(np.nanmean(dz_t[-11:-1]))
#mean of programs tip velocity values
print("mean of programs tip velocity values")
print(np.nanmean(Velocities[-11:-1]))
'''

#Grab volumes
vols=get_fracture_variable(Fr_list, 'volume', edge=3, return_time=False)
Vnrm=np.zeros(len(vols))
for j in range(0, len(vols)):
  Vnrm[j]=vols[j]/volume;

#source_coords
src=-c_crit2*(ylength*0.7)
#relative to src (0 is src)
uppertipheightnrm=uppertipheight-src
#normalised relative to ccrit
uppertipheightnrm=uppertipheightnrm/c_crit2
edgetipnrm=edgetip/c_crit2


plt.scatter(dz_t,uppertipheightnrm[0:-1],c='blue')
plt.scatter(Velocities,uppertipheightnrm,c='red')
num=len(uppertipheight)
plt.scatter(Vnrm[0:num],uppertipheightnrm[0:num],c='black')
plt.scatter(edgetipnrm[0:num],uppertipheightnrm[0:num],c='pink')

#Between 1.5* and 2* crit length
idx = (uppertipheightnrm>(2.5))#*(uppertipheightnrm<(2.5))
#mean of my calculated speed values
print("mean of my calculated speed values in interval 1.5:2")
print(np.nanmean(dz_t[idx[0:-1]]))
#mean of programs tip velocity values
print("mean of programs tip velocity values in interval 1.5:2")
print(np.nanmean(Velocities[idx]))


all=np.dstack((x, y, widths)).squeeze()
#all=numpy.asarray([ np.transpose(x), np.transpose(y), np.transpose(widths) ])
numpy.savetxt("all.csv", all, delimiter=",")

#plt.axvline(x=nunnvelocity)
#plt.axvline(x=nunnvelocity/np.pi,ls='--')
#plt.axvline(x=1/np.pi,ls='--')
plt.xlabel('Velocity [m/s]')
plt.ylabel('UpperTipHeight/C_crit2')
print("Volume In:")
print(volume)
plt.xlim(-1, 2)
plt.show()

#plt.xlim(-2, 15)
#plt.ylim(0, 4)

