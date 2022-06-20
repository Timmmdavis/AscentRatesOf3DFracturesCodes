# -*- coding: utf-8 -*-
"""
This file is an edited version of the PyFrac example files.

Tim Davis GFZ Potsdam, 16-Nov-2019
"""

import numpy as np
import math
import time

# local imports
from mesh import CartesianMesh
from properties import MaterialProperties, FluidProperties, InjectionProperties, SimulationProperties
from fracture import Fracture
from controller import Controller
from fracture_initialization import Geometry, InitializationParameters

# solid properties - Water and Gas Shale
nu = 0.25                           # Poisson's ratio       [dmlss] 
G= 8e9                             # Shear modulus         [pa] 
youngs_mod = (2*G)*(1+nu)           # Young's modulus       [pa] 
Eprime = youngs_mod / (1 - nu ** 2) # plain strain modulus  [pa]
K_IcReported = 2*1e6                        # fracture toughness    [pa.sqrt(m)]
K_Ic = 2*1e6                        # fracture toughness    [pa.sqrt(m)]
fluiddensity=2000                   # [kg/m3]
rockdensity=3000                    # [kg/m3]
deltagamma=(rockdensity-fluiddensity)*9.81 # gradient in weight [pa.m^-1]
eta=0.05;                            # fluidviscosity [pa.s] - water=~1.1e-3


Cl=0.;

CriticalVolume=((1-nu)/(16*G))*(((9*(math.pi**4)*(K_IcReported**8))/(deltagamma**5))**(1/3)) #[m3]
VolumeIn=1.95;

###################################################################
#################### Some analytical insights #####################

#Max velocity 
Max_v=(4/(27*eta*G*math.pi**2))*VolumeIn*deltagamma**2*(1-nu);
#Volume in radius - Davis
VolFractureRadius=((9*G*VolumeIn)/(16*deltagamma*(1-nu)))**(1/4)   
timeFor2c=(VolFractureRadius*2)/Max_v;

#Some approximate fracture length scales derived analytically based on the parameters above
critradius=((3*math.pi*K_IcReported)/(8*deltagamma))**(2/3)                          #[m]
pressure=(2*deltagamma*critradius)/3                                         #[pa]
maxaperture=((8*(1-(nu**2)))/(math.pi*youngs_mod))*pressure*critradius       #[m]
#CriticalVolume=((1-nu)/(16*G))*(((9*(math.pi**4)*(K_IcReported**8))/(deltagamma**5))**(1/3)) #[m3]
volscl=VolumeIn/CriticalVolume
#################### Some analytical insights #####################
###################################################################



fluidviscosity=eta;
#Scaling this to transistion between trapped and propagating
volume=VolumeIn 
print(volume)

#Length and rate of injection (based on volume and parameters above)
rate=0.015#(volume/timeFor2c)*2#volume/lengthofinjection_s                  #[m/s] 
lengthofinjection_s=volume/rate;#6000/(K_Ic/1e6)      


#Some scales for the limits of the grid (mesh)
ylength=14
xlength=2
# creating mesh
steprad=(VolFractureRadius*0.75)
Mesh = CartesianMesh(steprad*xlength, steprad*ylength, int(xlength*35),  int(ylength*35)) 

#Predicted time to get to top of mesh:
z=steprad*ylength
#time at given z:
t_end=(((z**6*(eta**2*math.pi**4*G))/(9*VolumeIn**(3)*deltagamma**(3)*(1-nu )))**(1/6))**3;

print("z")
print(z)
print("src")
print((steprad*ylength)+VolFractureRadius*1.25)

#vertical offset
offset=2000;
def sigmaO_func(x, y):
    """ This function provides the confining stress over the domain"""

    # only dependant on the depth
    density = rockdensity
    return (-offset-y) * density * 9.8 


# material properties
Solid = MaterialProperties(Mesh,
                           Eprime,
                           toughness=K_Ic,
                           confining_stress_func=sigmaO_func,
                           minimum_width=maxaperture/100, 
                           free_surf=0,
                           Carters_coef=Cl) #minimum_width=1e-6  
                           #Carters_coef=1e-6)

# injection parameters
#The simulation consists of a 100 minutes injection (6000s) of water into a rock
Q0 = np.asarray([[0, lengthofinjection_s], [rate, 0]])
Injection = InjectionProperties(Q0,
                                Mesh,
                                source_coordinates=[0, -(steprad*ylength)+VolFractureRadius*1.25])




# fluid properties
Fluid = FluidProperties(viscosity=fluidviscosity,density=fluiddensity) #, 

# simulation properties
simulProp = SimulationProperties()
totaltime= t_end*10; #Arbitary large value
simulProp.finalTime =  totaltime                      # the time at which the simulation stops.
simulProp.gravity = True                              # take the effect of gravity into account.
#simulProp.front_advancing = 'explicit'                # possible options include 'implicit', 'explicit' and default... 'predictor-corrector'.
simulProp.enable_remeshing = False  #kill once at end of domain   


# Formatting strings and saving file
critradiusf="%3g" % critradius 
volumef="%5.3f" % volume       
deltagammaf="%3g" % deltagamma 
Gf="%3g" % G 
K_Icf="%3g" % K_Ic      
Filename='_'.join(['./Data/PyFrac',
  'volscl',str(volscl),
  'vol',str(volumef),
  'visc',str(fluidviscosity),
  'c_crit',str(critradiusf),
  'G',str(Gf),
  'nu',str(nu),
  'K_Ic',str(K_Icf),
  'deltagamma',str(deltagammaf)])
simulProp.set_outputFolder(Filename)      # the disk address where the files are saved.


simulProp.tolFractFront = 3.0e-4          # increase the tolerance for fracture front iteration
simulProp.toleranceEHL  = 0.003           # tolerance for the elastohydrodynamic system solver.
simulProp.toleranceProjection = 0.025     # tolerance for the toughness iteration.
simulProp.maxFrontItrs = 25               # maximum iterations for the fracture front.
simulProp.maxSolverItrs = 250             # maximum iterations for the elastohydrodynamic solver.
simulProp.maxProjItrs = 10                # maximum projection iterations.
simulProp.maxReattempts = 16              # maximum reattempts in case of time step failure.
simulProp.reAttemptFactor = 0.9           # the factor by which time step is reduced on reattempts.
simulProp.outputEveryTS = 0#5               # the time after the output is generated (saving or plotting).
simulProp.saveRegime=1  # save tip regime
simulProp.plotVar = ['w']   #, 'v']          # plot fracture width and fracture front velocity
simulProp.plotFigure=True
##Check confining stress:
#simulProp.plotVar = ['footprint']   #'w', 'v']          # plot fracture width and fracture front velocity
#simulProp.bckColor = 'confining stress'

''' - default tolerances

toleranceFractureFront = 1.0e-3         # tolerance for the fracture front position solver.
toleranceEHL = 1.0e-4                   # tolerance for the elastohydrodynamic system solver.
tol_projection = 2.5e-3                 # tolerance for the toughness iteration.

# max iterations
max_front_itrs = 25                     # maximum iterations for the fracture front.
max_solver_itrs = 80                    # maximum iterations for the elastohydrodynamic solver.
max_proj_Itrs = 10                      # maximum projection iterations.

# time step re-attempt
max_reattemps = 8                       # maximum reattempts in case of time step failure.
reattempt_factor = 0.8                  # the factor by which time step is reduced on reattempts.
'''

# initialization parameters
Fr_geometry = Geometry('radial', radius=VolFractureRadius/6)
init_param = InitializationParameters(Fr_geometry, regime='K')

# creating fracture object
Fr = Fracture(Mesh,
              init_param,
              Solid,
              Fluid,
              Injection,
              simulProp)

# create a Controller
controller = Controller(Fr,
                        Solid,
                        Fluid,
                        Injection,
                        simulProp)

# run the simulation
controller.run()


####################
# plotting results #
####################

from visualization import *

# loading simulation results
time_srs = np.linspace(1, totaltime, totaltime/10)
Fr_list, properties = load_fractures(address="./Data/Filename",
                                     time_srs=time_srs)
time_srs = get_fracture_variable(Fr_list,
                                 variable='time')

# plot footprint
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

#plotting in 3D
Fig_Fr = plot_fracture_list(Fr_list,
                            variable='mesh',
                            projection='3D')
Fig_Fr = plot_fracture_list(Fr_list,
                            variable='width',
                            projection='3D',
                            fig=Fig_Fr)
Fig_Fr = plot_fracture_list(Fr_list,
                            variable='footprint',
                            projection='3D',
                            fig=Fig_Fr)

plt.show(block=True)
#  set block=True and comment last 2 lines if you want to keep the window open
#plt.show(block=False)
#plt.pause(5)
#plt.close()

