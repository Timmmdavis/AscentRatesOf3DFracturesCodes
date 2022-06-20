from visualization import *
from properties import PlotProperties
import numpy

# loading simulation results
time_srs = np.asarray([10, 40, 3e4, 5e4, 5.43e4])
#endtime=10*60 #3600
#time_srs = np.asarray(np.linspace(1, endtime, num=5))
time_srs=30*60#6.95e4#1e6#6.95e4
endtime=time_srs

fracturestring='Data/PyFrac_volscl_2.461903399796411_vol_1.950_visc_0.05_c_crit_38.6392_G_8e+09_nu_0.25_K_Ic_2e+06_deltagamma_9810'
Fr_list, properties = load_fractures(address=f'./{fracturestring}',
                                     time_srs=time_srs)
time_srs = get_fracture_variable(Fr_list,
                                variable='time')


'''
'''
# plot footprint
Fig_FP = None
plot_properties = PlotProperties(plot_FP_time=False,line_width=10,line_color='tomato')
Fig_FP = plot_fracture_list(Fr_list,
                            variable='footprint',
                            projection='2D',
                            fig=Fig_FP,
                            plot_prop=plot_properties)



Fr_list, properties = load_fractures(address=f'./{fracturestring}',
                                     time_srs=endtime)
time_srs = get_fracture_variable(Fr_list,
                                variable='time')


Fig_FP = plot_fracture_list(Fr_list,
                            variable='w', 
                            fig=Fig_FP, 
                            projection='2D_clrmap')#lk w pf


widths=get_fracture_variable(Fr_list, 'w', edge=4, return_time=False)
meshes=get_fracture_variable(Fr_list, 'mesh', edge=4, return_time=False)
fluidpressure=get_fracture_variable(Fr_list, 'pf', edge=4, return_time=False)
netpressure=get_fracture_variable(Fr_list, 'pn', edge=4, return_time=False)
print(meshes[0])
print(meshes[0].CenterCoor)
print(meshes[0].CenterCoor[:,0])
print(meshes[0].CenterCoor[:,1])
x=[meshes[0].CenterCoor[:,0]]
y=[meshes[0].CenterCoor[:,1]]


all=np.dstack((x, y, widths,fluidpressure,netpressure)).squeeze()
#all=numpy.asarray([ np.transpose(x), np.transpose(y), np.transpose(widths) ])
numpy.savetxt("all.csv", all, delimiter=",")

#numpy.savetxt("x.csv", x, delimiter=",")
#numpy.savetxt("y.csv", y, delimiter=",")
#numpy.savetxt("Widths.csv", widths, delimiter=",")

plt.show()


