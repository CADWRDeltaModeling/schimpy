from nodepy.runge_kutta_method import *
ssp104=loadRKM('SSP104')
ssp104.plot_stability_region(bounds=[-15,1,-10,10])