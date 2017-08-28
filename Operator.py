"""
This file is a user-interface to easily see results from the main code in OpRayTracer1.py
Input a desired bundle radius and choose whether to find RMS deviation or see a spot plot of rays at paraxial focus
The following example is for a planoconvex lens with curved surface (curvature 0.02/mm) facing the input

User will be able to see the diffraction limit, the RMS deviation and the focal point of the lens configuration.
"""
import OpRayTracer1 as op


choice1=float(raw_input("Bundle Radius - suggested value 6(mm): "))
choice2=str(raw_input("Type 'RMS' to see RMS deviation, or SPOT to see spot plot of rays: "))


bun=op.Bundles(6,choice1,6)
a=op.SphericalRefraction(50,1,1.5168,0.02,100)
b=op.SphericalRefraction(55,1.5168,1,0,100)
    
fp=op.focal_point_finder([a,b])
focalpoint=fp.focal_point()
print (focalpoint*(10**-3)), "Paraxial Focus (m)"
c=op.OutputPlane(focalpoint)

f=op.final(bun,[a,b,c])
if choice2 in ['SPOT', 'spot']:
    bun.initial_plot()
    f.spot_plot()
else:
    print (f.rms()*(10**-3)), "RMS (m)"

test =op.Test_Plotter(bun,[a,b,c])
test.go()
    
wavelength=(588)*(10**-9)
difscale=(wavelength*(focalpoint-50)*(10**-3))/(8*(10**-3))
print difscale, "Diffraction Scale, (m)"