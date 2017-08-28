# Optical_ray_tracer
Building and testing of optical ray tracer 


Dewi Gould
CID: 01051830
Optical Ray Tracer:  7/2/17

———————————————————————————————————

Modules included:
	•	OpRayTracer1.py
	•	Operator.py

———————————————————————————————————

GENERAL USAGE NOTES:
Module allows the production, propagation and visualisation of bundles of optical rays through a given system of optical lenses.
System operates in the x/y plane, with rays propagating in the z direction.

Operator.py is designed with a user-interface for ease of use. However this document will explain how to generate results from both the original code in OpRayTracer1 and the easy-to-use Operator file.

———————————————————————————————————

Inputs for Operator.py:

This has two user inputs:

Incident Beam Radius - suggested value of 6mm: user just types the number required. i.e. 6
Spot plot or RMS: If the user types ‘SPOT’ they will obtain a spot plot of the rays at the paraxial focal point. If the user types ‘RMS’ they will get the RMS deviation from the paraxial focal point.
This code will also plot the path of the rays and plot the initial positions of the rays, as well as giving the focal length of the optical configuration and the diffraction scale of the system.


The intention is for this code to make the collection of data quick and easy. However, the below documentation details how to obtain various pieces of information directly from the original code in OpRayTracer1.py

———————————————————————————————————
CREATING A BUNDLE OF RAYS AND PROPAGATING THROUGH TWO  REFRACTING OPTICAL ELEMENTS TO OUTPUT PLANE AT PARAXIAL FOCUS

import OpRayTracer1 as op
bun=op.Bundles(6,6,6)
lens1=op.SphericalRefraction(50,1,1.5168,0.02,100)
lens2=op.SphericalRefraction(55,1.5168,1,0,100)
focal_point=op.focal_point_finder([lens1,lens2]).focal_point()
outputplane=op.OutputPlane(focal_point)
f=op.final(bun, [lens1,lens2,outputplane])
f.spot_plot()
test=op.Test_Plotter(bun,[lens1,lens2,outputplane])
test.go()

This yields a spot plot and path plot for a bundle of rays propagating through a planoconvex lens.
Spot Plot is in the plane at z=paraxial focal point for the optical system, which is calculated using the class ‘focal_point_finder’
———————————————————————————————————
SPHERICAL REFLECTION EXAMPLE (convex reflecting surface)

import OpRayTracer1 as op
bun=op.Bundles(6,6,6)
a=op.SphericalReflection(10,-0.13,100)
b=op.OutputPlane(5)
test=op.Test_Plotter(bun,[a,b])
test.go()

Yields a plot of ray path after reflection off a convex reflecting surface - significant spherical aberration visible in this example.
———————————————————————————————————

TWO LENS OPTIMIZATION

import OpRayTracer1 as op
bun=op.Bundles(6,10,6)
m=op.mini(bun,0.02,0,5)
m.opt()

This minimises the RMS deviation for a given two-optical-surface system with respect to the two curvatures.
Generates the two curvatures which produce the minimum RMS deviation.

———————————————————————————————————
PLOT OF RMS AGAINST INCIDENT BEAM RADIUS

import OpRayTracer1 as op
r=op.rmsvarying(0.02,0)
r.rmsrelplot()

Yields a plot of RMS against incident beam radius for a two-optical-surface system to examine the relationship between the geometrical focus and the diffraction limit.

———————————————————————————————————
