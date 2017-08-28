import numpy as np
import pylab as pl
from scipy.optimize import fmin

class Ray:
    """
    Generates Optical Rays. 
    p0 = initial vector position of the ray
    k0 = initial vector direction of the ray
    A 'ray' is an instance of the class Ray, which will contain all the positions and directions of an optical ray.
    Functions p() and k() return the current position and direction respectively of a ray while vertices() returns all positions of a given 'ray'.
    """
    def __init__(self, p0=[0.0,0.0,0.0],k0=[0.0,0.0,0.0]):
        self.pp=np.array(p0,dtype=float)
        self.kk=np.array(k0,dtype=float)
        self.pall=[self.pp]
        self.kall=[self.kk]
    def p(self):
        return self.pall[-1]
    def k(self):
        return self.kall[-1]
    def append(self,otherp1,otherk1):
        self.pall.append(np.array(otherp1,dtype=float))
        self.kall.append(np.array(otherk1,dtype=float))
    def vertices(self):
        return self.pall
        
class Bundles:
    """
    'Bundles' Generator: generates a group of rays at initial positions,defined by n, rmax and m, in the x/y plane.
    Generates rigns of coordinates at radii r=irmax/n for i=0,1,2,....n. Each ring has n_i points uniformly distributed in polar angle.
    n_i = 1 for i =0, n_i=m*i for i>0.
    Generate Rays takes the initial positions created by rtuniform() and angles() and creates an array of 'rays' with a given initial direction: an array of instances of the class Ray
    """
    def __init__(self,nn,rrmax,mm):
        self.n=nn
        self.rmax=rrmax
        self.m=mm
    def angles(self,R,N):
        l=0
        for t in R:
            theta =0
            for i in range(N[l]):
                theta += ((2*np.pi)/N[l])
                yield t,theta
            l+=1
    def rtuniform(self):
        R=[]
        N=[]
        for i in range(self.n-1):
            r_i = (i*self.rmax)/self.n
            n_i=self.m*i
            R.append(r_i)
            N.append(n_i)
        N[0]=1
        return self.angles(R,N)
    def initial_plot(self):
        for r,t in self.rtuniform():
            pl.figure(1)
            pl.title("Spot Plot of Rays before propagation (z=0 plane)")
            pl.xlabel("x displacement (mm)")
            pl.ylabel("y displacement (mm)")
            pl.axis('scaled')
            pl.grid()
            pl.plot(r*np.cos(t),r*np.sin(t),'ro')
    def generate_rays(self):
        pos=[]
        for r,t in self.rtuniform():
            rr=np.array([0.0,0.0,0.0])
            xval=r*np.cos(t)
            yval=r*np.sin(t)
            rr[0]=xval
            rr[1]=yval
            rr[2]=0
            y=Ray(rr,[0,0,1])
            pos.append(y)
        return pos       
      
class OpticalElement:
    def propagate_ray(self,ray):
        "propagate a ray through the optical element"
        raise NotImplementedError()
        
class SphericalRefraction(OpticalElement):
    """
    Allows production of optical lenses of given position(z0), curvature(c), refractive index(n2) and aperture radius(ar). Positive curvature corresponds to centre of curvature at z>z0 
        """
    def __init__(self,z00,n10,n20,c0,ar0):
        self.z0=float(z00)
        self.n1=float(n10)
        self.n2=float(n20)
        self.c=float(c0)
        self.ar=float(ar0)
    def intercept(self,ray):
        "output is the surface normal at the point of intersection"
        pf1=ray.p()
        kf=ray.k()
        kfhat=kf/np.sqrt(np.dot(kf,kf))          
                 
        if self.c != 0:
            R=1.0/(self.c)
            z0=self.z0
            pf=pf1-(np.array((0,0,(z0+R)),dtype=float))
            magr=np.sqrt(np.dot(pf,pf))
            bracket= (np.dot(pf,kfhat)**2)-((magr**2)-(R**2))
            if bracket >=0:
                l_plus= float(-np.dot(pf,kfhat) + np.sqrt((np.dot(pf,kfhat)**2) - ((magr**2) - (R**2))))
                l_minus= float(-np.dot(pf,kfhat) - np.sqrt((np.dot(pf,kfhat)**2) - ((magr**2) - (R**2))))
            else:
                raise Exception("No real points of intersection")

            if l_plus>=0 and l_minus<0:
                vector = pf1 + l_plus*kfhat
                d=np.sqrt((vector[0]**2)+(vector[1]**2))
                if d<= self.ar:
                    l=vector-(np.array((0,0,(z0+R)),dtype=float))
                    return l
                else:
                    raise Exception( "Point of intersection outside aperture radius")
            elif l_minus>=0 and l_plus<0:
                vector = pf1 + l_minus*kfhat
                d=np.sqrt((vector[0]**2)+(vector[1]**2))
                if d<= self.ar:
                    l=vector-(np.array((0,0,(z0+R)),dtype=float))
                    return l
                else:
                    raise Exception( "Point of intersection outside aperture radius")
            elif  l_plus>=0 and l_minus>=0:
                if l_plus<l_minus:
                    vector = pf1 + l_plus*kfhat
                    d=np.sqrt((vector[0]**2)+(vector[1]**2))
                    if d<= self.ar:
                        l=vector-(np.array((0,0,(z0+R)),dtype=float))
                        return l
                    else:
                        raise Exception( "Point of intersection outside aperture radius")
                elif l_minus<l_plus:
                    vector = pf1+ l_minus*kfhat
                    d=np.sqrt((vector[0]**2)+(vector[1]**2))
                    if d<= self.ar:
                        l=vector-(np.array((0,0,(z0+R)),dtype=float))
                        return l
                    else:
                        raise Exception( "Point of intersection outside aperture radius")
                elif l_plus<0 and l_minus<0:
                    raise Exception( "No points of intersection in the forward direction")
        else:
            j= (self.z0-pf1[2])/(kfhat[2])
            if j>=0:
                vector = pf1 + j*kfhat
                d=np.sqrt((vector[0]**2)+(vector[1]**2))
                if d<= self.ar:
                    l=np.array((0,0,-1),dtype=float)
                    return l
                else:
                    raise Exception( "Point of intersection outside aperture radius")
            else:
                raise Exception( "No points of intersection in the forward direction")

    def refract(self,ray,norm):
        "Generates final direction of ray after refracting through lens"
        kf=ray.k()
        pp = (np.dot(kf,kf))
        kfhat=kf/np.sqrt(pp)
        n=norm
        nn=(np.dot(n,n))
        nhat=n/np.sqrt(nn)
        w = (np.dot(kfhat,nhat))
        theta1=np.arccos(w)
        stheta1=np.sin(theta1)
        a= float(self.n1/self.n2)
        m=1/a
        if stheta1<m:
            b=(np.cross(nhat,kfhat))
            c=-1*b
            d=np.cross(nhat,c)
            e=np.square(a)
            f=np.cross(nhat,kfhat)
            g=np.dot(f,f)
            h=1-(e*g)
            i=np.sqrt(h)
            r=i*nhat
            j=a*d
            kouth=j-r
            if self.c<0:
                kouth[2]=-kouth[2]
                return kouth
            else:
                return kouth
        else:
            raise Exception("Ray is totally internally reflected")
        
    def propagate_ray(self,ray):
        "Propagate ray to surface and refract it. If ray doesn't have point of intersection/ not in forward direction, error will be raised"
        pf11=ray.p()
        kf1=ray.k()
        kfhat1=kf1/np.sqrt(np.dot(kf1,kf1))
        if self.c != 0:
            R=1.0/(self.c)
            pnt= self.intercept(ray) + (np.array((0,0,(self.z0 + R)),dtype=float))
        else:
            j= (self.z0 - (pf11[2]))/(kfhat1[2])
            pnt = pf11 + (j*kfhat1)
        dnew=self.refract(ray,self.intercept(ray))
        ray.append(pnt,dnew)

class SphericalReflection(OpticalElement):
    """
    Creates a spherical reflecting surface.
    z0 is intersection of surface with optical axis
    c is curvature, and is a signed quantity
    ar is the aperture radius
    Spherical reflection can calculate the intercept of a ray with the surface and propagate the ray to that point and append the new reflected direction
    """
    def __init__(self,z0,c,ar):
        self.z0=z0
        self.c=c
        self.ar=ar
    def reflect(self,ray,norm):
        kf=ray.k()
        pp = (np.dot(kf,kf))
        kfhat=kf/np.sqrt(pp)
        n=norm
        nn=(np.dot(n,n))
        nhat=n/np.sqrt(nn)
        kout= kfhat - (2*nhat*(np.dot(kfhat,nhat)))
        return kout
    def propagate_ray(self,ray):
        pf11=ray.p()
        kf1=ray.k()
        kfhat1=kf1/np.sqrt(np.dot(kf1,kf1))
        a=SphericalRefraction(self.z0,1,1,self.c,self.ar)
        if self.c != 0:
            R=1.0/(self.c)
            pnt= a.intercept(ray) + (np.array((0,0,(self.z0 + R)),dtype=float))
        else:
            j= (self.z0 - (pf11[2]))/(kfhat1[2])
            pnt = pf11 + (j*kfhat1)
        dnew=self.reflect(ray,a.intercept(ray))
        ray.append(pnt,dnew)
        

class OutputPlane(OpticalElement):
    """
    Finds point of intersection with output plane and can propagate ray to the output plane"
    """
    def __init__(self,zz):
        self.z=zz
    def intercept(self,ray):
        df=ray.k()
        dd = np.sqrt((np.dot(df,df)))
        dfhat=df/dd
        pf=ray.p()
        a= (self.z - pf[2])/(dfhat[2])
        pt=pf + (a*df)
        return pt
    def propagate_ray(self,ray):
        ray.append((self.intercept(ray)),ray.k())  
 
class focal_point_finder:
    """
    Method to find focal point of given optical system
    """
    def __init__(self, Elementt=[]):
        self.element=Elementt
        self.ray=Ray([0.1,0,0],[0,0,1])
    def focal_point(self):
        for i in self.element:
            i.propagate_ray(self.ray)
        ppp=self.ray.p()
        ddd=self.ray.k()
        if ddd[0]==0.0:
            raise Exception("No Focal Point")
        lll=-(ppp[0])/(ddd[0])
        if lll<=0:
            raise Exception("No Focal Point")
        else:
            fp=ppp[-1] + (lll*ddd[-1])
            return fp  
class Test_Plotter:
    """
    Visualisation of paths of rays through an optical system
    """
    def __init__(self,bundle,ell=[]):
        self.bundles=bundle
        self.lenses=ell
        self.rayss=[]
    def go(self):
        self.rayss.append(self.bundles.generate_rays())
        l=[]
        for k in self.rayss:
            for i in k:
                for j in self.lenses:
                    j.propagate_ray(i)
                l.append(i.vertices())
 
        for j in l:
            x=[]
            y=[]
            z=[]
            for k in j:
                x.append(k[0])
                y.append(k[1])
                z.append(k[2])
                pl.figure(4)
                pl.title("Ray displayed in x/z plane")
                pl.xlabel("z displacement (mm)")
                pl.ylabel("x displacement (mm)")
                pl.plot(z,x,'r')
                pl.autoscale(enable=True, axis='y')
                pl.autoscale(enable=True, axis='x')
                pl.grid()
        pl.show()
       
class final:
    """
    Propagates a given group/ bundle of rays through an optical system. Generates spot plots and RMS deviation of rays from paraxial focus.
    """
    def __init__(self,Bundle, Elementss=[]):
        self.bundle=Bundle
        self.elements=Elementss
        self.rays=[]
    def create_rays(self):
        self.rays.append(self.bundle.generate_rays())
    def get_rays(self):
        return self.rays
    def propagate_rays(self):
        for i in self.elements:
            for j in self.rays:
                for k in j:
                    i.propagate_ray(k)
    def finalise(self):
        self.create_rays()
        self.propagate_rays()        
    def spot_plot(self):
        self.finalise()
        xx=[]
        yy=[]
        for j in self.rays:
            for i in j:
                p=i.p()
                yy.append(p[1])
                xx.append(p[0])
                pl.figure(3)
                pl.title("Spot Diagram in x/y plane (z=paraxial focus)")
                pl.xlabel("x displacement (mm)")
                pl.ylabel("y displacement (mm)")
                pl.axis('scaled')
                pl.grid()
                pl.plot(yy,xx,'ro')
        pl.show()
    def rms(self):
        self.finalise()
        rr=[]
        for i in self.rays:
            for j in i:
                rr.append(((j.p()[0])**2) + ((j.p()[1]**2)))
        h=0
        for i in rr:
           h+=i
        rms= np.sqrt((1./(len(rr)))*h)
        return rms

class mini:
    """ 
    Method to minimise RMS deviation from focal point of ray bundle through two lens system
    c1g and c2g are 'guesses' for values of the curvature - suggested values 0.02 and 0
    d is the separation in mm between the two optical elements - suggest value 5
    """
    def __init__(self,bun,c1g,c2g,d):
        self.bun=bun
        self.c1g=c1g
        self.c2g=c2g
        self.d=d
    def twolens(self,c):
        lens1=SphericalRefraction(50,1,1.5168,c[0],100)
        lens2=SphericalRefraction(50+self.d,1.5168,1,c[1],100)
        fp=focal_point_finder([lens1,lens2])
        try:
            focalpoint=fp.focal_point()
            lens3=OutputPlane(focalpoint)
            f=final(self.bun,[lens1,lens2,lens3])
            RMS=f.rms()*(10**-3)
            return RMS
        except:
            return 100000000

    def opt(self):
        res=fmin(self.twolens,(self.c1g,self.c2g))
        print res[0],res[1], "Curvature 1 and Curvature 2"

class rmsvarying:
    """
    Produces plot of relationship between RMS deviation around paraxial focal point and radius of incident beam
    c1 and c2 are curvatures of the two spherical elements
    """
    def __init__(self,c1,c2):
        self.c1=c1
        self.c2=c2
    def rmsrelplot(self):
        bundlers=[]
        Radius=[]
        DifScale=[]
        a=SphericalRefraction(50,1,1.5168,self.c1,100)
        b=SphericalRefraction(55,1.5168,1,self.c2,100)
        fp=focal_point_finder([a,b])
        focalpoint=fp.focal_point()
        c=OutputPlane(focalpoint)
        wavelength=(588)*(10**-9)
        for i in range(15):
            bundlers.append(Bundles(10,i,6))
            Radius.append(i)
            DifScale.append((wavelength*focalpoint*(10**-3))/(i*(10**-3)))

        rms=[]
        for j in bundlers:
            f=final(j,[a,b,c])
            rms.append((f.rms()*(10**-3)))
        pl.figure(5)
        pl.title("RMS deviation of ray from paraxial focus vs Radius of incident beam (mm), red line is Diffraction Scale")
        pl.scatter(Radius,rms)
        pl.plot(Radius,DifScale,'r')
        pl.xlabel("Radius of incident beam (mm)")
        pl.ylabel("RMS deviation of ray from paraxial focus")
        pl.ylim([-0.00005,0.0002])
        pl.show()      