import numpy as np
from sympy.geometry import Point,Circle,Segment,intersection
import matplotlib.pyplot as plt
from matplotlib.patches import Arc as Arc_patch
from matplotlib import collections  as mc

class Arc:
    def __init__(self,r,theta1,theta2,origin = (0,0)):
        x,y = origin
        self.x = x
        self.y = y
        self.r = r
        self.origin = origin
        self.theta1 = theta1
        self.theta2 = theta2
        self.endpoints = [Point(self.x+r*np.cos(theta1),self.y+r*np.sin(theta1)),Point(self.x+r*np.cos(theta2),self.y+r*np.sin(theta2))]
        self.circle = Circle(Point(x,y),r)  
        
    def encloses(self,point):
        xp,yp = float(point.x)-self.x,float(point.y)-self.y
        rot = np.array([[np.cos(self.theta1),-1*np.sin(self.theta1)],
                        [np.sin(self.theta1),np.cos(self.theta1)]])
        xp,yp = np.array([xp,yp]) @ rot
        theta = np.arctan2(yp,xp) if yp>=0 else 2*np.pi + np.arctan2(yp,xp)
        if (theta>=0 and theta<=self.theta2-self.theta1) or \
           (theta>=self.theta2-self.theta1 and theta<=0):
            return True
        else:
            return False

    def draw(self,ax,color):
        ax.add_patch(Arc_patch((self.x, self.y), 2*self.r,2*self.r, theta1=np.rad2deg(self.theta1), theta2=np.rad2deg(self.theta2), linewidth=1, color=color))       


class Interval:
    def __init__(self,r1,r2,theta1,theta2,origin=(0,0)):
        
        if r1>r2:
            r1,r2 = r2,r1
        if theta1 <0 and theta2<0:
            theta1 = 2*np.pi + theta1
            theta2 = 2*np.pi + theta2
#         if theta1>=2*np.pi:
#             theta1 = theta1 - 2*np.pi
#         if theta2>=2*np.pi:
#             theta2 = theta2 - 2*np.pi    
        #theta1,theta2 = self.correction([theta1,theta2])    
        if theta1>theta2:
            theta1,theta2 = theta2,theta1
            
        not_zero = lambda x: x if x != 0.0 else 1e-6    
        r1,r2,theta1,theta2 =  list(map(not_zero,[r1,r2,theta1,theta2]))  
        
        self.r1 = r1
        self.r2 = r2  
        self.theta1 = theta1
        self.theta2 = theta2
        x,y = origin
        self.origin = origin
        self.x = x
        self.y = y
        self.gamma1 = Arc(r1,theta1,theta2,origin=origin)
        self.gamma2 = Arc(r2,theta1,theta2,origin=origin)
        self.l1 = Segment((x+r1*np.cos(theta1),y+r1*np.sin(theta1)),(x+r2*np.cos(theta1),y+r2*np.sin(theta1)))
        self.l2 = Segment((x+r1*np.cos(theta2),y+r1*np.sin(theta2)),(x+r2*np.cos(theta2),y+r2*np.sin(theta2)))
        self.l = [self.l1,self.l2]
        self.gamma = [self.gamma1,self.gamma2]

    def __contains__(self,point):
        r = np.linalg.norm(np.array((float(point.x),float(point.y)))-np.array(self.origin))
        return (r<=self.r2 and r>=self.r1) and self.gamma1.encloses(point)    
    
    def __repr__(self):
        return f"Interval({self.r1},{self.r2},{self.theta1},{self.theta2},{self.origin})"
    
    def draw(self,ax,color):
        try:
            l11, l12 = tuple(self.l1.points[0]),tuple(self.l1.points[1])
            l21, l22 = tuple(self.l2.points[0]),tuple(self.l2.points[1])
            lc = mc.LineCollection([[l11,l12],[l21,l22]], colors = color, linewidths=1)
            self.gamma1.draw(ax,color)
            self.gamma2.draw(ax,color)
            ax.add_collection(lc) 
            
        except AttributeError:
            self.gamma1.draw(ax,color)
            self.gamma2.draw(ax,color)
    def fast_analysis(self,Iprime):
        def cart2pol(pt):
            x,y = float(pt.x),float(pt.y)
            xi,yi = np.array([x,y])-np.array(self.origin)
            rho = np.linalg.norm(np.array([xi,yi]))
            phi = np.arctan2(yi,xi) if yi>=0 else 2*np.pi + np.arctan2(yi,xi)
            return [rho, phi]
        def expts(**kwargs):
            if 'gamma' in kwargs and 'gammaprime' in kwargs:
                gamma = kwargs['gamma']
                gammaprime = kwargs['gammaprime']
                d = np.linalg.norm(np.array(gamma.origin)-np.array(gammaprime.origin))
                if d >= gamma.r:
                    rmin = d - gammaprime.r
                    return [[rmin],[]]
                else:
                    thetae = []
                    rmin = d - gammaprime.r
                    phi = cart2pol(Point(Iprime.origin))[1]
                    thetalol = gammaprime.r/d
                    if thetalol>1 or thetalol<-1:
                        thi = 0
                    else:
                        thi = np.arcsin(thetalol)
                    if (phi+thi>=self.theta1 and phi+thi<=self.theta2):
                        thetae.append((phi + thi )) 
                    if (phi-thi>=self.theta1 and phi-thi<=self.theta2):
                        thetae.append((phi - thi)) 
              
                    return [[rmin],thetae]
                    
                return [[],[]]
            if 'l' in kwargs and 'gammaprime' in kwargs:
                l = kwargs['l']
                gammaprime = kwargs['gammaprime']
                d = np.linalg.norm(np.array(self.origin)-np.array(gammaprime.origin))
                phi = cart2pol(Point(Iprime.origin))[1]
                thetalol = gammaprime.r/d
                if thetalol>1 or thetalol<-1:
                    thi = 0
                else:
                    thi = np.arcsin(thetalol)
                thetae = phi + thi*(1 if (phi+thi>=self.theta1 and phi+thi<=self.theta2) else -1)
                return [thetae]
            if 'gamma' in kwargs and 'lprime' in kwargs:
                lprime = kwargs['lprime']
                gamma = kwargs['gamma']
                seg = Segment(*intersection(lprime,gamma.circle))
                return [seg.midpoint]
            


        rs = []
        thetas = []
        test_pts = []
        arc_cut = []

        for gmi in self.gamma:
            for gmp in Iprime.gamma:
                pts = arc_arc(gmi,gmp)
                if len(pts) == 2:
                    expt = expts(gamma=gmi,gammaprime=gmp)
                    rs += expt[0]
                    thetas += expt[1]
                test_pts += pts  
            for lp in Iprime.l:
                pts = line_arc(lp,gmi)
                if len(pts) == 2:
                    test_pts += expts(lprime=lp,gamma=gmi)
                test_pts += pts 
                arc_cut.append(1)
        for li in self.l:
            for gmp in Iprime.gamma:
                pts = line_arc(li,gmp)
                if len(pts)==2:
                    thetas += expts(l=li,gammaprime=gmp) 
                test_pts += pts    
            for lp in Iprime.l:
                pts = line_line(li,lp) 
                test_pts += pts
                if len(pts)!=0 and len(arc_cut)!=0:
                    ptx,pty = float(pts[0].x),float(pts[0].y)
                    d = np.linalg.norm(np.array([ptx,pty])-np.array(self.origin))
                    inta = intersection(lp,Circle(self.origin,d))
                    if len(inta) == 2:
                        midpt = Segment(*inta).midpoint
                    if len(inta) == 1: 
                        midpt = inta 
                    else:
                        midpt = []      
                    test_pts += midpt
                
        if self.gamma2.encloses(Point(Iprime.origin)) and Iprime.gamma2.encloses(Point(self.origin)):
            expt = expts(gamma=self.gamma2,gammaprime=Iprime.gamma2)
            rs += expt[0]
            thetas += expt[1]
        for point in Iprime.gamma2.endpoints+Iprime.gamma1.endpoints:
            if point in self:
                test_pts+= [point]

        rad_pts = np.array(list(map(cart2pol,test_pts)))
        if len(rad_pts) == 0:
            return Interval(0,0,0,0,self.origin)
        
        radii = list(rad_pts[:,0]) + rs
        
        thetas = self.correction(list(rad_pts[:,1])) + thetas
        rmax,rmin = max(radii),min(radii)
        thetamax,thetamin = max(thetas),min(thetas)
        if thetamax-thetamin>np.pi:
            thetamax = thetamax-2*np.pi
        return Interval(rmin,rmax,thetamin,thetamax,self.origin) 
    
    @staticmethod
    def correction(thetas):
        tuts = []
        tuts1 = []
        for theta in thetas:
            if theta>=0 and theta<= np.pi/2:
                tuts.append(0)
            if theta>=3*np.pi/2 and theta<= 2*np.pi:
                tuts.append(1)
                tuts1.append(1)
            if theta<0 and theta>= -1*np.pi/2:
                tuts1.append(0)    

        if any(tuts) and not all(tuts):
            return [(theta-2*np.pi if theta>=3*np.pi/2 and theta<= 2*np.pi else theta) for theta in thetas]
        if any(tuts1) and not all(tuts1):
            return [(theta+2*np.pi if theta>-1*np.pi/2 and theta< 0 else theta) for theta in thetas]    
        else:
            return thetas
        
    def complement(self,J):
        if J.r2 == 1e-6 and J.theta1 - J.theta2 == 0:
            return [self]
        Is = []
        if abs(self.theta1 - J.theta1)>1e-4 and abs(self.theta1 - J.theta1)<2*np.pi-1e-4:
            I1 = Interval(self.r1,self.r2,*self.correction([self.theta1,J.theta1]),self.origin)
            Is.append(I1)
            
        if abs(J.theta1 - J.theta2)>1e-4 and abs(J.theta1 - J.theta2)<2*np.pi-1e-4:    
            I2 = Interval(self.r1,J.r1,*self.correction([J.theta1,J.theta2]),self.origin)
            Is.append(I2)

        if abs(J.theta2 - self.theta2)>1e-4 and abs(J.theta2 - self.theta2)<2*np.pi-1e-4:    
            I3 = Interval(self.r1,self.r2,*self.correction([J.theta2,self.theta2]),self.origin)
            Is.append(I3)

        return Is   

def arc_arc(gamma2,gamma2prime):
    eta,etaprime = gamma2.circle,gamma2prime.circle
    intersect = intersection(eta,etaprime)
    pts = []
    if len(intersect) > 0:
        if gamma2.encloses(intersect[0]) and gamma2prime.encloses(intersect[0]):
            pts.append(intersect[0])
    if len(intersect) > 1:
        if gamma2.encloses(intersect[1]) and gamma2prime.encloses(intersect[1]):
            pts.append(intersect[1])    
        
    return pts

def line_arc(l1prime,gamma2):
    intersect = intersection(l1prime,gamma2.circle)
    pts = []
    if len(intersect)>0:
        if gamma2.encloses(intersect[0]):
            pts.append(intersect[0])
    if len(intersect)>1:
        if gamma2.encloses(intersect[1]):
            pts.append(intersect[1]) 
    return pts

def line_line(l1,l2):
    return intersection(l1,l2)            