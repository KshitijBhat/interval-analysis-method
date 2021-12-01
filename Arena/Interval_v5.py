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
        
    def __contains__(self,point):
        xp,yp = float(point.x)-self.x,float(point.y)-self.y
        rot = np.array([[np.cos(self.theta1),-1*np.sin(self.theta1)],
                        [np.sin(self.theta1),np.cos(self.theta1)]])
        xp,yp = np.array([xp,yp]) @ rot
        theta = np.arctan2(yp,xp) if yp>=0 else 2*np.pi + np.arctan2(yp,xp)
        return (theta>=0 and theta<=self.theta2-self.theta1) or (theta>=self.theta2-self.theta1 and theta<=0)

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
        #if theta1>theta2:
            #theta1,theta2 = theta2,theta1
            
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
        return (r<=self.r2 and r>=self.r1) and (point in self.gamma1) 
    
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
        test_pts = []
        O = Point(self.origin)
        D = np.linalg.norm(np.array(self.origin)-np.array(Iprime.origin))
        if D>self.r2 + Iprime.r2:
            return Interval(0,0,0,0,self.origin)

        def cart2pol(pt):
            x,y = float(pt.x),float(pt.y)
            xi,yi = np.array([x,y])-np.array(self.origin)
            rho = np.linalg.norm(np.array([xi,yi]))
            phi = np.arctan2(yi,xi) if yi>=0 else 2*np.pi + np.arctan2(yi,xi)
            return [rho, phi]

        for gmp in Iprime.gamma:
            for gmi in self.gamma:
                pts = arc_arc(gmi,gmp)
                test_pts += pts    
            for li in self.l:
                pts = line_arc(li,gmp)
                test_pts += pts 
            if D<=self.r2 + Iprime.r2:
                for line in gmp.circle.normal_lines(O):
                    base_normal = line.points[0]
                    if base_normal in gmp:
                        if base_normal in self:
                            test_pts += [base_normal]
            if D**2<=self.r2**2 + Iprime.r2**2:                
                for line in gmp.circle.tangent_lines(O):
                    tangent_pt = line.points[1]
                    if tangent_pt in gmp:
                        if tangent_pt in self:
                            test_pts += [tangent_pt]    
        for lp in Iprime.l:
            for gmi in self.gamma:
                pts = line_arc(lp,gmi)
                test_pts += pts    
            for li in self.l:
                test_pts += line_line(li,lp) 
            base_perpendicular = lp.perpendicular_segment(O).points[1]
            if base_perpendicular in self:
                test_pts += [base_perpendicular]    


        for point in Iprime.gamma2.endpoints+Iprime.gamma1.endpoints:
            if point in self:
                test_pts+= [point]
        if not test_pts:
            return Interval(0,0,0,0,self.origin)
        rad_pts = np.array(list(map(cart2pol,test_pts)))
        R = rad_pts[:,0]
        THETA = rad_pts[:,1]
        rmax,rmin = max(R),min(R)
        thetamax,thetamin = max(THETA),min(THETA)
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

#intersections
def arc_arc(gamma2,gamma2prime):
    eta,etaprime = gamma2.circle,gamma2prime.circle
    intersect = intersection(eta,etaprime)
    pts = []
    if len(intersect) > 0:
        if intersect[0] in gamma2 and intersect[0] in gamma2prime:
            pts.append(intersect[0])
    if len(intersect) > 1:
        if intersect[1] in gamma2 and intersect[1] in gamma2prime:
            pts.append(intersect[1])    
        
    return pts

def line_line(l1,l2):
    return intersection(l1,l2)

def line_arc(l1prime,gamma2):
    intersect = intersection(l1prime,gamma2.circle)
    pts = []
    if len(intersect)>0:
        if intersect[0] in gamma2:
            pts.append(intersect[0])
    if len(intersect)>1:
        if intersect[1] in gamma2:
            pts.append(intersect[1]) 
    return pts  
