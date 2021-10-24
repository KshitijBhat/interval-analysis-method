#!/usr/bin/env python3
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
        self.theta1 = theta1
        self.theta2 = theta2
        self.endpoints = [Point(self.x+r*np.cos(theta1),self.y+r*np.sin(theta1)),Point(self.x+r*np.cos(theta2),self.y+r*np.sin(theta2))]
        self.circle = Circle(Point(x,y),r)  
        
    def encloses(self,point):
        xp,yp = float(point.x)-self.x,float(point.y)-self.y
        theta = np.arctan2(yp,xp) if yp>=0 else 2*np.pi + np.arctan2(yp,xp)
        if (theta>=self.theta1 and theta<=self.theta2) or \
           (theta>=self.theta2 and theta<=self.theta1):
            return True
        else:
            return False

    def draw(self,ax,color):
        ax.add_patch(Arc_patch((self.x, self.y), 2*self.r,2*self.r, theta1=np.rad2deg(self.theta1), theta2=np.rad2deg(self.theta2), linewidth=1, color=color))       

def arc_circle(gamma2,etaprime):
    if len(intersection(etaprime,gamma2.circle))>0:
        [endpt1,endpt2] = gamma2.endpoints
        incl = [etaprime.encloses(endpt1),etaprime.encloses(endpt2)]
        thetas = (gamma2.theta1,gamma2.theta2) if gamma2.theta2>gamma2.theta1 else (gamma2.theta2,gamma2.theta1)
        x,y = float(etaprime.center.x),float(etaprime.center.y)
        center_angle = np.arctan2(y,x) if y>=0 else 2*np.pi + np.arctan2(y,x)

        if any(incl) and not all(incl):
            return 1
        elif (center_angle<=thetas[1] and center_angle>=thetas[0]):
            if not any(incl) and center_angle<=np.pi:
                return 2
            if all(incl) and center_angle>=np.pi:
                return 2
    else:
        return 0

def arc_arc(gamma2,gamma2prime):
    """
    Returns True if gamma1 includes gamma2 partially or wholly
    Else returns False
    """
    eta,etaprime = gamma2.circle,gamma2prime.circle
    # if len(intersection(eta,etaprime))==0:
    #     return False
    # Lemma 3.3 and 3.4
    if (arc_circle(gamma2,etaprime)==2 and arc_circle(gamma2prime,eta)==2) or \
       (arc_circle(gamma2,etaprime)==2 and arc_circle(gamma2prime,eta)==1) or \
       (arc_circle(gamma2,etaprime)==1 and arc_circle(gamma2prime,eta)==2):
        return True
    # if (arc_circle(gamma2,etaprime)==1 and arc_circle(gamma2prime,eta)==1) and \
    #    (gamma2.theta1+gamma2.theta2>=0 and gamma2.theta1+gamma2.theta2<=2*np.pi) and \
    #    (gamma2prime.theta1+gamma2prime.theta2>=0 and gamma2prime.theta1+gamma2prime.theta2<=2*np.pi):
    #     return True
    

    else:
        intersect = intersection(eta,etaprime)
        if len(intersect) > 1:
            if gamma2.encloses(intersect[0]) and gamma2prime.encloses(intersect[0]):
                return True
        return False     

def line_arc(l1prime,gamma2):
    """
    Returns True if gamma2 includes l1prime partially or wholly
    Else returns False
    """
    intersect = intersection(l1prime,gamma2.circle)
    if len(intersect)>0:
        if gamma2.encloses(intersect[0]):
            return True
        try:
            if gamma2.encloses(intersect[1]):
                return True
            else:
                return False
        except IndexError:
            return False
    return False

def line_line(l1,l2):
    """
    Returns True if l1 includes l2 partially or wholly
    Else returns False
    """
    intersect = intersection(l1,l2)
    if len(intersect)>0:
        return True
    else:
        return False                            

class Interval:
    def __init__(self,r1,r2,theta1,theta2,origin=(0,0)):
        
        if r1>r2:
            r1,r2 = r2,r1
        
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
        r = np.linalg.norm(np.array((float(point[0]),float(point[1])))-np.array(self.origin))
        return (r<=self.r2 and r>=self.r1) and self.gamma1.encloses(point)    
        
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
            
            
        
    def bisect(self,i):
        if i == 0:
            L = Interval(self.r1, (self.r1+self.r2)/2, self.theta1, self.theta2)
            R = Interval((self.r1+self.r2)/2 ,self.r2, self.theta1, self.theta2)
        else:
            L = Interval(self.r1, self.r2, self.theta1, (self.theta1+self.theta2)/2)
            R = Interval(self.r1, self.r2, (self.theta1+self.theta2)/2, self.theta2) 
        return L,R
    
    def inclusion_test(self,Iprime):
        for i in [1,0]:
            for j in [1,0]:
                if arc_arc(self.gamma[i],Iprime.gamma[j]):
                    return True
                if line_arc(Iprime.l[i],self.gamma[j]):
                    return True
                if line_arc(self.l[i],Iprime.gamma[j]):
                    return True
                if line_line(self.l[i],Iprime.l[j]):
                    return True
        return False
    
    def interval_analysis(self,Iprime,Nr,Ntheta):
        if not self.inclusion_test(Iprime):
            J = Interval(0,0,0,0)
        else:
            N = Nr
            n = 0
            for i in [0,1]:
                L,R = self.bisect(i)
                while n<N:
                    LL,LR = L.bisect(i)
                    RL,RR = R.bisect(i)
                    if not R.inclusion_test(Iprime):
                        L = LL
                        R = LR
                    elif not L.inclusion_test(Iprime):
                        L = RL
                        R = RR
                    else:
                        if LL.inclusion_test(Iprime):
                            L = LL
                        else:
                            L = LR
                        if RR.inclusion_test(Iprime):
                            R = RR
                        else:
                            R = RL
                    n += 1
                if i == 0:    
                    Jr = [L.r1,R.r2]
                else:
                    Jtheta = [L.theta1,R.theta2]
                N = Ntheta
                n = 0
            J = Interval(*Jr,*Jtheta)    
        return J        