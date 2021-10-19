#!/usr/bin/env python3
import numpy as np
from sympy.geometry import Point,Circle,Segment,intersection
import matplotlib.pyplot as plt
from matplotlib.patches import Arc as Arc_patch

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