from numpy import sin, cos, pi, sqrt, fabs, arctan2
import numpy as np
import matplotlib.pyplot as plt

d = 1000

def genEll(r1,r2,theta=0):
# generates x, y point ellipse based from two radii and rotation angle theta
	t=np.linspace(0, 2*pi-2*pi/d, d-1)
	x=r1*np.cos(t)
	y=r2*np.sin(t)


	if(theta!=0):

		tx=x
		ty=y

		x=tx*cos(theta)-ty*sin(theta)
		y=tx*sin(theta)+ty*cos(theta)
	
	return x, y

def genCirc(r):

	return genEll(r, r)

def genSemiCirc(r):

	t=np.linspace(-1*pi,0, d)
	x=r*cos(t)
	y=r*sin(t)
	return x, y

def genSemiEll(r1,r2):

	t=np.linspace(-1*pi,0, d)
	x=r1*np.cos(t)
	y=r2*np.sin(t)

	return x, y

