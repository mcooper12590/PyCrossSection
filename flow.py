from numpy import log, sin, cos, pi, sqrt, roll, fabs, hypot
from CrossSection import *

def findH(g_h, cs, Q, S, rh):

	L, R = cs.findLR(g_h)
	wx = cs.x[L:R]
	wy = cs.y[L:R]
	wcs = CrossSection(wx, wy)
	wcs.calcShapeParams()
	Pw = wcs.pp[-1]
	A = wcs.A
	H = A/Pw

	C = 2.5*sqrt(9.81)*log(0.37*H/rh)

	u_bar = C*sqrt(H*S)
	Qcalc = u_bar*A

	return fabs(Q - Qcalc)

class flow:

#	def __init__(self, cs, Q, S, rh):
	def __init__(self, cs, Q, rh):

		self.cs = cs

		self.Q = Q
		self.rh = rh
#		self.S = S

		self.calcFlowParams()

	def calcFlowParams(self):

		self.u_bar = self.Q/self.cs.A
		self.P = self.cs.pp[-1]
		self.R = self.cs.A/self.P

#		self.calcChezy()
		self.calcUmax()
		self.genVGrad()
		self.calcphi()

		self.genT_b()

	def calcChezy(self):

		C = 2.5*sqrt(9.81)*log(0.37*self.R/self.rh)
		self.S = self.u_bar**2/(C**2*self.R)

	def calcUmax(self):

		ri = self.cs.r_l
		rh = self.rh
		cx = self.cs.umx
		cy = self.cs.umy

		f = (ri/(ri-rh))*(1+(rh/(ri*log(ri/rh)))-(1/log(ri/rh)))
		fr = roll(f, 1)

		fa = (f + fr)/2

		Ai = fabs( ( cx*(self.cs.y-self.cs.ym)+self.cs.x*(self.cs.ym-cy)+self.cs.xm*(cy-self.cs.y) ) / 2 )

		fw = fa*Ai

		self.umax = self.u_bar*(Ai.sum()/fw.sum())

	def genVGrad(self):

		self.vgrad = (self.umax/self.rh)*(1/log(self.cs.r_l/self.rh))

	def calcphi(self):

		sum = ((self.vgrad**2)*self.cs.l).sum()
		self.phi = (9.81*self.S)/sum

	def genT_b(self):

		self.T_b = self.phi*998.2*self.cs.A*self.vgrad**2
