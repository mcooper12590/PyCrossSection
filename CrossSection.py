from numpy import sin, cos, pi, fabs, sign, roll, arctan2, diff, cumsum, hypot, logical_and, where, linspace
from scipy import interpolate

class CrossSection:

	d = 1000

	def __init__(self, x, y):
		
		self.x = x
		self.y = y
		self.roll()

	def setUMPoint(self, umx, umy):

		self.umx = umx
		self.umy = umy
	
	def roll(self):

		self.xm = roll(self.x, 1)
		self.ym = roll(self.y, 1)
		self.xp = roll(self.x, self.x.size-1)
		self.yp = roll(self.y, self.y.size-1)

	def calcShapeParams(self):

		self.genL()
		self.calcA()

	def genL(self):

		self.l = hypot(self.x - self.xp, self.y - self.yp)
		self.pp = cumsum(self.l)

	def calcA(self):

		self.sA = (self.xm*self.y - self.x*self.ym).sum() * 0.5
		self.A = fabs(self.sA)

	def genRL(self):

		self.r_l = hypot(self.x-self.umx, self.y-self.umy)


	def findLR(self, h):

		ymin = self.y.min()
		a_h = ymin + h

		condL = logical_and(self.y > a_h, a_h > self.yp)
		condR = logical_and(self.y < a_h, a_h < self.yp)

		L = where(condL)[0][0] + 1
		R = where(condR)[0][0]
		return L,R

	def findCentroid(self):

		m = self.xm*self.y-self.x*self.ym
		cx = (1/(6*self.sA))*((self.x + self.xm)*m).sum()
		cy = (1/(6*self.sA))*((self.y + self.ym)*m).sum()

		return cx, cy
				
	def redraw(self, rl):

		alpha = arctan2(self.xp-self.xm, self.yp-self.ym)

		nx = self.x + sign(self.x)*rl*cos(alpha)
		ny = self.y - sign(self.x)*rl*sin(alpha)

		c = ccw(self.x, self.y, self.xm, self.ym, nx, ny)

		nx[c] = (self.x - sign(self.x)*rl*cos(alpha))[c]
		ny[c] = (self.y + sign(self.x)*rl*sin(alpha))[c]

		tck, u = interpolate.splprep([nx, ny], u=None, s=0.0)
		un = linspace(u.min(),u.max(),CrossSection.d)
		nx, ny = interpolate.splev(un, tck, der=0)

		self.x = nx
		self.y = ny
		self.roll()

def ccw(x, y, xm, ym, nx, ny):

	        return (x - xm) * (ny - ym) > (y - ym) * (nx - xm)
