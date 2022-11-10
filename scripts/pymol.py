from pymol import cmd, stored
from math import pi, cos, sin

view = cmd.get_view()
setup = list(view[9:])

# t is the angle in radians
def Rx(t):
	t = float(t)*pi
	return [1., 0.,         0., 
		0., cos(t), -sin(t),
		0., sin(t), cos(t)]

# t is the angle in radians
def Ry(t):
	t = float(t)*pi
	return [cos(t),  0., sin(t), 
		0., 1.,  0.,
		-sin(t), 0., cos(t)]

# t is the angle in radians
def Rz(t):
	t = float(t)*pi
	return [cos(t), -sin(t), 0., 
		sin(t), cos(t),  0.,
		0.,      0.,     0.]

@cmd.extend
def viewx(angle):
	cmd.set_view(Rx(angle) + setup)

@cmd.extend
def viewy(angle):
	cmd.set_view(Ry(angle) + setup)
	
@cmd.extend
def viewz(angle):
	cmd.set_view(Rz(angle) + setup)
