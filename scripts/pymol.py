from pymol import cmd
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
		0.,      0.,     1.]

@cmd.extend
def viewx(angle):
	cmd.set_view(Rx(angle) + setup)

@cmd.extend
def viewy(angle):
	cmd.set_view(Ry(angle) + setup)

@cmd.extend
def viewz(angle):
	cmd.set_view(Rz(angle) + setup)

counter = 0
@cmd.extend
def snapshot():
	cmd.set("mesh_radius", 0.02) 	# thinner mesh lines
	cmd.set("ray_shadows", "off")	# disable shadows
	cmd.set("antialias", 2)		# best quality
	cmd.set("hash_max", 4196)	# improve performance
	cmd.bg_color("white")		# make it a bit easier to see

	# guess file name
	objects = cmd.get_object_list("(all)")
	filename = ""
	for o in objects:
		if "SAS" in o:
			filename = o
			break
	if filename == "":
		print("Couldn't determine filename. Saving as 'temp'.")
		filename = "temp"

	global counter
	counter += 1
	cmd.png(filename+"_"+str(counter)+".png", width=1000, height=1000, dpi=600, ray=1)

# pymol doesn't seem to be able to parse arguments correctly, so we need a separate method for when we cannot use ray
@cmd.extend
def drawshot():
	global counter
	counter += 1
	cmd.set("opaque_background", "on")
#	cmd.draw()
	cmd.png("temp_"+str(counter)+".png", width=1000, height=1000, dpi=600, ray=0)
