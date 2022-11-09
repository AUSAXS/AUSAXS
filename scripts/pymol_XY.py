from pymol import cmd, stored

Rx =   [1., 0., 0., 
		0., 0., 0., 
		0., 0., 0.]

Ry =   [0, 0, 0, 
		0, 1, 0, 
		0, 0, 0]

Rz =   [0, 0, 0, 
		0, 0, 0, 
		0, 0, 1]

camera_orig = [0., 0., 0.]
model_orig = [0., 0., 0.]
dfront = [400.]
drear = [700.]
orthoscopic = [0.]

cmd.set_view(Rx + camera_orig + model_orig + dfront + drear + orthoscopic)