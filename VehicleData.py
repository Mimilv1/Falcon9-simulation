
class Rocket():
    def __init__(self, mass, heigth, radius):
        self.mass = mass
        self.heigth = heigth
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0 # heigth
        self.phi = 0.0 # rotation around x axis
        self.theta = 0.0 # rotation around y axis
        self.psi = 0.0 #rotation around z axis
        self.radius = radius

        

Faclon9 = Rocket(549000, 70,1.85)