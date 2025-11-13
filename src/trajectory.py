import numpy as np
import matplotlib.pyplot as plt
import mpld3
from ss3d import SolarSystem, SolarSystemBody

#app = Flask(__name__)
solar_system = SolarSystem(400)
mpld3.show()
body = SolarSystemBody(solar_system, 100, velocity=(100, 100, 100))

for _ in range(100):
    solar_system.update_all()
    solar_system.draw_all()