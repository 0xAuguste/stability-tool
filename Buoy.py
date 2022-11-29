import pymesh
import numpy as np
import time

class Buoy():
    """Buoy class which contains a 3D mesh and physical parameters"""
    
    def __init__(self, buoy_file:str, mass:float, Cg:"list[float]"):
        self.mesh = pymesh.load_mesh(buoy_file)
        self.mass = mass
        self.CG = Cg # center of gravity
        self.CB = {} # center of buoyancy TODO: write these values to a csv
        self.water_line = {}