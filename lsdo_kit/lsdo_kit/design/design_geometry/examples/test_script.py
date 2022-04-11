import numpy as np
import pyiges

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 

from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
import scipy.sparse as sps

import matplotlib.pyplot as plt

from vedo import Points, Plotter, LegendBox

# for debugging in vscode only
import os
os.chdir("../lsdo_geo/lsdo_kit/design_geometry/examples")

''' Rib creation script - demonstrates how to import a step file and process of creating a component rib'''

# importing geometry and creating geometry object
path_name = 'CAD/'
file_name = 'eVTOL.stp'
geo = DesignGeometry(path_name + file_name)

print('hi')