import numpy as np
import pyiges

from lsdo_kit.design.design_geometry.design_geometry import DesignGeometry 
from lsdo_kit.old_files.mesh import Mesh

import matplotlib.pyplot as plt
import time

import vedo
from vedo import Points, Plotter, LegendBox

from lsdo_kit.cython.basis_matrix_surface_py import get_basis_surface_matrix
import scipy.sparse as sps

time_start = time.time()

''' Spars and ribs creation script '''
path_name = 'CAD/'
file_name = 'uCRM-9_wingbox.stp' #_wing_with_tip #_wing_structure_PENGoLINS
geo = DesignGeometry(path_name + file_name)

# path_name = 'CAD/'
# file_name = 'test_uCRM-9_wingbox.igs'
# geo.write_iges(path_name + file_name, plot = True)

# # Concatenates vertically all the linear matrices 
# geo.assemble()
# # Evaluate the physical coornadites of points to be fitted
# points_to_be_fitted = geo.evaluate()


