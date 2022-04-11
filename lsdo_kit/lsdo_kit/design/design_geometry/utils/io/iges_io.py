import numpy as np
from vedo import Points, Plotter, colors, LegendBox, show

from lsdo_kit.design.design_geometry.bsplines.bspline_surface import BSplineSurface

import pandas as pd
import re
import numpy as np
import copy
from vedo import Points, Plotter, colors, LegendBox, show


def read_iges(geo,file_name):
    """
    Read and import a geometry that's in an IGES format.
    Parameters
    ----------
    fileName : str
        File name of iges file. Should have .igs extension.
    """

    # print('The IGES reader still needs to be implemented. Please use the step reader.')
    # f = open(file_name, 'r')
    # Ifile = []
    # for line in f:
    #     line = line.replace(';', ',')  #This is a bit of a hack...
    #     Ifile.append(line)
    # f.close()

    # start_lines   = int((Ifile[-1][1:8]))
    # general_lines = int((Ifile[-1][9:16]))
    # directory_lines = int((Ifile[-1][17:24]))
    # parameter_lines = int((Ifile[-1][25:32]))

    # # Now we know how many lines we have to deal with
    # dir_offset  = start_lines + general_lines
    # para_offset = dir_offset + directory_lines

    # surf_list = []
    # # Directory lines is ALWAYS a multiple of 2
    # for i in range(directory_lines//2):
    #     # 128 is bspline surface type
    #     if int(Ifile[2*i + dir_offset][0:8]) == 128:
    #         start = int(Ifile[2*i + dir_offset][8:16])
    #         num_lines = int(Ifile[2*i + 1 + dir_offset][24:32])
    #         surf_list.append([start, num_lines])

    # nSurf = 1 #

    # print('Found %d surfaces in Iges File.'%(len(surf_list)))

    # surfs = []

    # for isurf in range(nSurf):  # Loop over our patches
    #     print(isurf)
    #     data = []
    #     # Create a list of all data
    #     # -1 is for conversion from 1 based (iges) to python
    #     para_offset = surf_list[isurf][0]+dir_offset+directory_lines-1

    #     for i in range(surf_list[isurf][1]):
    #         aux = Ifile[i+para_offset][0:69].split(',')
    #         for j in range(len(aux)-1):
    #             data.append(float(aux[j]))

    #     # Now we extract what we need
    #     Nctlu = int(data[1]+1)
    #     Nctlv = int(data[2]+1)
    #     ku    = int(data[3]+1)
    #     kv    = int(data[4]+1)

    #     counter = 10
    #     tu = data[counter:counter+Nctlu+ku]
    #     counter += (Nctlu + ku)

    #     tv = data[counter:counter+Nctlv+kv]
    #     counter += (Nctlv + kv)
    #     #print('tu',tu)
    #     #print('tv',tv)
    #     weights = data[counter:counter+Nctlu*Nctlv]
    #     weights = np.array(weights)
    #     if weights.all() != 1:
    #         print('WARNING: Not all weight in B-spline surface are 1. A NURBS surface CANNOT be replicated exactly')
    #     counter += Nctlu*Nctlv

    #     coef = np.zeros([Nctlu, Nctlv, 3])
    #     for j in range(Nctlv):
    #         for i in range(Nctlu):
    #             coef[i, j, :] = data[counter:counter +3]
    #             counter += 3

    #     # Last we need the ranges
    #     prange = np.zeros(4)

    #     prange[0] = data[counter    ]
    #     prange[1] = data[counter + 1]
    #     prange[2] = data[counter + 2]
    #     prange[3] = data[counter + 3]

    #     # Re-scale the knot vectors in case the upper bound is not 1
    #     tu = np.array(tu)
    #     tv = np.array(tv)
    #     '''
    #     if not tu[-1] == 1.0:
    #         tu /= tu[-1]

    #     if not tv[-1] == 1.0:
    #         tv /= tv[-1]
    #     '''
    # #return tu, tv, Nctlu, Nctlv, ku, kv #coef, 
    pass


def write_iges(geo, file_name, plot = False):
        """
        Write the surface to IGES format
        Parameters
        ----------
        fileName : str
            File name of iges file. Should have .igs extension.
        """
        if plot == True:
            vp_out = Plotter()
            vps = []
            for surf, color in zip(geo.output_bspline_entity_dict.values(), colors.colors.values()):
                vps.append(Points(surf.control_points, r=8, c = color).legend(surf.name))
            #TODO legend
            #lb = LegendBox(vps, nmax=i, width = 0.2, pad = 0, pos = "top-left")
            vp_out.show(vps, 'Control points', axes=1, viewup="z", interactive = False)
        f = open(file_name, 'w')
        print('Exporting', file_name)
        #TODO Change to correct information
        f.write('                                                                        S      1\n')
        f.write('1H,,1H;,7H128-000,11H128-000.IGS,9H{unknown},9H{unknown},16,6,15,13,15, G      1\n')
        f.write('7H128-000,1.,6,1HM,8,0.016,15H19970830.165254, 0.0001,0.,               G      2\n')
        f.write('21Hdennette@wiz-worx.com,23HLegacy PDD AP Committee,11,3,               G      3\n')
        f.write('13H920717.080000,23HMIL-PRF-28000B0,CLASS 1;                            G      4\n')
        Dcount = 1
        Pcount = 1
        for surf in geo.output_bspline_entity_dict.values():
            paraEntries = 13 + (len(surf.knots_u)) + (len(surf.knots_v)) + surf.shape[0] * surf.shape[1] + 3 * surf.shape[0] * surf.shape[1] + 1
            paraLines = (paraEntries - 10) // 3 + 2
            if np.mod(paraEntries - 10, 3) != 0:
                paraLines += 1
            f.write("     128%8d       0       0       1       0       0       000000001D%7d\n" % (Pcount, Dcount))
            f.write(
            "     128       0       2%8d       0                               0D%7d\n" % (paraLines, Dcount + 1)
            )
            Dcount += 2
            Pcount += paraLines
        Pcount  = 1
        counter = 1
        for surf in geo.output_bspline_entity_dict.values():
            f.write(
                "%10d,%10d,%10d,%10d,%10d,          %7dP%7d\n"
                % (128, surf.shape[0] - 1, surf.shape[1] - 1, surf.order_u - 1, surf.order_v - 1, Pcount, counter)
            )
            counter += 1
            f.write("%10d,%10d,%10d,%10d,%10d,          %7dP%7d\n" % (0, 0, 1, 0, 0, Pcount, counter))

            counter += 1
            pos_counter = 0

            for i in range(len(surf.knots_u)):
                pos_counter += 1
                f.write("%20.12g," % (np.real(surf.knots_u[i])))
                if np.mod(pos_counter, 3) == 0:
                    f.write("  %7dP%7d\n" % (Pcount, counter))
                    counter += 1
                    pos_counter = 0

            for i in range(len(surf.knots_v)):
                pos_counter += 1
                f.write("%20.12g," % (np.real(surf.knots_v[i])))
                if np.mod(pos_counter, 3) == 0:
                    f.write("  %7dP%7d\n" % (Pcount, counter))
                    counter += 1
                    pos_counter = 0

            for i in range(surf.shape[0] * surf.shape[1]):
                pos_counter += 1
                f.write("%20.12g," % (1.0))
                if np.mod(pos_counter, 3) == 0:
                    f.write("  %7dP%7d\n" % (Pcount, counter))
                    counter += 1
                    pos_counter = 0

            for j in range(surf.shape[1]):
                for i in range(surf.shape[0]):
                    for idim in range(3):
                        pos_counter += 1
                        cntrl_pts = np.reshape(surf.control_points, (surf.shape[0], surf.shape[1],3))
                        f.write("%20.12g," % (np.real(cntrl_pts[i, j, idim])))
                        if np.mod(pos_counter, 3) == 0:
                            f.write("  %7dP%7d\n" % (Pcount, counter))
                            counter += 1
                            pos_counter = 0

            for i in range(4):
                pos_counter += 1
                if i == 0:
                    f.write("%20.12g," % (np.real(surf.knots_u[0])))
                if i == 1:
                    f.write("%20.12g," % (np.real(surf.knots_u[1])))
                if i == 2:
                    f.write("%20.12g," % (np.real(surf.knots_v[0])))
                if i == 3:
                    f.write("%20.12g;" % (np.real(surf.knots_v[1])))
                if np.mod(pos_counter, 3) == 0:
                    f.write("  %7dP%7d\n" % (Pcount, counter))
                    counter += 1
                    pos_counter = 0
                else:  
                    if i == 3:
                        for j in range(3 - pos_counter):
                            f.write("%21s" % (" "))
                        pos_counter = 0
                        f.write("  %7dP%7d\n" % (Pcount, counter))
                        counter += 1

            Pcount += 2 
        f.write('S%7dG%7dD%7dP%7d%40sT%6s1\n'%(1, 4, Dcount-1, counter-1, ' ', ' '))
        f.close()  
        print('Complete export')