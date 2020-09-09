#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
from schimpy.schism_mesh import read_mesh


def create_arg_parser():
    """ Create ArgumentParser
    """
    parser = argparse.ArgumentParser(description='Identify small areas')
    parser.add_argument('--input_mesh',type=str,default=None,help='Input mesh')
    parser.add_argument('--warn',type=float,default=10.0,help = 'Threshold for warning (areas smaller)')
    parser.add_argument('--fail',type=float,default=4.,help = 'Threshold for failure (areas smaller)')
    return parser

def small_areas(mesh,warn=10.,fail=1.,logger=None):
    if logger is not None:
        logger.info("Checking for elements with small areas. Thresholds: warn={}, fail={}".format(warn,fail))
    if isinstance(mesh,str):
        mesh = read_mesh(mesh)
    areas = mesh.areas()
    centroids = mesh.centroids()

    sorted_order = np.argsort(areas)
    areas_sorted = areas[sorted_order]
    centroids_sorted = centroids[sorted_order]
    small, = np.where(areas_sorted < warn)
    warnings = []    
    for s in small:
        x,y = centroids_sorted[s]
        warnings.append("Global element: {:>9d} Area: {:.3f} Centroid: {:.1f},{:.1f}".format(sorted_order[s]+1,areas_sorted[s],x,y))

    for item in warnings: 
        if logger is None:
            print(item)
        else:
            logger.warning(item)
        
    if areas_sorted[0] < fail:
        raise ValueError("Mesh contains areas smaller than the failure threshold. Consult the log or printout above for areas and warnings")
        
    return
        


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    print(args.input_mesh)
    mesh = read_mesh(args.input_mesh)
    warn = args.warn
    fail = args.fail
    small_areas(mesh,warn,fail)
    
    
if __name__=="__main__":
    main()    
    
