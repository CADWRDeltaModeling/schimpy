# -*- coding: utf-8 -*-
""" 3D Version of schism_mesh
"""

from .triquadmesh import TriQuadMesh, BoundaryType, EdgeType
from .schism_vertical_mesh import read_vmesh
import osgeo.ogr
import osgeo.osr
import numpy as np
import os
import math
from copy import deepcopy
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, PatchCollection
import matplotlib.patches as patches
from shapely.geometry import Point, Polygon, MultiLineString
import rtree.index

__all__ = ['SchismMesh', 'BoundaryType', 'read_mesh', 'write_mesh',]


def find_intersection(s1, s2):
    """ Calculate an intersection point from two line segments
        defined by two end points
    """
    det = (s1[0, 0] - s1[1, 0]) * (s2[0, 1] - s2[1, 1]) - \
          (s1[0, 1] - s1[1, 1]) * (s2[0, 0] - s2[1, 0])
    if det == 0.:  # parallel
        return None
    x = (s1[0, 0] * s1[1, 1] - s1[0, 1] * s1[1, 0]) * \
        (s2[0, 0] - s2[1, 0]) - \
        (s1[0, 0] - s1[1, 0]) * \
        (s2[0, 0] * s2[1, 1] - s2[0, 1] * s2[1, 0])
    y = (s1[0, 0] * s1[1, 1] - s1[0, 1] * s1[1, 0]) * \
        (s2[0, 1] - s2[1, 1]) - \
        (s1[0, 1] - s1[1, 1]) * \
        (s2[0, 0] * s2[1, 1] - s2[0, 1] * s2[1, 0])
    intersection = np.array([x, y]) / det
    v1 = s1[0, :] - intersection
    v2 = s1[1, :] - intersection
    sign1 = np.dot(v1, v2)
    v1 = s2[0, :] - intersection
    v2 = s2[1, :] - intersection
    sign2 = np.dot(v1, v2)
    return None if sign1 > 0. or sign2 > 0. else intersection


class SchismBoundary(object):
    """ A class describing an individual boundary.
    """

    def __init__(self, nodes, btype, comment=None):
        self._nodes = nodes
        self._btype = btype
        self._comment = comment

    def __repr__(self):
        return self._nodes

    @property
    def btype(self):
        """ Boundary type.

            :getter: Get the type of the boundary.
        """
        return self._btype

    @property
    def nodes(self):
        """ A sequence of node indexes in the boundary.

            :getter: Get the sequence of nodes in the boundary.
        """
        return self._nodes

    def n_nodes(self):
        """ Get the total number of nodes in the boundary.

            Return
            ------
            integer
                the total number of nodes in the boundary.
        """
        return len(self._nodes)

    @property
    def comment(self):
        """ Get a comment or name of the boundary.
        """
        return self._comment

    @comment.setter
    def comment(self, value):
        """ Set a comment or name of the boundary
        """
        self._comment = value


class SchismMesh(TriQuadMesh):
    """ Memory model of 3D SCHISM mesh
    """

    def __init__(self):
        super(SchismMesh, self).__init__()
        self._boundaries = []
        self._metadata = {}
        self._vmesh = None
        self._areas = None
        self._points = None
        self._polygons = None
        self._edge_lengths = None
        self._centroids = None
        self._merged_mesh = None

    @property
    def vmesh(self):
        return self._vmesh

    @property
    def boundaries(self):
        """ An array of the boundary information

            :getter: Get the array of the boundary information
        """
        return self._boundaries

    def n_boundaries(self, btype=None):
        """ Get the number of boundaries. If a boundary type is given,
            it counts only boundaries with the corresponding type.

            Parameters
            ----------
            btype: integer, optional
                Type of the boundary.
        """
        if btype is None:
            return len(self._boundaries)
        else:
            return sum(b.btype == btype for b in self._boundaries)

    def n_total_boundary_nodes(self, btype):
        """ Get the total node boundary of a given type

            Parameters
            ----------
            btype: integer
                Type of the boundary.
        """
        return sum(b.n_nodes() for b in self._boundaries
                   if b.btype == btype)

    def add_boundary(self, nodes, btype, comment=None):
        """ Add one boundary.

            Parameters
            ----------
            nodes: an array of integer
                array of inter node indexes in the boundary path

            btype: integer
                an integer constant for boundary types.

            comment: optional, string
        """
        nodes = self._rearrange_boundary_nodes_in_ccw(nodes)
        # Call the super class method first
        super(SchismMesh, self).add_boundary(nodes, btype)
        self._boundaries.append(SchismBoundary(nodes, btype, comment))

    def clear_boundaries(self):
        """ Delete all the current boundary information
        """
        del self._boundaries[:]
        self._clear_edge_types()

    def _check_if_beginning_of_boundary(self, node_i, btype=None):
        """ Check if the given node is any of starting nodes of boundary
            node strings.
        """
        return any(
            btype is None
            and node_i == boundary.nodes[0]
            or btype is not None
            and boundary.btype == btype
            and node_i == boundary.nodes[0]
            for boundary in self._boundaries
        )

    def _get_next_node_on_boundary_and_remove_edge(self, node_i,
                                                   not_assigned,
                                                   ccw=True):
        edges_i = self.get_edges_from_node(node_i)
        for edge_i in edges_i:
            edge = self._edges[edge_i]
            if ccw:
                if edge[2] != EdgeType.INTERNAL and edge[0] == node_i:
                    try:
                        not_assigned.remove(edge_i)
                        return edge[1]
                    except Exception as exc:
                        print(
                            f"Attempted to remove edge {edge_i} based on node {node_i} which is not in not_assigned list"
                        )
                        raise
            elif edge[2] != EdgeType.INTERNAL and edge[1] == node_i:
                not_assigned.remove(edge_i)
                return edge[0]
        return None

    def _fill_land_boundaries(self):
        """
        This function creates land boundaries that are adjacent to
        the open boundaries. The direction of the filling is CCW.
        """
        if len(self._boundaries) == 0:
            return
        for open_boundary in self._boundaries:
            last = open_boundary.nodes[-1]
                # Check if there is another boundary right next this
            if not self._check_if_beginning_of_boundary(last):
                ns = [last]
                done = False
                while not done:
                    next_node = self._get_next_node_on_boundary(last)
                    ns.append(next_node)
                    if self._check_if_beginning_of_boundary(next_node):
                        done = True
                        self.add_boundary(ns, BoundaryType.LAND)
                    else:
                        last = next_node

    def _fill_island_boundaries(self, not_assigned):
        """ This function fills missing island boundaries.
        """
        first_node_i = self._edges[not_assigned[0]][0]
        ns = [first_node_i]
        last_node_i = first_node_i
        done = False
        while not done:
            next_node = self._get_next_node_on_boundary_and_remove_edge(
                last_node_i, not_assigned)
            ns.append(next_node)
            if next_node == first_node_i:
                done = True
                self.add_boundary(ns, BoundaryType.ISLAND)
            else:
                last_node_i = next_node

    def _get_not_assigned_boundary_edges(self):
        """ Get all edges that are not assigned as boundaries
        """
        if self._edges is None:
            self.build_edges_from_elems()
        return [
            i for i, edge in enumerate(self._edges) if edge[2] == EdgeType.BOUNDARY
        ]

    def fill_land_and_island_boundaries(self):
        """ Fill land and island boundaries for boundary edges not assigned
        to boundaries.
        """
        self._fill_land_boundaries()
        not_assigned = self._get_not_assigned_boundary_edges()
        while True:
            if len(not_assigned) > 0:
                self._fill_island_boundaries(not_assigned)
            else:
                break

    def fill_open_boundaries(self):
        """ Fill open boundaries for boundary edges not assigned
        to boundaries.
        """
        if len(self._boundaries) == 0:
            raise ValueError('No boundaries are defined')
        for boundary in self._boundaries:
            last = boundary.nodes[-1]
                # Check if there is another boundary right next this
            if not self._check_if_beginning_of_boundary(last):
                ns = [last]
                while True:
                    next_node = self._get_next_node_on_boundary(last)
                    ns.append(next_node)
                    if self._check_if_beginning_of_boundary(next_node):
                        self.add_boundary(ns, BoundaryType.OPEN)
                        break
                    else:
                        last = next_node

    def find_two_neighboring_node_paths(self, line_segment):
        """
        Find two neighboring node paths around a line_segment

        Parameters
        ----------
        line_segment: array-like
            two end points of the line segment, start_x, start_y, end_x, end_y

        Returns
        ------
        up_path: array of int
            upstream side of node paths
        down_path: array of int
            downstream side of node paths
        """
        # Accessing parent class data directly
        if self._elem_index is None:
            self._build_elem_index()

        _line = np.array(line_segment).reshape(2, 2)
        p1 = _line[0, ]
        p2 = _line[1, ]

        # normal = np.array((x[1, 1] - x[0, 1], x[0, 0] - x[1, 0]))
        box = self._box_from_points(_line)
        hits = self._elem_index.intersection(box)
        # Test which one is actually intersect
        intersected_edges = []
        intersections = []
        for hit in hits:
            nodes_i = self.elem(hit)
            n_nodes = len(nodes_i)
            nodes = self._nodes[nodes_i, :2]
            for i in range(n_nodes):
                edges = np.roll(nodes, -i, axis=0)[:2]
                intersection = find_intersection(_line, edges)
                if intersection is not None:
                    edge = {nodes_i[i], nodes_i[(i + 1) % n_nodes]}
                    if edge not in intersected_edges:
                        intersected_edges.append(edge)
                        intersections.append(intersection)
        dists = [np.linalg.norm(x - p1) for x in intersections]
        edges_sorted = [intersected_edges[i] for i in np.argsort(dists)]

        up_path = []
        down_path = []
        tan = (p2 - p1)
        norm = np.array((-tan[1], tan[0]))
        for edge in edges_sorted:
            n1 = edge.pop()
            n2 = edge.pop()
            dir = np.dot(norm, self.nodes[n1, :2] - p1)
            if dir > 0.:
                if n1 not in down_path:
                    down_path.append(n1)
                if n2 not in up_path:
                    up_path.append(n2)
            else:
                if n1 not in up_path:
                    up_path.append(n1)
                if n2 not in down_path:
                    down_path.append(n2)

        # Check if node strings are continuous and if not fix it
        # Up to quad elements are assumed.
        up_path_updated = deepcopy(up_path)
        n_inserted_nodes = 0
        for i in range(len(up_path) - 1):
            edge = self.find_edge((up_path[i], up_path[i + 1]))
            if edge is None:
                elem_up = self.get_elems_i_from_node(up_path[i]).intersection(
                    self.get_elems_i_from_node(up_path[i + 1]))
                nodes_in_elem = self.elem(elem_up.pop())
                for node_i in nodes_in_elem:
                    if node_i != up_path[i] and node_i != up_path[i + 1] and node_i not in down_path:
                        up_path_updated.insert(
                            i + 1 + n_inserted_nodes, node_i)
                        n_inserted_nodes += 1
                        break

        down_path_updated = deepcopy(down_path)
        n_inserted_nodes = 0
        for i in range(len(down_path) - 1):
            edge = self.find_edge((down_path[i], down_path[i + 1]))
            if edge is None:
                elem_up = self.get_elems_i_from_node(down_path[i]).intersection(
                    self.get_elems_i_from_node(down_path[i + 1]))
                nodes_in_elem = self.elem(elem_up.pop())
                for node_i in nodes_in_elem:
                    if node_i != down_path[i] and node_i != down_path[i + 1] and node_i not in up_path:
                        down_path_updated.insert(
                            i + 1 + n_inserted_nodes, node_i)
                        n_inserted_nodes += 1
                        break

        return up_path_updated, down_path_updated

    def _order_up_nodes_from_point(self, nodes, x):
        """ Order up nodes based on the distance from x

            Parameters
            ----------
            nodes:
                list of node indexes
        """
        dist = self._distance(nodes, x)
        sorted_indexes = np.argsort(dist)
        return [nodes[i] for i in sorted_indexes]

    def _distance(self, nodes_i, x):
        """ Calculate a distance from a node[node_i] to a point x(x1, x2)

            Returns
            -------
            float
                distance
        """
        nodes = tuple(self._nodes[i] for i in nodes_i)
        diffs = tuple(np.subtract(nodes[i][:2], x)
                      for i in range(len(nodes)))
        return tuple(np.linalg.norm(diff) for diff in diffs)

    def _rearrange_boundary_nodes_in_ccw(self, nodes):
        """ Make sure the boundary nodes are in CCW ordering.
            Assume that nodes are in row at least.
            nodes = the list of node indexes
            return = reordered list of the node indexes
        """
        edge_i = self.find_edge(nodes[:2], True)
        if edge_i is not None:
            return nodes
        new_nodes = list(reversed(nodes))
        del nodes
        return new_nodes

    def trim_to_left_of_mesh(self, line_segments):
        """ Trim mesh using line_segments.
            This function trims the mesh on the left sides
            of the line segments. The left side here means left when you look
            at the second end point of a line segment from the first one.
            An actual path to trimming is a nodal path that is on the
            right side of the line segment.
            To manage torus like mesh topology, the argument takes an array
            of line segments. The user need to provides a set of line segments
            to make sure the left side of the line segments does not
            cover the whole mesh.
            To trim multiple parts of a mesh, use this function
            multiple times.
            This function clears up any boundary definitions associated with
            the original grid. It is user's responsibility to create them
            again.

            line_segments = array of line segments defined by two end points
        """
        # self._logger.info("Trimming the mesh...")
        paths = []
        for l in line_segments:
            p = self.find_two_neighboring_node_paths(l)
            paths.append(p[0])
        self.trim_elems(paths)
        self.clear_boundaries()
        # self._logger.info("Removed the old boundary information.")

    def build_z(self, elev=0.):
        """ Build vertical coordinates

            Returns
            -------
            numpy.ndarray
        """
        if self.vmesh is None:
            raise ValueError("No vertical mesh information")
        # NOTE: Possibility of cyclic reference
        return self.vmesh.build_z(self, elev)

    @property
    def z(self):
        """ Get the elevation of levels at all the nodes

            Returns
            -------
            numpy.ndarray
                Depth array of each levels at the nodes.
                The shape of the array is (n_nodes, n_vert_levels)
                The levels below the most bottom level will be filled with
                the elevation of the bottom level
        """
        if self.vmesh._z is None:
            self.build_z()
        return self.vmesh._z

    @property
    def n_vert_levels(self):
        if self._vmesh is not None:
            return self._vmesh.n_vert_levels()
        else:
            raise ValueError("No vertical mesh is provided")

    def get_coordinates_3d(self):
        """ Get the array of each 3-d nodes

            Returns
            -------
            numpy.ndarray
                Array of 3-d coordinates
                The shape of the array is (n_nodes * n_vert_levels, 3)
                The points below the most bottom level will be filled with
                the elevation of the the bottom level
        """
        n_nodes = self.n_nodes()
        n_vert_levels = self._vmesh.n_vert_levels()
        points = np.empty((n_nodes * n_vert_levels, 3))
        for i in range(n_nodes):
            z = self.z[i]
            for j in range(n_vert_levels):
                r = i * n_vert_levels + j
                points[r, :2] = self.nodes[i, :2]
                points[r, 2] = z[j]
        return points

    def get_centers_of_elements(self):
        """ Get the array of the centers of the element

            Returns
            -------
            numpy.ndarray
                Array of 2-d coordinates of the centers of the elements
                The shape of the array is (n_elems, 2)
        """
        centers = np.empty((self.n_elems(), 2))
        for i, elem in enumerate(self.elems):
            centers[i] = np.mean(self.nodes[elem, :2], axis=0)
        return centers

    def get_centers_of_sides(self):
        """ Get the array of the centers of the sides

            Returns
            -------
            numpy.ndarray
                Array of 3-d coordinates of the centers of the elements
                The shape of the array is (n_elems, 2)
        """
        centers = np.empty((self.n_edges(), 3))
        for i, edge in enumerate(self.edges):
            centers[i, :2] = np.mean(self.nodes[self.edges[i, :2], :2], axis=0)
        return centers

    def areas(self):
        """ Get the array of element areas

            Returns
            -------
            numpy.ndarray
        """
        if self._areas is None:
            self._calculate_areas()
        return self._areas

    def _create_2d_points(self):
        self._points = [Point(n[0], n[1]) for n in self._nodes]

    def _create_2d_polygons(self):
        n = self._nodes
        self._polygons = [Polygon(n[e, :]) for e in self.elems]

    def _calculate_areas(self):
        """ Calculate the areas of the elements
        """
        if self._polygons is None:
            self._create_2d_polygons()
        self._areas = np.array([p.area for p in self._polygons])

    def _calculate_edge_lens(self):
        """ Calculate the lengths of the edges
        """
        edges = self.edges
        self._edge_lengths = np.linalg.norm(
            self._nodes[edges[:, 0], :2] - self._nodes[edges[:, 1], :2],
            axis=1)

    def _calculate_centroids(self):
        """ Calculate the centroids of all the elements
            NOTE: Assume that the polygons are simple polygons
            that do not consist of multiple parts.
        """
        if self._polygons is None:
            self._create_2d_polygons()
        self._centroids = np.array(
            [p.centroid.coords[0] for p in self._polygons])

    def edge_len(self):
        """ Get a NumPy array of edge lengths
            The ordering of the edges is decided from triquadmesh.
            Look edges properties to find out connectivity

            Returns
            -------
            numpy.ndarray
        """
        if self._edge_lengths is None:
            self._calculate_edge_lens()
        return self._edge_lengths

    def centroids(self):
        """ Get a Numpy array of element centroids
            Returns
            -------
            numpy.ndarray
        """
        if self._centroids is None:
            self._calculate_centroids()
        return self._centroids

    def merged_mesh(self):
        """
        merge all the meshes to create a boundary polygon
        """
        if self._merged_mesh is None:
            from shapely.ops import cascaded_union
            self._create_2d_polygons()
            self._merged_mesh = cascaded_union(self._polygons)
        return self._merged_mesh

    def to_geopandas(self,feature_type='polygon',crs=None,shp_fn=None,
                     node_values=None,elem_values=None,edge_values=None,
                     value_name=None,create_gdf=True):
        """
        Create mesh polygon and points as geopandas dataframe and shapefiles if
        Shap_fn file name is supplied.
        """
        import geopandas as gpd
        import pandas as pd
        df = pd.DataFrame()
        if feature_type == 'polygon':
            self._create_2d_polygons()
            features = self._polygons
            if elem_values is not None:
                df[value_name] = elem_values
        elif feature_type == 'point':
            self._create_2d_points()
            features = self._points
            if node_values is not None:
                df[value_name] = node_values
        elif feature_type == 'edge':
            features = self.get_centers_of_sides()[:, :2]
            if edge_values is not None:
                df[value_name] = edge_values
        df['geometry'] = features
        gdf = gpd.GeoDataFrame(df,geometry='geometry')

        gdf.crs = crs or None
        if shp_fn:
            gdf.to_file(shp_fn)
        if create_gdf:
            return gdf

    def plot_elems(self,var=None,ax=None,inpoly=None,plot_nan=False, **kwargs):
        """
        Plot variables (1D array on element) based on SCHISM mesh grid.
        if var is None,just plot the grid.
        """
        xy = [self._nodes[e, :2] for e in self.elems]
        if len(xy) != len(var):
            raise ValueError("input var has different len compared to input elem")
        if inpoly is not None:
            xy = np.asarray(xy)[inpoly]
            if var is not None and not plot_nan:
                vpoly = np.isnan(var) # do not plot nan data
                inpoly = np.logical_and(inpoly, ~vpoly)
            coll = PolyCollection(xy,array=var[inpoly],**kwargs)
        elif plot_nan:
            coll = PolyCollection(xy,array=var,**kwargs)
        else:
            inpoly = ~np.isnan(var)
            xy = np.asarray(xy)[inpoly]
            coll = PolyCollection(xy,array=var[inpoly],**kwargs)
        if not ax:
            fig, ax = plt.subplots()
        ax.add_collection(coll)
        ax.axis('equal')
        return coll

    def plot_nodes(self,var,ax=None,inpoly=None,plot_nan=False,**kwargs):
        if len(self.nodes) != len(var):
            raise ValueError("input var has different len compared to input node")
        velem = np.asarray([var[el].mean(axis=0) for el in self.elems])
        if inpoly is not None:
            inpoly = np.asarray([np.all(inpoly[el]) for el in self.elems])
        return self.plot_elems(
            var=velem, ax=ax, inpoly=inpoly, plot_nan=plot_nan, **kwargs
        )

    def plot_edges(self,var,ax=None, size=500,inpoly=None,**kwargs):
        xy = np.array(self.get_centers_of_sides())
        if len(xy) != len(var):
            raise ValueError("input var has different len compared to input edge")
        if not ax:
            fig, ax = plt.subplots()
        if inpoly is not None:
            xy = self.get_centers_of_sides()[inpoly][:2]
            return ax.scatter(xy[:, 0], xy[:, 1], c=var[inpoly], s=size, **kwargs)
        else:
            xy = xy[:,:2]
            return ax.scatter(xy[:,0], xy[:,1], c=var, s=size,**kwargs)

    def plot_mesh_boundary(self,ax=None,edgecolor='k',**kwargs):
        """
        Plot the edge of the computational grid.
        """
        boundary_mesh = self.merged_mesh()
        patch = []
        if isinstance(boundary_mesh, Polygon):
            multi_poly = [boundary_mesh]
        else:
            multi_poly = list(boundary_mesh)
        for pi in multi_poly:
            if isinstance(pi.boundary, MultiLineString):
                boundary = list(pi.boundary)
                for b in boundary:
                    cxy = np.asarray(b.xy).T
                    poly = patches.Polygon(cxy,True)
                    patch.append(poly)
            else:
                cxy = np.asarray(pi.boundary.xy).T
                poly = patches.Polygon(cxy,True)
                patch.append(poly)

        p = PatchCollection(patch,facecolor='none',edgecolor=edgecolor,**kwargs)
        if not ax:
            fig, ax = plt.subplots()
            bounds = boundary_mesh.bounds
            ax.set_xlim((bounds[0],bounds[2]))
            ax.set_ylim((bounds[1],bounds[3]))
        ax.add_collection(p)
        plt.axis('equal')
        return p

class SchismMeshReader(object):
    """ Schism Mesh Reader abstract class
    """

    def read(self, *args, **kwargs):
        raise NotImplementedError()


class SchismMeshWriter(object):
    """ Schism Mesh Writer abstract class
    """

    def write(self, *args, **kwargs):
        raise NotImplementedError()


class SchismMeshGr3Reader(SchismMeshReader):
    """ Read a mesh from a GR3 file
    """

    def __init__(self, logger=None):
        """ Constructor
        """
        self._logger = logger
        self._mesh = None

    def read_header(self, f):
        # First line
        first_line = f.readline()
        self._mesh._metadata['name'] = first_line.strip()
        # Second line: # of elements and # of nodes
        tkns = f.readline().split()
        if len(tkns) < 2:
            raise ValueError('Not enough items is the second line')
        n_elems, n_nodes = list(map(int, tkns[:2]))
        if n_elems < 1 or n_nodes < 1:
            raise ValueError('Invalid # of elems or nodes')
        self._mesh.allocate(n_elems, n_nodes)  # Allocate memory

    def read_nodes(self, f):
        for node_counter, _ in enumerate(range(self._mesh.n_nodes())):
            line = f.readline()
            tkns = line.split()
            if len(tkns) < 4:
                self._logger.error("Error reading: %s", line)
                raise ValueError("Node block is corrupt.")
            node_coords = list(map(float, tkns[1:4]))
            self._mesh.set_node(node_counter, node_coords)

    def read_elems(self, f):
        for elem_i in range(self._mesh.n_elems()):
            line = f.readline()
            tkn = line.split()
            if len(tkn) < 5:
                self._logger.error("Error reading: %s", line)
                raise ValueError("Element block is corrupt")
            type_elem = int(tkn[1])
            if type_elem < 3 or type_elem > 4:
                self._logger.error("Error reading: %s", line)
                self._logger.error(
                    "Only triangular or quadrilateral are supported.")
                raise ValueError("Element block is corrupt")
            # Zero-based connectivities
            try:
                if type_elem == 3:
                    connectivities = np.subtract(
                        np.array(list(map(int, tkn[2:5]))), 1)
                elif type_elem == 4:
                    connectivities = np.subtract(
                        np.array(list(map(int, tkn[2:6]))), 1)
                else:
                    self._logger.error("Error reading: %s", line)
                    raise ValueError("Element block is corrupt")
            except ValueError:
                self._logger.error("Error reading: %s", line)
                raise ValueError("Element block is corrupt")
            self._mesh.set_elem(elem_i, connectivities)

    def read_boundaries(self, f):
        """ Read boundary information
        """
        # Open boundaries
        # # of open boundaries
        line = f.readline()
        tkns = line.split()
        if len(tkns) < 1:
            if self._logger is not None:
                self._logger.warning("No boundary information is present?")
            return
        try:
            n_open_boundaries = int(tkns[0])
            line = f.readline()
            tkns = line.split()
            if len(tkns) < 1:
                self._logger.error("Error reading: %s", line)
                raise ValueError("Boundary block is corrupt")
            n_open_boundary_nodes = int(tkns[0])
            if len(tkns) < 1:
                self._logger.error("Error reading: %s", line)
                raise ValueError("Boundary block is corrupt")
            for _ in range(n_open_boundaries):
                line = f.readline()
                tkns = line.split()
                if len(tkns) < 1:
                    self._logger.error("Error reading: %s", line)
                    raise ValueError("Boundary block is corrupt")
                n_nodes = int(tkns[0])
                nodes = []
                for _ in range(n_nodes):
                    line = f.readline()
                    tkns = line.split()
                    if len(tkns) < 1:
                        self._logger.error("Error reading: %s", line)
                        raise ValueError("Boundary block is corrupt")
                    node = int(tkns[0]) - 1  # Zero based
                    nodes.append(node)
                self._mesh.add_boundary(nodes,
                                        BoundaryType.OPEN)

            line = f.readline()
            tkns = line.split()
            if len(tkns) < 1:
                self._logger.error("Error reading: %s", line)
                raise ValueError("Boundary block is corrupt")
            n_land_boundaries = int(tkns[0])
            line = f.readline()
            tkns = line.split()
            if len(tkns) < 1:
                self._logger.error("Error reading: %s", line)
                raise ValueError("Boundary block is corrupt")
            n_land_boundary_nodes = int(tkns[0])
            for _ in range(n_land_boundaries):
                line = f.readline()
                tkns = line.split()
                if len(tkns) < 1:
                    self._logger.error("Error reading: %s", line)
                    raise ValueError("Boundary block is corrupt")
                n_nodes = int(tkns[0])
                nodes = []
                for _ in range(n_nodes):
                    line = f.readline()
                    tkns = line.split()
                    if len(tkns) < 1:
                        self._logger.error("Error reading: %s", line)
                        raise ValueError("Boundary block is corrupt")
                    node = int(tkns[0]) - 1  # Zero based
                    nodes.append(node)
                self._mesh.add_boundary(nodes,
                                        BoundaryType.LAND)
        except ValueError:
            raise ValueError("Boundary block is corrupt")

    def read(self, *args, **kwargs):
        """ Read in a hgrid.gr3 file.
            If mode is 1, it does not read in boundary information.

            Returns
            -------
            SchismMesh
        """
        fpath_mesh = kwargs.get('fpath_mesh', 'hgrid.gr3')
        read_boundary = kwargs.get('read_boundary', True)
        for i, arg in enumerate(args):
            if i == 0:
                fpath_mesh = arg
            elif i == 1:
                read_boundary = arg

        if not os.path.exists(fpath_mesh):
            raise ValueError(f'File not found:{fpath_mesh}')
        if self._logger is not None:
            self._logger.debug(f"Reading in a gr3 file: {fpath_mesh} ...")

        self._mesh = SchismMesh()
        # Horizontal mesh
        with open(fpath_mesh, 'r') as f:
            self.read_header(f)
            self.read_nodes(f)
            self.read_elems(f)

            # Boundary info
            if read_boundary:
                self.read_boundaries(f)
        return self._mesh


class SchismMeshSmsReader(SchismMeshReader):
    """ Read a mesh from a 2dm file
    """

    def __init__(self, logger=None):
        """ Constructor
        """
        self._logger = logger
        self.mesh = SchismMesh()

    def read_header(self, fpath):
        with open(fpath, 'r') as file_object:
            line = file_object.readline()
            row = line.strip().split()
            if row[0] != 'MESH2D':
                raise ValueError("Not supported 2dm format")
            line = file_object.readline()
            row = line.strip().split()
            if row[0] == 'MESHNAME':
                name = line.split(' ', 1)[1].strip().strip('"')
                self.mesh._metadata['name'] = name
            else:
                raise ValueError("No mesh name in 2dm file")
            row = file_object.readline().strip().split()

    def prepare_mesh(self, fpath):
        n_nodes = self.read_n_nodes(fpath)
        n_elems = self.read_n_elements(fpath)
        self.mesh.allocate(n_elems, n_nodes)  # Allocate memory

    def read_n_elements(self, fpath):
        n_elems = 0
        with open(fpath, 'r') as file_object:
            for line in file_object:
                if len(line.strip()) < 1:
                    continue
                row = line.strip().split()
                if row[0] in ('E3T', 'E4Q'):
                    n_elems += 1
        return n_elems

    def read_n_nodes(self, fpath):
        n_nodes = 0
        with open(fpath, 'r') as file_object:
            for line in file_object:
                if len(line.strip()) < 1:
                    continue
                row = line.strip().split()
                if row[0] == 'ND':
                    n_nodes += 1
        return n_nodes

    def read_elements(self, fpath):
        with open(fpath, 'r') as file_object:
            elem_i = 0
            for line in file_object:
                if len(line.strip()) < 1:
                    continue
                row = line.strip().split()
                if row[0] in ('E3T', 'E4Q'):
                    # elem_i = int(row[1]) - 1
                    if row[0] == 'E3T':
                        connectivities = np.subtract(
                            np.array(list(map(int, row[2:5]))), 1)
                    elif row[0] == 'E4Q':
                        connectivities = np.subtract(
                            np.array(list(map(int, row[2:6]))), 1)
                    else:
                        raise ValueError("Element block is corrupt")
                    self.mesh.set_elem(elem_i, connectivities)
                    elem_i += 1

    def read_nodes(self, fpath):
        with open(fpath, 'r') as file_object:
            for line in file_object:
                if len(line.strip()) < 1:
                    continue
                row = line.strip().split()
                if row[0] == 'ND':
                    node_i = int(row[1]) - 1
                    node_coords = list(map(float, row[2:5]))
                    try:
                        self.mesh.set_node(node_i, node_coords)
                    except:
                        print("Bad node: . This can be caused by deleting nodes in SMS and failing to renumber:".format(node_i))
                        raise

    def read_boundary(self, fpath, nodestring_option=None):
        """ Read in nodestrings as boundaries
        """
        if nodestring_option is None:
            return
        with open(fpath, 'r') as file_object:
            bound_nodes = []
            for line in file_object:
                row = line.strip().split()
                if row[0] == 'NS':
                    nodes = list(map(int, row[1:]))
                    if nodes[0] < 0:
                        bound_nodes.append(-nodes[0] - 1)
                        if nodestring_option == 'open':
                            self.mesh.add_boundary(
                                bound_nodes, BoundaryType.OPEN)
                        elif nodestring_option == 'land':
                            self.mesh.add_boundary(
                                bound_nodes, BoundaryType.LAND)
                        else:
                            print('boundary_assignment', nodestring_option)
                            raise ValueError('Not supported option')
                        bound_nodes = []
                    else:
                        bound_nodes.extend([x - 1 for x in nodes])
        if nodestring_option == 'open':
            self.mesh.fill_land_and_island_boundaries()
        else:
            self.mesh.fill_open_boundaries()

    def read(self, *args, **kwargs):
        fpath_mesh = kwargs.get('fpath_mesh', 'hgrid.gr3')
        nodestring_option = kwargs.get('nodestring_option')
        for i, arg in enumerate(args):
            if i == 0:
                fpath_mesh = arg

        self.read_header(fpath_mesh)
        self.prepare_mesh(fpath_mesh)
        self.read_elements(fpath_mesh)
        self.read_nodes(fpath_mesh)
        self.read_boundary(fpath_mesh, nodestring_option)
        return self.mesh


class SchismMeshGr3Writer(SchismMeshWriter):

    def _write_boundary(self):
        mesh = self.mesh
        f = self.f
        # Open
        buf = "%d = Number of open boundaries\n" % mesh.n_boundaries(
            BoundaryType.OPEN)
        f.write(buf)
        buf = "%d = Total number of open boundary nodes\n" % mesh.n_total_boundary_nodes(
            BoundaryType.OPEN)
        f.write(buf)
        openbound_count = 0
        for bndry in mesh._boundaries:
            if bndry.btype == BoundaryType.OPEN:
                openbound_count += 1
                if bndry.comment is None:
                    buf = "%d = Number of nodes for open boundary %d\n" % \
                              (bndry.n_nodes(), openbound_count)
                else:
                    buf = "%d %s\n" % (bndry.n_nodes(), bndry.comment)
                f.write(buf)
                buf = "".join("%d\n" % (node_i + 1) for node_i in bndry.nodes)
                f.write(buf)
                # else:
                #     raise ValueError("Unsupported boundary type.")

        # Land & Island
        buf = "%d = Number of land boundaries\n" % (
            mesh.n_boundaries(BoundaryType.LAND) +
            mesh.n_boundaries(BoundaryType.ISLAND))
        f.write(buf)
        buf = "%d = Total number of land boundary nodes\n" % (
            mesh.n_total_boundary_nodes(BoundaryType.LAND) +
            mesh.n_total_boundary_nodes(BoundaryType.ISLAND))
        f.write(buf)
        landbound_count = 0
        for bndry in mesh._boundaries:
            if bndry.btype in [BoundaryType.LAND, BoundaryType.ISLAND]:
                landbound_count += 1
                island_flag = 0 if bndry.btype == BoundaryType.LAND else 1
                buf = (
                    "%d %d = Number of nodes for land boundary %d\n"
                    % (bndry.n_nodes(), island_flag, landbound_count)
                    if bndry.comment is None
                    else "%d %d %s\n"
                    % (bndry.n_nodes(), island_flag, bndry.comment)
                )
                f.write(buf)
                f.write("".join("%d\n" % (node_i + 1) for node_i in bndry.nodes))

    def write(self, *args, **kwargs):
        """
        Write a GR3 format grid.

        Parameters
        ----------
        mesh: SchismMesh
            a mesh to write
        fpath: str
            output file name
        node_attr: array-like, optional
            a list of node attribute. If provided, this will replace depth.
        boundary: boolean, optional
            default = True
            If true, boundary information will be written.
            Otherwise, it will not be appended.
        """
        mesh = kwargs.get('mesh')
        fpath = kwargs.get('fpath')
        node_attr = kwargs.get('node_attr')
        boundary = kwargs.get('write_boundary', True)
        for i, arg in enumerate(args):
            if i == 0:
                mesh = arg
            elif i == 1:
                fpath = arg
            elif i == 2:
                node_attr = arg
            elif i == 3:
                boundary = arg
        self.mesh = mesh

        with open(fpath, 'w') as f:
            self.f = f
            # Header
            name = mesh._metadata.get('name')
            if name is None or name == '':
                name = os.path.basename(fpath)
            buf = "%s\n" % name
            f.write(buf)

            n_elems = mesh.n_elems()
            n_nodes = mesh.n_nodes()
            buf = "%d %d ! # of elements and nodes\n" % (n_elems, n_nodes)
            f.write(buf)

            padding = 2
            maxnum = int(math.log10(max(n_elems, n_nodes))) + padding
            ifmtj = "%-" + ("%ii" % maxnum)
            ffmt = "%18.8f"
            nfmt = ifmtj + ffmt * 3 + "\n"
            # Nodes
            if node_attr is not None:
                for i in range(n_nodes):
                    buf = nfmt % (i + 1,
                                  mesh.nodes[i, 0],
                                  mesh.nodes[i, 1],
                                  node_attr[i])
                    f.write(buf)
            else:
                for i in range(n_nodes):
                    buf = nfmt % (i + 1,
                                  mesh.nodes[i, 0],
                                  mesh.nodes[i, 1],
                                  mesh.nodes[i, 2])
                    f.write(buf)

            # Elements
            for elem_i in range(n_elems):
                elem = mesh.elem(elem_i) + 1
                n_nodes_in_elem = len(elem)
                buf = '%d %d' % (elem_i + 1, n_nodes_in_elem)
                fmt = ' %d' * n_nodes_in_elem + '\n'
                buf += fmt % tuple(elem)
                f.write(buf)

            # Boundaries
            if boundary:
                self._write_boundary()


#class SchismMeshShapefileWriter(SchismMeshWriter):
#
#    def write(self,  *args, **kwargs):
#        """ Write a shapefile of the mesh
#            Support only UTM10N at the moment
#
#            Parameters
#            ----------
#            mesh: SchismMesh
#            fpath: str
#                output file name
#            node_attr: numpy.array
#                node values to overwrite
#        """
#        mesh = kwargs.get('mesh')
#        fpath = kwargs.get('fpath')
#        node_attr = kwargs.get('node_attr')
#        for i, arg in enumerate(args):
#            if i == 0:
#                mesh = arg
#            elif i == 1:
#                fpath = arg
#            elif i == 2:
#                node_attr = arg
#
#        if os.path.exists(fpath):
#            print("File with the output file name exists already", fpath)
#            raise RuntimeError("File exists already")
#
#        spatial_reference = osgeo.osr.SpatialReference()
#        crs = kwargs.get('crs')
#        if crs is None:
#            crs = 'EPSG:26910'
#        spatial_reference.ImportFromProj4(crs)
#        driver_name = 'ESRI Shapefile'
#        driver = osgeo.ogr.GetDriverByName(driver_name)
#        if driver is None:
#            print('%s is not available.' % driver_name)
#            raise RuntimeError()
#        datasource = driver.CreateDataSource(str(fpath))
#        if datasource is None:
#            raise RuntimeError()
#        layer = datasource.CreateLayer('mesh',
#                                       spatial_reference,
#                                       osgeo.ogr.wkbPolygon)
#        for i in range(4):
#            field = osgeo.ogr.FieldDefn('node%d' % i, osgeo.ogr.OFTInteger)
#            layer.CreateField(field)
#        field = osgeo.ogr.FieldDefn('cell', osgeo.ogr.OFTInteger)
#        layer.CreateField(field)
#        feature_defn = layer.GetLayerDefn()
#        feature = osgeo.ogr.Feature(feature_defn)
#        ring = osgeo.ogr.Geometry(osgeo.ogr.wkbLinearRing)
#        polygon = osgeo.ogr.Geometry(osgeo.ogr.wkbPolygon)
#
#        nodes = mesh.nodes
#        elems = mesh.elems
#        for cell_i, cell in enumerate(elems):
#            for i, node_index in enumerate(cell):
#                node = nodes[node_index]
#                ring.AddPoint(*node)
#                feature.SetField(i, int(node_index))
#            feature.SetField(4, int(cell_i))
#            ring.AddPoint(*nodes[cell[0]])
#            polygon.AddGeometry(ring)
#            feature.SetGeometry(polygon)
#            layer.CreateFeature(feature)
#            ring.Empty()
#            polygon.Empty()
#        datasource.Destroy()

class SchismMeshShapefileWriter(SchismMeshWriter):

    def write(self,  *args, **kwargs):
        """ Write a shapefile of the mesh
            Support only UTM10N at the moment

            Parameters
            ----------
            mesh: SchismMesh
            fpath: str
                output file name
            node_attr: numpy.array
                node values to overwrite
        """
        import pandas as pd
        import geopandas as gpd
        import logging

        mesh = kwargs.get('mesh')
        fpath = kwargs.get('fpath')
        node_attr = kwargs.get('node_attr')
        for i, arg in enumerate(args):
            if i == 0:
                mesh = arg
            elif i == 1:
                fpath = arg
            elif i == 2:
                node_attr = arg

        if os.path.exists(fpath):
            print("File with the output file name exists already", fpath)
            raise RuntimeError("File exists already")

        crs = kwargs.get('crs')
        if crs is None:
            crs = 'EPSG:26910'

        nodes = mesh.nodes
        node_values = nodes[:,2]
        value_name = os.path.basename(fpath).split('.')[0]

        # Two shape files will be generated, one for the polygons and the other
        # for the nodes.

        node_values = nodes[:, 2]
        value_name = os.path.basename(fpath).split('.')[0]

        # Two shape files will be generated, one for the polygons and the other
        # for the nodes.
        fpath_str = os.path.split(fpath)
        fdir = fpath_str[0]
        fname = fpath_str[1].replace('.shp', '')
        fpath_poly = os.path.join(fdir, f'{fname}_polygon.shp')
        fpath_point = os.path.join(fdir, f'{fname}_point.shp')

        mesh.to_geopandas('polygon', crs, fpath_poly, create_gdf=False)
        mesh.to_geopandas('point', crs, fpath_point,
                          node_values, value_name=value_name, create_gdf=False)
        logging.info(f"{fpath} generated")

class SchismMeshNetcdfWriter(SchismMeshWriter):

    def write(self,  *args, **kwargs):
        """ Write a netcdf of the mesh
            Support only UTM10N at the moment

            Parameters
            ----------
            mesh: SchismMesh generated by read_mesh(gr3_fn)
            fpath: str
                output file name
            node_attr: numpy.array
                node values to overwrite
        """
        import logging

        mesh = kwargs.get('mesh')
        fpath = kwargs.get('fpath')
        node_attr = kwargs.get('node_attr')

        for i, arg in enumerate(args):
            if i == 0:
                mesh = arg
            elif i == 1:
                fpath = arg
            elif i == 2:
                node_attr = arg

        if os.path.exists(fpath):
            print("File with the output file name exists already", fpath)
            raise RuntimeError("File exists already")

        nodes = mesh.nodes
        elems = mesh.elems
        edges = mesh.edges[:, 0:2]
        #n_node = range(mesh.n_nodes())
        #n_face = range(mesh.n_elems())
        #n_edge = range(mesh.n_edges())
        # if less than 4 nodes, give a negative value
        contour = np.zeros([len(elems), 4, 2])*-99.99
        for i in range(len(elems)):
            mesh_nodes = mesh.nodes[elems[i]][:, 0:2]
            # need to check if the mesh_nodes are consistenly in anti-clockwise direction
            contour[i, 0:len(mesh_nodes), :] = mesh_nodes
        depth = nodes[:, 2]
        value_name = os.path.basename(fpath).replace('.nc', '')
        n_face_node = range(4)  # maximum 4 nodes

        # the nc file generated may not exactly follow ugrid convention
        ds_node_x = xr.DataArray(nodes[:, 0],  # coords=[n_node],
                                 dims=['n_hgrid_node'],
                                 name='hgrid_node_x')
        ds_node_x.attrs['long_name'] = "node x-coordinate"
        ds_node_x.attrs['standard_name'] = 'longitude'
        # alternatively, degrees lat/lon
        ds_node_x.attrs['units'] = "degrees_east"
        ds_node_x.attrs['mesh'] = "SCHISM_hgrid"
        ds_node_x.attrs['location'] = "node"

        ds_node_y = xr.DataArray(nodes[:, 1],  # coords=[n_node],
                                 dims=['n_hgrid_node'],
                                 name='hgrid_node_y')
        ds_node_y.attrs['long_name'] = "node x-coordinate"
        ds_node_x.attrs['standard_name'] = 'latitude'
        # alternatively, degrees lat/lon
        ds_node_y.attrs['units'] = "degrees_north"
        ds_node_y.attrs['mesh'] = "SCHISM_hgrid"
        ds_node_y.attrs['location'] = "node"

        ds_values = xr.DataArray(depth,  # coords=[n_node],
                                 dims=['n_hgrid_node'],
                                 name='depth')
        ds_values.attrs['long_name'] = "Bathymetry"
        ds_values.attrs['units'] = "meters"
        ds_values.attrs['positive'] = "down"
        ds_values.attrs['mesh'] = "SCHISM_hgrid"
        ds_values.attrs['location'] = "node"

        ds_elem_contour_x = xr.DataArray(contour[:, :, 0],  # coords=[n_face,n_face_node],
                                         dims=['n_hgrid_face', 'n_face_node'], name='hgrid_contour_x')
        ds_elem_contour_x.attrs['long_name'] = "x_coordinates of 2D mesh contour"
        # alternatively, degrees lat/lon
        ds_elem_contour_x.attrs['units'] = "degree_east"
        ds_elem_contour_x.attrs['mesh'] = "SCHISM_hgrid"
        ds_elem_contour_x.attrs['location'] = "contour"

        ds_elem_contour_y = xr.DataArray(contour[:, :, 1],  # coords=[n_face,n_face_node],
                                         dims=['n_hgrid_face', 'n_face_node'], name='hgrid_contour_y')
        ds_elem_contour_y.attrs['long_name'] = "y_coordinates of 2D mesh contour"
        # alternatively, degrees lat/lon
        ds_elem_contour_y.attrs['units'] = "degree_north"
        ds_elem_contour_y.attrs['mesh'] = "SCHISM_hgrid"
        ds_elem_contour_y.attrs['location'] = "contour"

        # change from 0 based to 1 based.
        schism_face_node = np.asarray([np.append(e.astype('int32'), np.nan) if len(
            e) < 4 else e.astype('int32') for e in mesh.elems]) + 1
        ds_face_nodes = xr.DataArray(schism_face_node,  # coords=[n_face, n_face_node],
                                     dims=['n_hgrid_face', 'n_face_node'], name="hgrid_face_nodes").astype('int32')
        ds_face_nodes.attrs['long_name'] = "Horizontal Element Table"
        ds_face_nodes.attrs['cf_role'] = "face_node_connectivity"
        ds_face_nodes.attrs['start_index'] = 1

        # edge to nodes; change from 0 based to 1 based.
        schism_edge_nodes = mesh.edges[:, :2] + 1
        ds_edge_nodes = xr.DataArray(schism_edge_nodes,  # coords=[n_edge,np.arange(2)],
                                     dims=['n_hgrid_edge', 'two'], name='hgrid_edge_nodes').astype('int32')
        ds_edge_nodes.attrs['long_name'] = "Map every edge to the two nodes that it connects"
        ds_edge_nodes.attrs['cf_role'] = "edge_node_connectivity"
        ds_edge_nodes.attrs['start_index'] = 1

        face_x = np.array([np.mean(nodes[e, 0]) for e in elems])
        ds_face_x = xr.DataArray(face_x,  # coords=[n_face],
                                 dims=['n_hgrid_face'], name='hgrid_face_x')
        ds_face_x.attrs['long_name'] = "x_coordinate of 2D mesh face"
        ds_face_x.attrs['standard_name'] = "longitude"
        ds_face_x.attrs['units'] = "degrees_east"
        ds_face_x.attrs['mesh'] = "SCHISM_hgrid"

        face_y = np.array([np.mean(nodes[e, 1]) for e in elems])
        ds_face_y = xr.DataArray(face_y,  # coords=[n_face],
                                 dims=['n_hgrid_face'], name='hgrid_face_y')
        ds_face_y.attrs['long_name'] = "y_coordinate of 2D mesh face"
        ds_face_y.attrs['standard_name'] = "latitude"
        ds_face_y.attrs['units'] = "degrees_north"
        ds_face_y.attrs['mesh'] = "SCHISM_hgrid"

        edge_x = np.array([np.mean(nodes[e, 0]) for e in edges])
        ds_edge_x = xr.DataArray(edge_x,  # coords=[n_edge],
                                 dims=['n_hgrid_edge'], name='hgrid_edge_x')
        ds_edge_x.attrs['long_name'] = "x_coordinate of 2D mesh edge"
        ds_edge_x.attrs['standard_name'] = "longitude"
        ds_edge_x.attrs['units'] = "degrees_east"
        ds_edge_x.attrs['mesh'] = "SCHISM_hgrid"

        edge_y = np.array([np.mean(nodes[e, 1]) for e in edges])
        ds_edge_y = xr.DataArray(edge_y,  # coords=[n_edge],
                                 dims=['n_hgrid_edge'], name='hgrid_edge_y')
        ds_edge_y.attrs['long_name'] = "y_coordinate of 2D mesh edge"
        ds_edge_y.attrs['standard_name'] = "latitude"
        ds_edge_y.attrs['units'] = "degrees_north"
        ds_edge_y.attrs['mesh'] = "SCHISM_hgrid"

        ds = xr.merge([ds_node_x, ds_node_y, ds_elem_contour_x, ds_elem_contour_y, ds_values,
                       ds_face_nodes, ds_edge_nodes, ds_face_x, ds_face_y, ds_edge_x, ds_edge_y])

        if hasattr(mesh, 'vmesh'):
            vmesh = mesh.vmesh
            n_vert_levels = range(vmesh.n_vert_levels())
            node_bottom_index = vmesh.kbps + 1  # from 0 based to 1 based.
            ds_node_bottom = xr.DataArray(node_bottom_index,  # coords=[n_node],
                                            dims=['n_hgrid_node'], name='node_bottom_index')
            ds_node_bottom.attrs['long_name'] = "bottom level index at each node"
            ds_node_bottom.attrs['units'] = "non-dimensional"
            ds_node_bottom.attrs['mesh'] = "SCHISM_hgrid"
            ds_node_bottom.attrs['location'] = "node"
            ds_node_bottom.attrs['start_index'] = 1

            ele_bottom_index = np.array(
                [np.min(node_bottom_index[e]) for e in elems])
            ds_ele_bottom = xr.DataArray(ele_bottom_index,  # coords=[n_face],
                                            dims=['n_hgrid_face'], name='ele_bottom_index')
            ds_ele_bottom.attrs['long_name'] = "bottom level index at each element"
            ds_ele_bottom.attrs['units'] = "non-dimensional"
            ds_ele_bottom.attrs['mesh'] = "SCHISM_hgrid"
            ds_ele_bottom.attrs['location'] = "elem"
            ds_ele_bottom.attrs['start_index'] = 1

            edge_bottom_index = np.array(
                [np.min(node_bottom_index[e]) for e in edges[:, :2]])
            ds_edge_bottom = xr.DataArray(edge_bottom_index,  # coords=[n_edge],
                                            dims=['n_hgrid_edge'], name='edge_bottom_index')
            ds_edge_bottom.attrs['long_name'] = "bottom level index at each edge"
            ds_edge_bottom.attrs['units'] = "non-dimensional"
            ds_edge_bottom.attrs['mesh'] = "SCHISM_hgrid"
            ds_edge_bottom.attrs['location'] = "edge"
            ds_edge_bottom.attrs['start_index'] = 1

            z = mesh.z
            ds_z = xr.DataArray(z,  # coords=[n_node,n_vert_levels],
                                dims=['n_hgrid_node', 'n_vgrid_layers'], name='z')
            ds_z.attrs['mesh'] = 'SCHISM_hgrid'
            ds_z.attrs['data_horizontal_center'] = 'node'
            ds_z.attrs['data_vertical_center'] = 'full'
            ds_z.attrs['i23d'] = 2
            ds_z.attrs['ivs'] = vmesh.ivcor

            ds = xr.merge(
                [ds, ds_node_bottom, ds_ele_bottom, ds_edge_bottom, ds_z])

        ds.to_netcdf(fpath)
        logging.info(f"{fpath} generated")

class SchismMeshIoFactory(object):
    """ A factory class for SchismMeshIo
    """
    registered_readers = {'gr3': 'SchismMeshGr3Reader',
                          'sms': 'SchismMeshSmsReader'}
    registered_writers = {'gr3': 'SchismMeshGr3Writer',
                          'shp': 'SchismMeshShapefileWriter',
                          'nc': 'SchismMeshNetcdfWriter'}

    def get_reader(self, name):
        if name in self.registered_readers:
            return globals()[self.registered_readers[name]]()
        else:
            raise ValueError('Not in SchismMeshIoFactory')

    def show_registered_readers(self):
        print(self.registered_readers)

    def get_writer(self, name):
        if name in self.registered_writers:
            return globals()[self.registered_writers[name]]()
        else:
            raise ValueError('Not in SchismMeshIoFactory')


def read_mesh(fpath_mesh, fpath_vmesh=None, old_vgrid=True,**kwargs):
    """ Read a mesh data

        Returns
        -------
        SchismMesh
    """

    if fpath_mesh.endswith('.gr3'):
        reader = SchismMeshIoFactory().get_reader('gr3')
        mesh = reader.read(fpath_mesh)
        vmesh = read_vmesh(fpath_vmesh,old_vgrid) if fpath_vmesh is not None else None
        mesh._vmesh = vmesh
        if 'crs' in kwargs:
            mesh.crs = kwargs['crs']
        return mesh
    elif fpath_mesh.endswith('.2dm'):
        reader = SchismMeshIoFactory().get_reader('sms')
        mesh = reader.read(fpath_mesh, **kwargs)
        if 'crs' in kwargs:
            mesh.crs = kwargs['crs']
        return mesh
    else:
        raise ValueError("Unsupported extension")

def write_mesh(mesh, fpath_mesh, node_attr=None, write_boundary=False,
               crs=None, **kwargs):
    """ Write a horizontal mesh

        Parameters
        ----------
        mesh: SchismMesh
            Mesh
        fpath_mesh: str
    """
    if fpath_mesh.endswith(('gr3', 'll', 'ic')):
        writer = SchismMeshIoFactory().get_writer('gr3')
        writer.write(mesh=mesh,
                     fpath=fpath_mesh,
                     node_attr=node_attr,
                     write_boundary=write_boundary,
                     **kwargs)
    elif fpath_mesh.endswith(('shp')):
        writer = SchismMeshIoFactory().get_writer('shp')
        writer.write(mesh=mesh,
                     fpath=fpath_mesh,
                     node_attr=node_attr,
                     crs=crs,**kwargs)
    elif fpath_mesh.endswith(('nc')):
        writer = SchismMeshIoFactory().get_writer('nc')
        writer.write(mesh=mesh,
                     fpath=fpath_mesh,
                     node_attr=node_attr,**kwargs)
    else:
        raise ValueError("Unsupported extension")

def compare_mesh(mesh1, mesh2, dist_threshold=None):
    """
    For mesh2, find the nearest nodes in mesh1.
    Parameters
    ----------
    mesh1 : schism_mesh
        source mesh
        Generated by read_mesh(hgrid_fn)
    mesh2 : schism_mesh
        output_mesh
        Generated by read_mesh(hgrid_fn).

    Returns
    -------
    indices: indices of the nearest nodes in mesh1
    dist: distances between the nearest nodes

    """
    mesh1_idx = rtree.index.Rtree()
    for i, p in enumerate(mesh1.nodes[:,0:2]):
        mesh1_idx.insert(i, np.append(p,p))

    dist = []
    indices = []

    for n2 in mesh2.nodes[:,0:2]:
        idx = list(mesh1_idx.nearest(tuple(n2)))[0]
        dist2 = (n2[0]-mesh1.nodes[idx,0])**2 + (n2[1]-mesh1.nodes[idx,1])**2
        dist.append(np.sqrt(dist2))
        indices.append(idx)

    if dist_threshold is not None:
        dist = np.array(dist)
        print(
            f"The overlapping nodes account for {len(dist[dist <= dist_threshold]) / len(dist) * 100} percent of total nodes"
        )

    return indices, dist
