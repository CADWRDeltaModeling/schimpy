#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Line String data based on Shapely LineStrings"""
from . import schism_yaml
import shapely.geometry
from shapely.wkb import loads
from osgeo.osr import SpatialReference
from osgeo.ogr import (
    Open,
    GetDriverByName,
    FieldDefn,
    wkbLineString,
    Feature,
    CreateGeometryFromWkb,
)
import os
import copy


class LineString(shapely.geometry.LineString):
    """ """

    # Updated for Shapely 2 following this issue comment:
    # https://github.com/shapely/shapely/issues/1233#issuecomment-977837620
    _id_to_attrs = {}
    __slots__ = shapely.geometry.LineString.__slots__

    def __new__(cls, coordinates=None, prop=None):
        p = super().__new__(cls, coordinates)
        p.__class__ = cls
        return p

    def __del__(self):
        del self._id_to_attrs[id(self)]

    def __init__(self, coordinates=None, prop=None):
        """
        Constructor

        Parameters
        ----------
        Returns
        -------
        """
        self._id_to_attrs[id(self)] = dict(prop=prop)

    @property
    def prop(self):
        return self._id_to_attrs[id(self)]["prop"]

    @prop.setter
    def prop(self, val):
        self._id_to_attrs[id(self)]["prop"] = val


class LineStringIo(object):
    def read(self, **kwargs):
        raise NotImplementedError()

    def write(self, **kwargs):
        raise NotImplementedError()


class LineStringYamlReader(LineStringIo):
    def read(self, fpath, **kwargs):
        with open(fpath, "r") as f:
            data = schism_yaml.load(f)["linestrings"]
            linestrings = []
            for row in data:
                linestrings.append(
                    LineString(
                        coordinates=row["coordinates"],
                        prop=dict([(k, row[k]) for k in row if k != "coordinates"]),
                    )
                )
            return linestrings


class LineStringShapefileReader(LineStringIo):
    def read(self, fpath, **kwargs):
        """
        Parameters
        ----------
        fpath: str
            input file name

        Returns
        -------
        lines
            list of LineStrings
        """
        if os.path.exists(fpath):
            datasource = Open(fpath)
            layer = datasource.GetLayer(0)
            feat = layer.GetFeature(0)
            field_names = [
                feat.GetFieldDefnRef(i).GetName() for i in range(feat.GetFieldCount())
            ]
            lines = []
            for feature in layer:
                geom = feature.GetGeometryRef()
                name_geom = geom.GetGeometryName()
                if name_geom in ("LINESTRING",):
                    line = LineString(
                        loads(bytes(geom.ExportToWkb())),
                        dict(
                            [
                                (k, feature.GetField(i))
                                for i, k in enumerate(field_names)
                            ]
                        ),
                    )
                    lines.append(line)
            return lines
        else:
            raise ValueError("File not found")


class LineStringYamlWriter(LineStringIo):
    """Write line strings from a YAML file"""

    def write(self, fpath, lines):
        with open(fpath, "w") as f:
            data = {}
            data["linestrings"] = []
            for line in lines:
                leaf = {}
                leaf["coordinates"] = list([list(xy) for xy in line.coords])
                for k in line.prop:
                    leaf[k] = line.prop[k]
                data["linestrings"].append(leaf)
            schism_yaml.dump(data, f)


class LineStringShapefileWriter(LineStringIo):
    def write(self, fpath, lines, spatial_reference=None, driver_name=None, **kwargs):
        """
        Parameters
        ----------
        fpath: str
            output file name
        lines: array of schism_linestring.LineString
            list of LineStrings
        spatial_reference: osgeo.osr.SpatialReference
        default: NAD83, UTM zone 10N, meter

        todo: reference needs to be in api right now hard to use
        """
        # Boilerplate to create a SHP file
        if spatial_reference is None:
            spatial_reference = SpatialReference()
            spatial_reference_str = (
                "+proj=utm +zone=10 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m"
                " +no_defs"
            )
            spatial_reference.ImportFromProj4(spatial_reference_str)
        if driver_name is None:
            driver_name = "ESRI Shapefile"
        driver = GetDriverByName(driver_name)
        if driver is None:
            print("%s is not available." % driver_name)
            raise RuntimeError()
        datasource = driver.CreateDataSource(fpath)
        if datasource is None:
            raise RuntimeError("Cannot create a GIS file")
        layer = datasource.CreateLayer("layer", spatial_reference, wkbLineString)
        fields = []
        for l in lines:
            if l.prop is not None:
                for k in l.prop:
                    if k not in fields and k != "coordinates":
                        fields.append(k)
        for field in fields:
            layer.CreateField(FieldDefn(field))
        feature_defn = layer.GetLayerDefn()
        feature = Feature(feature_defn)

        for i, line in enumerate(lines):
            feature.SetGeometry(CreateGeometryFromWkb(line.wkb))
            for j, f in enumerate(fields):
                val = line.prop.get(f)
                if val is not None:
                    feature.SetField(j, val)
            layer.CreateFeature(feature)
        datasource.Destroy()


class LineStringIoFactory(object):
    registered_readers = {
        "yaml": "LineStringYamlReader",
        "shp": "LineStringShapefileReader",
    }
    registered_writers = {
        "yaml": "LineStringYamlWriter",
        "shp": "LineStringShapefileWriter",
    }

    def get_reader(self, name):
        if name in self.registered_readers:
            return globals()[self.registered_readers[name]]()

    def get_writer(self, name):
        if name in self.registered_writers:
            return globals()[self.registered_writers[name]]()


def read_linestrings(fpath):
    fpath = str(fpath)  # Convert PosixPath to string
    if fpath.endswith(".yaml"):
        return LineStringIoFactory().get_reader("yaml").read(fpath)
    elif fpath.endswith(".shp"):
        return LineStringIoFactory().get_reader("shp").read(fpath)
    else:
        raise ValueError("Not supported file type")


def write_linestrings(fpath, lines):
    """
    Parameters
    ----------
    fpath: str
        output file name
    lines: array of schism_linestring.LineString
        list of LineStrings
    """
    fpath = str(fpath)  # Convert PosixPath to string
    if fpath.endswith(".yaml"):
        return LineStringIoFactory().get_writer("yaml").write(fpath, lines)
    if fpath.endswith(".shp"):
        return LineStringIoFactory().get_writer("shp").write(fpath, lines)
    else:
        raise ValueError("Not supported file type")
