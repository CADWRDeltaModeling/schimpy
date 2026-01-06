#!/usr/bin/env python
import sys
import click
from osgeo import ogr


def read_keyfile(keyfile):
    """Reads a file pairing labels and values"""
    keys = {}
    with open(keyfile, "r") as kf:
        for line in kf:
            if line and len(line) > 2:
                key, val = line.strip().split()
                keys[key] = float(val)
    return keys


def create_poly(shapefile, dsetname, keyfile, polyfile, type, default=None):
    """Converts a polygon shapefile to the yaml input of the preprocessor"""
    keys = read_keyfile(keyfile)
    names = set()

    f = open(polyfile, "w")
    ds = ogr.Open(shapefile)
    if ds is None:
        print("Open failed.\n")
        sys.exit(1)
    print("Opening dataset: %s" % ds.GetLayer(0).GetName())

    lyr = ds.GetLayerByName(dsetname)

    lyr.ResetReading()
    if default:
        f.write("default: %s\n\n" % default)
    f.write("polygons:\n")

    for feat in lyr:
        feat_defn = lyr.GetLayerDefn()
        # name = "poly_" + feat.GetFieldAsString(3).strip()
        label = feat.GetFieldAsString(2).strip()
        print("Name: %s" % (label))

        # print "mindepth %s" %rough
        # for i in range(feat_defn.GetFieldCount()):
        # field_defn = feat_defn.GetFieldDefn(i)

        # # Tests below can be simplified with just :
        # # print feat.GetField(i)
        # if field_defn.GetType() == ogr.OFTInteger:
        # print "int %d" % feat.GetFieldAsInteger(i)
        # elif field_defn.GetType() == ogr.OFTReal:
        # print "real %.3f" % feat.GetFieldAsDouble(i)
        # elif field_defn.GetType() == ogr.OFTString:
        # print "string %s" % feat.GetFieldAsString(i)
        # else:
        # print "string %s" % feat.GetFieldAsString(i)

        geom = feat.GetGeometryRef()
        if geom is not None and geom.GetGeometryType() == ogr.wkbPolygon:
            val = keys[label]
            f.write("  %s:\n" % label)
            f.write("    attribute: %s \n" % val)
            f.write("    type: %s\n" % type)
            f.write("    vertices:\n")
            for ring in geom:
                # f.write("%s %s %s None\n" %(label,ring.GetPointCount(),val))
                points = ring.GetPointCount()
                for p in range(points):
                    x, y, z = ring.GetPoint(p)
                    f.write("      %s %s\n" % (x, y))
    ds = None
    f.close()


@click.command()
@click.argument("shapefile", type=click.Path(exists=True))
@click.argument("dset")
@click.option(
    "--keyfile",
    required=True,
    type=click.Path(exists=True),
    help="File mapping material property names to numerical values (space separated).",
)
@click.option(
    "--out",
    required=True,
    help="Output polygon file.",
)
@click.option(
    "--type",
    default="min",
    help="Polygon attribute type (min/max).",
)
@click.option(
    "--default",
    default=None,
    help="Global default in polygon specs.",
)
def material_poly_cli(shapefile, dset, keyfile, out, type, default):
    """Convert shapefile from SMS into polygon format for preprocessor.
    
    SHAPEFILE: Name of shapefile.
    
    DSET: Name of dataset inside shapefile.
    """
    create_poly(shapefile, dset, keyfile, out, type, default)


if __name__ == "__main__":
    material_poly_cli()
