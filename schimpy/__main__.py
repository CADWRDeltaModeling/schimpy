import click
from schimpy.batch_metrics import batch_metrics_cli
from schimpy.clip_dems import main as clip_dems
from schimpy.contour_smooth import contour_smooth

# from schimpy.convert_mesh import convert_mesh
# from schimpy.convert_polygons import convert_polygons
# from schimpy.convert_linestrings import convert_linestrings
# from schimpy.convert_points import convert_points
# from schimpy.combine_consume import combine_consume
# from schimpy.prepare_schism import prepare_schism
# from schimpy.hotstart_inventory import hotstart_inventory
# from schimpy.create_vgrid_lsc2 import create_vgrid_lsc2
# from schimpy.schism_hotstart import schism_hotstart
# from schimpy.split_quad import split_quad
# from schimpy.model_time import model_time
# from schimpy.gen_elev2d import gen_elev2d
# from schimpy.small_areas import small_areas
# from schimpy.station import station
# from schimpy.nudging import create_nudging
# from schimpy.interpolate_structure import interpolate_structure
# from schimpy.merge_th import merge_th
# from schimpy.archive_ts import archive_ts


@click.group(
    help="schimpy CLI tools for pre- and post- processing for SCHISM and data processing."
)
@click.help_option("-h", "--help")  # Add the help option at the group level
def cli():
    """Main entry point for schimpy commands."""
    pass


# Register the commands
cli.add_command(batch_metrics_cli, "batch_metrics")
cli.add_command(clip_dems, "clip_dems")
cli.add_command(contour_smooth, "contour_smooth")
# cli.add_command(convert_mesh, "convert_mesh")
# cli.add_command(convert_polygons, "convert_polygons")
# cli.add_command(convert_linestrings, "convert_linestrings")
# cli.add_command(convert_points, "convert_points")
# cli.add_command(combine_consume, "combine_consume")
# cli.add_command(prepare_schism, "prepare_schism")
# cli.add_command(hotstart_inventory, "hotstart_inventory")
# cli.add_command(create_vgrid_lsc2, "create_vgrid_lsc2")
# cli.add_command(schism_hotstart, "schism_hotstart")
# cli.add_command(split_quad, "split_quad")
# cli.add_command(model_time, "model_time")
# cli.add_command(gen_elev2d, "gen_elev2d")
# cli.add_command(small_areas, "small_areas")
# cli.add_command(station, "station")
# cli.add_command(create_nudging, "create_nudging")
# cli.add_command(interpolate_structure, "interpolate_structure")
# cli.add_command(merge_th, "merge_th")
# cli.add_command(archive_ts, "archive_ts")


if __name__ == "__main__":
    cli()
