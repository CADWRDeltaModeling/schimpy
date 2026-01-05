import click
from schimpy.batch_metrics import batch_metrics_cli
from schimpy.clip_dems import clip_dems_cli
from schimpy.contour_smooth import contour_smooth_cli
from schimpy.convert_mesh import convert_mesh_cli
from schimpy.convert_polygons import convert_polygons_cli
from schimpy.convert_linestrings import convert_linestrings_cli
from schimpy.convert_points import convert_points_cli
from schimpy.combine_consume import combine_consume_cli
from schimpy.model_time import model_time_cli
from schimpy.sms2gr3 import convert_2dm_cli
from schimpy.prepare_schism import prepare_schism_cli
from schimpy.hotstart_inventory import hotstart_inventory_cli
from schimpy.create_vgrid_lsc2 import create_vgrid_lsc2_cli
from schimpy.schism_hotstart import create_hotstart_cli

# from schimpy.split_quad import split_quad
# from schimpy.small_areas import small_areas
# from schimpy.station import station
# from schimpy.nudging import create_nudging
# from schimpy.interpolate_structure import interpolate_structure
# from schimpy.merge_th import merge_th


@click.group(
    help="schimpy CLI tools for pre- and post- processing for SCHISM and data processing."
)
@click.help_option("-h", "--help")  # Add the help option at the group level
def cli():
    """Main entry point for schimpy commands."""
    pass


# Register the commands
cli.add_command(batch_metrics_cli, "batch_metrics")
cli.add_command(clip_dems_cli, "clip_dems")
cli.add_command(contour_smooth_cli, "contour_smooth")
cli.add_command(convert_mesh_cli, "convert_mesh")
cli.add_command(convert_polygons_cli, "convert_polygons")
cli.add_command(convert_linestrings_cli, "convert_linestrings")
cli.add_command(convert_points_cli, "convert_points")
cli.add_command(combine_consume_cli, "combine_consume")
cli.add_command(model_time_cli, "model_time")
cli.add_command(convert_2dm_cli, "convert_2dm")
cli.add_command(prepare_schism_cli, "prepare_schism")
cli.add_command(hotstart_inventory_cli, "hotstart_inventory")
cli.add_command(create_vgrid_lsc2_cli, "create_vgrid_lsc2")
cli.add_command(create_hotstart_cli, "schism_hotstart")

# cli.add_command(split_quad, "split_quad")
# cli.add_command(small_areas, "small_areas")
# cli.add_command(station, "station")
# cli.add_command(create_nudging, "create_nudging")
# cli.add_command(interpolate_structure, "interpolate_structure")
# cli.add_command(merge_th, "merge_th")


if __name__ == "__main__":
    cli()
