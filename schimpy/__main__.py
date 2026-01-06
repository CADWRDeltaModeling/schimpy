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
from schimpy.split_quad import split_quads_cli
from schimpy.small_areas import small_areas_cli
from schimpy.station import convert_station_cli
from schimpy.nudging import create_nudge_cli
from schimpy.interpolate_structure import interpolate_structure_cli
from schimpy.merge_th import merge_th_cli
from schimpy.create_station_output import create_station_output_cli
from schimpy.param import set_param_cli
from schimpy.check_mesh_skewness import check_skewness_cli
from schimpy.combine_flux import combine_flux_cli
from schimpy.convert_sav_class_to_number import convert_sav_class_to_number_cli
from schimpy.create_mesh_n_levels import create_mesh_n_levels_cli
from schimpy.cruise import cruise_plot_cli
from schimpy.cut_mesh import cut_mesh_cli
from schimpy.download_hrrr import download_hrr_clip_cli
from schimpy.embed_raster import embed_raster_cli
from schimpy.schism_yaml import schism_yaml_cli
from schimpy.grid_opt import grid_opt_cli
from schimpy.material_poly import material_poly_cli
from schimpy.simulation_timing import simulation_timing_cli
from schimpy.stacked_dem_fill import stacked_fill_cli
from schimpy.mesh_volume_tvd import mesh_volume_tvd_cli
from schimpy.subset_schism_output import subset_schism_output_cli
from schimpy.scripts.archive_ts import archive_ts_cli


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
cli.add_command(split_quads_cli, "split_quad")
cli.add_command(small_areas_cli, "small_areas")
cli.add_command(convert_station_cli, "station")
cli.add_command(create_nudge_cli, "create_nudging")
cli.add_command(interpolate_structure_cli, "interpolate_structure")
cli.add_command(merge_th_cli, "merge_th")
cli.add_command(create_station_output_cli, "create_station_output")
cli.add_command(set_param_cli, "set_param")
cli.add_command(check_skewness_cli, "check_skewness")
cli.add_command(combine_flux_cli, "combine_flux")
cli.add_command(convert_sav_class_to_number_cli, "convert_sav_class_to_number")
cli.add_command(create_mesh_n_levels_cli, "create_mesh_n_levels")
cli.add_command(cruise_plot_cli, "cruise_plot")
cli.add_command(cut_mesh_cli, "cut_mesh")
cli.add_command(download_hrr_clip_cli, "download_hrrr")
cli.add_command(embed_raster_cli, "embed_raster")
cli.add_command(schism_yaml_cli, "schism_yaml")
cli.add_command(grid_opt_cli, "grid_opt")
cli.add_command(material_poly_cli, "material_poly")
cli.add_command(simulation_timing_cli, "simulation_timing")
cli.add_command(mesh_volume_tvd_cli, "mesh_volume_tvd")
cli.add_command(subset_schism_output_cli, "subset_schism_output")
cli.add_command(archive_ts_cli, "archive_ts")

if __name__ == "__main__":
    cli()
