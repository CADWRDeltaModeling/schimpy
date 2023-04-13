SCHISM hotstart

Examples
--------

Examples to create hotstart.nc for the following cases can be found 
in `BayDeltaSCHISM <https://github.com/CADWRDeltaModeling/BayDeltaSCHISM>`_
All test cases are under the directory `examples/hotstart/ <https://github.com/CADWRDeltaModeling/BayDeltaSCHISM/tree/master/examples/hotstart>`
All test cases start from the same date and has the same destination grid. 

========  ==========================================  =======================
Cases      Directory                                    Tracers
========  ==========================================  =======================
basic      examples/basic/                              temp, salt
hotstart   examples/hotstart_from_previous_hotstart     temp, salt
flooding   examples/flooded_island                      temp, salt
age        examples/tracer_age                          temp, salt, gen, age 
sed        examples/sed                                 temp, salt, sed
bio        examples/bio                                 temp, salt, cos or icm
========  ==========================================  =======================


Steps
-----    

	#. Edit yaml file "hotstart.yaml": date, hgrid_input_file, and vgrid_input_file. 
	#. Create domain division poplygons in a shapefile. 
	    For most applications, users can directly use hotstart_regions.shp in BayDeltaSCHISM/examples/hotstart/shapefile without any modification. 
        However, if your domain is larger than the area covered by the polygon, an error will be generated: "Orphaned nodes or cells found at ........"  
	    Most of time, polygon overlapping should not be allowed (allow_overlap=False), except in "island_flooding" test case where "allow_overlap" is set to True. 
	    Note that it is the responsibility of the users to create polygons that are non overlapping and complete. 
        Hint: enable "snapping options" in qgis (normally under "project" menu) to create non overlaying polygons or us SMS instead. 
	#. Prepare extrude_casts data for San Francisco Bay. 	
	    It is the user's responsiblity to create the cast data file
        An `example file <https://github.com/CADWRDeltaModeling/BayDeltaSCHISM/blob/master/examples/hotstart/data_in/polaris_transect_2021-04-20.csv>` 
        is provided to demonstrate the format of the data file. 
	    USGS data is available from \\cnrastore-bdo\Modeling_Data\usgs_cruise. 
	#. Decide if suisun_marsh region should use "extrude_casts" or "obs_points". 
        The general rule is that "extrude_casts" should be used when the USGS cruise data is available for the region. 
	#. Prepare input data from observations. 
	    It is the user's responsiblity to create the observational input files. 
        `Example input files <https://github.com/CADWRDeltaModeling/BayDeltaSCHISM/blob/master/examples/hotstart/data_in>` are provided to demonstrate the format of the file. 
	#. run python script "create_hotstart.py". 



