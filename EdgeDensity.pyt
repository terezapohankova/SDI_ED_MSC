# -*- coding: utf-8 -*-

import arcpy
import math



class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Edge Density"
        self.alias = "ED"

        # List of tool classes associated with this toolbox
        self.tools = [Tool]


class Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Density"
        self.description = "Calculation of Edge Density based on https://fragstats.org/index.php/fragstats-metrics/patch-based-metrics/area-and-edge-metrics/l4-edge-density"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
       
       
        # Input SHP layer with intersected grid and land patches
        param0 = arcpy.Parameter(
            displayName="Input Layer",
            name="input_layer",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Input")
        
     
        # Output SHP Layer
        param1 = arcpy.Parameter(
            displayName="Output Layer",
            name="out_features",
            datatype="DEShapefile",
            parameterType="Required",
            direction="Output")
        
    
        # Name of field with information about gridcell ID
        param2 = arcpy.Parameter(
            displayName="GRID ID",
            name="GRID_id",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        

        # Number of sides of grid (square -> 4, octagon ->8")
        param3 = arcpy.Parameter(
            displayName="Number of Grid Sides (4, 8 etc.)",
            name="grid_sides",
            datatype="GPLong",
            parameterType="Required",
            direction="Input")
        
        # Name of output Table 
        param4 = arcpy.Parameter(
            displayName="Name of Output Table",
            name="output_table",
            datatype="GPString",
            parameterType="Required",
            direction="Output")
        
        
        param1.parameterDependencies = [param0.name]
        param1.schema.clone = True
        
        params = [param0, param1, param2, param3, param4]
        
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
        
        # get text value of parameters
        INPUT_LAYER = parameters[0].valueAsText
        OUTPUT_LAYER = parameters[1].valueAsText
    
        GRID_ID = parameters[2].valueAsText
        GRID_SIDES = parameters[3].valueAsText
        
        OUTPUT_TABLE = parameters[4].valueAsText
        
        # Define new field names
        LAND_LENGTH_FIELD_NAME = "L_LAND"  # perimeter of landpatch (all sides)
        LAND_AREA_FIELD_NAME = "A_LAND"
        
        GRID_AREA_FIELD_NAME = "A_GRID"
        GRID_LENGTH_FIELD_NAME = "L_GRID" # perimeter of gridcell (all sides)
                
        ED_FIELD_NAME = "EDGE_DENS"
        
        # Define datype for fields
        FIELD_TYPE = "FLOAT"
      
        # Create copy of input layer to preserve original data
        input_layer_copy = arcpy.Copy_management(
            INPUT_LAYER, 
            OUTPUT_LAYER)
        
        # add field for calculation for perimeter of the land patch    
        arcpy.management.AddField(
            input_layer_copy, 
            LAND_LENGTH_FIELD_NAME, 
            FIELD_TYPE)
        
        # add field for calculation for gridcell area
        arcpy.management.AddField(
            input_layer_copy, 
            GRID_AREA_FIELD_NAME, 
            FIELD_TYPE)
        
        # add field for calculation for land patch area
        arcpy.management.AddField(
            input_layer_copy, 
            LAND_AREA_FIELD_NAME, 
            FIELD_TYPE)
        
        # add field for calculation of edge density
        arcpy.management.AddField(
            input_layer_copy, 
            ED_FIELD_NAME, 
            FIELD_TYPE)
        
        # add field for calculation gridcell length
        arcpy.management.AddField(
            input_layer_copy, 
            GRID_LENGTH_FIELD_NAME, 
            FIELD_TYPE)
    
        # calculate area of land patch and write it into column
        arcpy.management.CalculateGeometryAttributes(
            input_layer_copy,
            [[LAND_LENGTH_FIELD_NAME, "PERIMETER_LENGTH"], [LAND_AREA_FIELD_NAME, "AREA"]],
            length_unit = "METERS",
            area_unit = "SQUARE_METERS"
        )
        
        # sum patch perimeter and patch area for each ID
        # write it into separate table and store in memory
        summary_table = r"in_memory/summary_table"
        arcpy.analysis.Statistics(
            OUTPUT_LAYER,
            summary_table,
            [[LAND_LENGTH_FIELD_NAME, "SUM"], [LAND_AREA_FIELD_NAME, "SUM"]],
            GRID_ID
        )
            
        # Join the memory-help summary statistics back to the output layer
        # based on Grid ID
        arcpy.management.JoinField(
            OUTPUT_LAYER,
            GRID_ID,
            summary_table,
            GRID_ID,
            [f"SUM_{LAND_LENGTH_FIELD_NAME}", f"SUM_{LAND_AREA_FIELD_NAME}"],
        )
        
        # calulate area of grid patches based on sum of areas of land patches (which are always cut by grid)
        # aka how much space land patches occupy in each gridcell
        # should be equal to area of griddcell (unless the are gaps in the layer) 
        arcpy.management.CalculateField(
           OUTPUT_LAYER,
           GRID_AREA_FIELD_NAME,
           expression=f"!SUM_{LAND_AREA_FIELD_NAME}!",
           expression_type="PYTHON"           
        )
        
        # calclaute perimeter of gridcell
        # calculate the square root of SUM_{LAND_AREA_FIELD_NAME} -> get the length of one side of the grid
        # multiply by the number of sides of the grid -> get the perimeter of the grid
        arcpy.management.CalculateField(
            OUTPUT_LAYER,
            GRID_LENGTH_FIELD_NAME,
            expression =   "math.sqrt(!{}!) * {}".format(f"SUM_{LAND_AREA_FIELD_NAME}", GRID_SIDES),
            expression_type = "PYTHON"        
        )
        
        # calulate EDGE INDEX
        # (total length of land patches - perimeter of grid cell) / grid area  
        arcpy.management.CalculateField(
           OUTPUT_LAYER,
           ED_FIELD_NAME,
           expression=f"(!SUM_{LAND_LENGTH_FIELD_NAME}! - !{GRID_LENGTH_FIELD_NAME}!) / !{GRID_AREA_FIELD_NAME}!",
           expression_type="PYTHON"           
        )
        
        # prepare to export table of ED
        # dissolve based on grid ID
        # table structure : ID Grid -> ED
        dissolve_table = r"in_memory/DISSOLVE"
        arcpy.management.Dissolve(
            OUTPUT_LAYER,
            dissolve_table,
            GRID_ID,
            [[ED_FIELD_NAME, "FIRST"]]
        )
        
        # export table in excel with the name from user parameter
        
        arcpy.conversion.TableToExcel(
            dissolve_table,
            f"{OUTPUT_TABLE}.xlsx"
        )
        
        return


