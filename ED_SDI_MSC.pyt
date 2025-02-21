# -*- coding: utf-8 -*-

import arcpy
import math
import sys
import os


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Edge Density & Shanon Diversity Index & Mean Complixity Index"
        self.alias = "ED_SDI_MSC_"

        
        # List of tool classes associated with this toolbox
        self.tools = [Tool1, Tool2, Tool3]


class Tool1(object):
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
            datatype=["DEShapeFile","DEFeatureClass"],
            parameterType="Required",
            direction="Input"
            )
        
     
        # Output SHP Layer
        param1 = arcpy.Parameter(
            displayName="Output Layer",
            name="out_features",
            datatype=["DEShapeFile","DEFeatureClass"],
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
        """param4 = arcpy.Parameter(
            displayName="Path to Output Table Name",
            name="output_table",
            datatype="GPString",
            parameterType="Required",
            direction="Output")"""
        
        
        param1.parameterDependencies = [param0.name]
        param1.schema.clone = True
        
        params = [param0, param1, param2, param3]
        
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
        
        #OUTPUT_TABLE = parameters[4].valueAsText
        
        # Define new field names
        LAND_LENGTH_FIELD_NAME = "L_LAND"  # perimeter of landpatch (all sides)
        LAND_AREA_FIELD_NAME = "A_LAND"
        
        GRID_AREA_FIELD_NAME = "A_GRID"
        GRID_LENGTH_FIELD_NAME = "L_GRID" # perimeter of gridcell (all sides)
        
        #names for tool polygon neighbours to detect shared borders
        SHARED_BORDER_FIELD_NAME = "LENGTH"
        SHARED_BORDER_ID = "src_Id"
                
        ED_FIELD_NAME = "EDGE_DENS"
        
        # Define datype for fields
        FIELD_TYPE = "DOUBLE"
      
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
        
        # get length of shared border between the patches
        polygon_neighb =  r"in_memory/polygon_neib"

        
        arcpy.analysis.PolygonNeighbors(
            OUTPUT_LAYER,
            polygon_neighb,
            GRID_ID,
            area_overlap = "NO_AREA_OVERLAP",
            both_sides = "NO_BOTH_SIDES",
            out_area_units = "SQUARE_METERS"
        )

        

        arcpy.management.AlterField(
            polygon_neighb, 
            "src_" + GRID_ID, 
            GRID_ID,
            GRID_ID, 
            )

        arcpy.conversion.TableToExcel(
         polygon_neighb,
         "summaryy.xlsx"
        )

        arcpy.management.JoinField(
            OUTPUT_LAYER,
            GRID_ID,
            polygon_neighb,
            GRID_ID,
            [SHARED_BORDER_FIELD_NAME]
        )

        

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
        
        ## SUM_{LAND_LENGTH_FIELD_NAME} - sum of lengths of border in the grid 
        ## GRID_LENGTH_FIELD_NAME       - length of perimeter of the grid cell (for squares -> 4 * grid_length)
        ## SHARED_BORDER_FIELD_NAME     - length of inside borders between patches; length of shared border within one cell
        ## GRID_AREA_FIELD_NAME         - area of the grid cell
        arcpy.management.CalculateField(
           OUTPUT_LAYER,
           ED_FIELD_NAME,
           expression=f"(!SUM_{LAND_LENGTH_FIELD_NAME}! - !{GRID_LENGTH_FIELD_NAME}! -!{SHARED_BORDER_FIELD_NAME}!) / !{GRID_AREA_FIELD_NAME}!",
           expression_type="PYTHON"           
        )

        arcpy.management.DeleteField(
            OUTPUT_LAYER, 
            [LAND_LENGTH_FIELD_NAME, LAND_AREA_FIELD_NAME, GRID_AREA_FIELD_NAME, 
             GRID_LENGTH_FIELD_NAME, f"SUM_{LAND_AREA_FIELD_NAME}", f"SUM_{LAND_LENGTH_FIELD_NAME}", SHARED_BORDER_FIELD_NAME], 
            "DELETE_FIELDS")
        
        """# prepare to export table of ED
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
        ) """
        return
    
class Tool2(object): 
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shanon Diversity Index"
        self.description = "Calculation of Shanon Diversity Index on https://fragstats.org/index.php/fragstats-metrics/patch-based-metrics/diversity-metrics/l4-shannons-diversity-index"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
       
       
        # Input SHP layer with intersected grid and land patches
        param0 = arcpy.Parameter(
            displayName="Input Layer",
            name="input_layer",
            datatype=["DEShapeFile","DEFeatureClass"],
            parameterType="Required",
            direction="Input"
            )
     
        # Output SHP Layer
        param1 = arcpy.Parameter(
            displayName="Output Layer",
            name="out_features",
            datatype=["DEShapeFile","DEFeatureClass"],
            parameterType="Required",
            direction="Output")
        
    
        # Name of field with information about gridcell ID
        param2 = arcpy.Parameter(
            displayName="GRID ID",
            name="GRID_id",
            datatype="Field",
            parameterType="Required",
            direction="Input")
        
        
        param1.parameterDependencies = [param0.name]
        param1.schema.clone = True
        
        params = [param0, param1, param2]
        
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
        
        # Define new field names
        SHAPE_AREA_FIELD_NAME = "AREA_M"  # shape area of one patch in sq m
        LANDSCAPE_PROPORTION_FIELD_NAME = "pi" #proportion of landscape defined by each occurence
        LN_MULTIPLY_PROPORTION = "pilnpi"
        SHANON_INDEX = "SDI" 
        # Define datype for fields
        FIELD_TYPE = "FLOAT"
    
        # Create copy of input layer to preserve original data
        input_layer_copy = arcpy.Copy_management(
            INPUT_LAYER, 
            OUTPUT_LAYER)
        
        # add field for calculation for perimeter of the land patch    
        arcpy.management.AddField(
            input_layer_copy, 
            SHAPE_AREA_FIELD_NAME, 
            FIELD_TYPE)
        
        # add field for proportion of landscaoe by class 
        arcpy.management.AddField(
            input_layer_copy, 
            LANDSCAPE_PROPORTION_FIELD_NAME, 
            FIELD_TYPE)
        
        # add field for proportion of landscaoe by class 
        arcpy.management.AddField(
            input_layer_copy, 
            LN_MULTIPLY_PROPORTION, 
            FIELD_TYPE)
        
        # add field for proportion of landscaoe by class 
        arcpy.management.AddField(
            input_layer_copy, 
            SHANON_INDEX, 
            FIELD_TYPE)
        
        # calculate area of land patch and write it into column
        arcpy.management.CalculateGeometryAttributes(
            input_layer_copy,
            [[SHAPE_AREA_FIELD_NAME, "AREA"]],
            length_unit = "METERS",
            area_unit = "SQUARE_METERS"
        )
        
      # sum patch perimeter and patch area for each ID
        # write it into separate table and store in memory
        summary_table = r"in_memory/summary_table"
        arcpy.analysis.Statistics(
            OUTPUT_LAYER,
            summary_table,
            [[SHAPE_AREA_FIELD_NAME, "SUM"]],
            GRID_ID
        )
            
        # Join the memory-help summary statistics back to the output layer
        # based on Grid ID
        arcpy.management.JoinField(
            OUTPUT_LAYER,
            GRID_ID,
            summary_table,
            GRID_ID,
            [f"SUM_{SHAPE_AREA_FIELD_NAME}"],
        )
        
        arcpy.conversion.TableToExcel(
         summary_table,
         "summary.xlsx"
        )
        
        # calculate pi - a proportion of landscape class by dividing the are of each row by sum of all areas of id
        # id 173_a -> 5698
        # id 173_b -> 2569
        # id 173_c -> 2585
        # id 173 = sum(173_a + 173_b + 173_c)
        arcpy.management.CalculateField(
            OUTPUT_LAYER,
            LANDSCAPE_PROPORTION_FIELD_NAME,
            expression =   "!{}! / !{}!".format(SHAPE_AREA_FIELD_NAME ,f"SUM_{SHAPE_AREA_FIELD_NAME}"),
            expression_type = "PYTHON"        
        )
        
        arcpy.management.CalculateField(
            OUTPUT_LAYER,
            LN_MULTIPLY_PROPORTION,
            expression =   "!{}! * math.log(!{}!)".format(LANDSCAPE_PROPORTION_FIELD_NAME, LANDSCAPE_PROPORTION_FIELD_NAME),
            expression_type = "PYTHON"        
        )
        
         # sum patch perimeter and patch area for each ID
        # write it into separate table and store in memory
        summary_table_2 = r"in_memory/summary_table_2"
        arcpy.analysis.Statistics(
            OUTPUT_LAYER,
            summary_table_2,
            [[LN_MULTIPLY_PROPORTION, "SUM"]],
            GRID_ID
        )
        
        # Join the memory-help summary statistics back to the output layer
        # based on Grid ID
        arcpy.management.JoinField(
            OUTPUT_LAYER,
            GRID_ID,
            summary_table_2,
            GRID_ID,
            [f"SUM_{LN_MULTIPLY_PROPORTION}"],
        )
        
        
        arcpy.management.CalculateField(
            OUTPUT_LAYER,
            SHANON_INDEX,
            expression =   "!{}! * (-1)".format(f"SUM_{LN_MULTIPLY_PROPORTION}") ,
            expression_type = "PYTHON"        
        )


        arcpy.management.DeleteField(
            OUTPUT_LAYER, 
            [SHAPE_AREA_FIELD_NAME, LANDSCAPE_PROPORTION_FIELD_NAME, LN_MULTIPLY_PROPORTION, 
             f"SUM_{LN_MULTIPLY_PROPORTION}", f"SUM_{SHAPE_AREA_FIELD_NAME}"], 
            "DELETE_FIELDS")


        return

class Tool3(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mean Shape Complexity"
        self.description = "Based on https://info.undp.org/docs/pdc/Documents/ECU/MetricasFragstats-English.pdf (C19)"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
       
       
        # Input SHP layer with intersected grid and land patches
        param0 = arcpy.Parameter(
            displayName="Input Layer",
            name="input_layer",
            datatype=["DEShapeFile","DEFeatureClass"],
            parameterType="Required",
            direction="Input"
            )
        
     
        # Output SHP Layer
        param1 = arcpy.Parameter(
            displayName="Output Layer",
            name="out_features",
            datatype=["DEShapeFile","DEFeatureClass"],
            parameterType="Required",
            direction="Output")
        
    
        # Name of field with information about gridcell ID
        param2 = arcpy.Parameter(
            displayName="GRID ID",
            name="GRID_id",
            datatype="Field",
            parameterType="Required",
            direction="Input")
    
        
        """# Name of output Table 
        param4 = arcpy.Parameter(
            displayName="Path to Output Table Name",
            name="output_table",
            datatype="GPString",
            parameterType="Required",
            direction="Output")"""
        
        
        param1.parameterDependencies = [param0.name]
        param1.schema.clone = True
        
        #params = [param0, param1, param2, param3]
        params = [param0, param1, param2]
        
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
        #OUTPUT_TABLE = parameters[3].valueAsText
        
        AREA_FIELD_NAME = "AREA"
        PERIM_FIELD_NAME = "PERIM" # perimeter of gridcell (all sides)
        
        MSC_FIELD_NAME = "MSC"

        # Define datype for fields
        FIELD_TYPE = "DOUBLE"
      
        # Create copy of input layer to preserve original data
        input_layer_copy = arcpy.Copy_management(
            INPUT_LAYER, 
            OUTPUT_LAYER)
        

        # calculate area of land patch and write it into column
        arcpy.management.CalculateGeometryAttributes(
            input_layer_copy,
            [[PERIM_FIELD_NAME, "PERIMETER_LENGTH"], [AREA_FIELD_NAME, "AREA"]],
            length_unit = "METERS",
            area_unit = "SQUARE_METERS"

        )

        # sum patch perimeter and patch area for each ID
        # write it into separate table and store in memory
        summary_table = r"in_memory/summary_table"
        arcpy.analysis.Statistics(
            OUTPUT_LAYER,
            summary_table,
            [[PERIM_FIELD_NAME, "SUM"], [AREA_FIELD_NAME, "SUM"]],
            GRID_ID
        )

        arcpy.management.CalculateField(
            OUTPUT_LAYER,
            MSC_FIELD_NAME,
            expression="round((!{}! / (2 * math.sqrt(math.pi * !{}!))), 3)".format(PERIM_FIELD_NAME, AREA_FIELD_NAME),
            expression_type="PYTHON"
        )
       
        arcpy.management.DeleteField(
            OUTPUT_LAYER, 
            [PERIM_FIELD_NAME, AREA_FIELD_NAME], 
            "DELETE_FIELDS")
        
        """arcpy.management.MakeTableView(
            OUTPUT_LAYER, 
            OUTPUT_TABLE, 
            )
        
        arcpy.conversion.TableToExcel(
            OUTPUT_TABLE,
            f"{OUTPUT_TABLE}.xlsx
        )"""

        return
