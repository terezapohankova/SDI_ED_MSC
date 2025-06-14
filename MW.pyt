# -*- coding: utf-8 -*-

import arcpy
import math
import sys
import statistics
import os


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Edge Density & Shanon Diversity Index & Mean Complixity Index"
        self.alias = "ED_SDI_MSC"

        
        # List of tool classes associated with this toolbox
        self.tools = [ED, SDI, MSC, MSC_MW, ED_MW, SDI_MW]

def get_parameter_info(input_datatype, output_datatype, column_datatype):
    """Shared function for getting parameter definitions with customizable data types."""
    # Define the first parameter (input layer or table)
    param0 = arcpy.Parameter(
        displayName="Input Layer",  # Label for the input parameter
        name="input_layer",  # Name of the input parameter
        datatype=input_datatype,  # Data type for input (customizable for each class)
        parameterType="Required",  # This parameter is required
        direction="Input"  # Direction of data flow (input)
    )
    
    # Define the second parameter (output layer or table)
    param1 = arcpy.Parameter(
        displayName="Output Layer (FeatureClass or Shapefile)",  # Label for the output parameter
        name="out_features",  # Name of the output parameter
        datatype=output_datatype,  # Data type for output (customizable for each class)
        parameterType="Required",  # This parameter is required
        direction="Output"  # Direction of data flow (output)
    )

    # Define the third parameter (column to be normalized)
    param2 = arcpy.Parameter(
        displayName="Column ID",  # Label for the column parameter
        name="column",  # Name of the column parameter
        datatype=column_datatype,  # Data type for the column (customizable for each class)
        parameterType="Required",  # This parameter is required
        direction="Input"  # Direction of data flow (input)
    )

    # Number of sides of grid (square -> 4, octagon ->8")
    param3 = arcpy.Parameter(
        displayName="Number of Grid Sides (4, 8 etc.)",
        name="grid_sides",
        datatype="GPLong",
        parameterType="Required",
        direction="Input")

    # Set parameter dependencies (output depends on input)
    param1.parameterDependencies = [param0.name]
    param1.schema.clone = True  # Cloning schema to ensure output matches input format

    # Combine the parameters into a list
    params = [param0, param1, param2, param3]
    
    # Set a filter to accept only valid field types for the 'column' parameter
    params[2].filter.list = ['Short', 'Long', 'Float', 'Double']  # Field data types
    params[2].parameterDependencies = [params[0].name]  # Column depends on input layer

    return params  # Return the list of parameters

class ED(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Density"
        self.description = "Calculation of Edge Density based on https://fragstats.org/index.php/fragstats-metrics/patch-based-metrics/area-and-edge-metrics/l4-edge-density"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
       
       # Get the parameter definitions with custom data types for this class
        return get_parameter_info(
            input_datatype=["DEShapeFile","GPFeatureLayer"],  # Input is a Feature Layer
            output_datatype=["GPFeatureLayer", "DEShapeFile"],  # Output is a Feature Layer
            column_datatype="Field"  # Column is a field
        )     
        

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
            
        # Define new field names
        LAND_LENGTH_FIELD_NAME = "L_LAND"  # perimeter of landpatch (all sides)
        LAND_AREA_FIELD_NAME = "A_LAND"
        
        GRID_AREA_FIELD_NAME = "A_GRID"
        GRID_LENGTH_FIELD_NAME = "L_GRID" # perimeter of gridcell (all sides)
        
        #names for tool polygon neighbours to detect shared borders
        SHARED_BORDER_FIELD_NAME = "LENGTH"                
        ED_FIELD_NAME = "EDGE_DENS"
        
        # Define datype for fields
        FIELD_TYPE = "DOUBLE"
      
        # Create copy of input layer to preserve original data
        input_layer_copy = arcpy.Copy_management(
            INPUT_LAYER, 
            OUTPUT_LAYER)
    
        
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

        with arcpy.da.UpdateCursor(OUTPUT_LAYER, SHARED_BORDER_FIELD_NAME) as cursor:
            for row in cursor:
                if (row[0] == None):
                    row[0] = 0
                cursor.updateRow(row)

        

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
           expression_type="PYTHON",
           field_type=FIELD_TYPE           
        )
        
        # calclaute perimeter of gridcell
        # calculate the square root of SUM_{LAND_AREA_FIELD_NAME} -> get the length of one side of the grid
        # multiply by the number of sides of the grid -> get the perimeter of the grid
        arcpy.management.CalculateField(
            OUTPUT_LAYER,
            GRID_LENGTH_FIELD_NAME,
            expression =   "math.sqrt(!{}!) * {}".format(f"SUM_{LAND_AREA_FIELD_NAME}", GRID_SIDES),
            expression_type = "PYTHON",
            field_type=FIELD_TYPE       
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
           expression_type="PYTHON",
           field_type=FIELD_TYPE           
        )

        arcpy.management.DeleteField(
            OUTPUT_LAYER, 
            [LAND_LENGTH_FIELD_NAME, LAND_AREA_FIELD_NAME, GRID_AREA_FIELD_NAME, 
             GRID_LENGTH_FIELD_NAME, f"SUM_{LAND_AREA_FIELD_NAME}", f"SUM_{LAND_LENGTH_FIELD_NAME}", SHARED_BORDER_FIELD_NAME], 
            "DELETE_FIELDS")  
        
        return
    
class SDI(object): 
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shanon Diversity Index"
        self.description = "Calculation of Shanon Diversity Index on https://fragstats.org/index.php/fragstats-metrics/patch-based-metrics/diversity-metrics/l4-shannons-diversity-index"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
       
       
        # Get the parameter definitions with custom data types for this class
        return get_parameter_info(
            input_datatype=["DEShapeFile","GPFeatureLayer"],  # Input is a Feature Layer
            output_datatype=["GPFeatureLayer", "DEShapeFile"],  # Output is a Feature Layer
            column_datatype="Field"  # Column is a field
        )    
        
        
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
        FIELD_TYPE = "DOUBLE"
    
        # Create copy of input layer to preserve original data
        input_layer_copy = arcpy.Copy_management(
            INPUT_LAYER, 
            OUTPUT_LAYER)
        
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
            expression_type = "PYTHON",
            field_type=FIELD_TYPE        
        )
        
        arcpy.management.CalculateField(
            OUTPUT_LAYER,
            LN_MULTIPLY_PROPORTION,
            expression =   "!{}! * math.log(!{}!)".format(LANDSCAPE_PROPORTION_FIELD_NAME, LANDSCAPE_PROPORTION_FIELD_NAME),
            expression_type = "PYTHON",
            field_type=FIELD_TYPE        
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
            expression_type = "PYTHON",
            field_type=FIELD_TYPE         
        )

        arcpy.management.DeleteField(
            OUTPUT_LAYER, 
            [SHAPE_AREA_FIELD_NAME, LANDSCAPE_PROPORTION_FIELD_NAME, LN_MULTIPLY_PROPORTION, 
             f"SUM_{LN_MULTIPLY_PROPORTION}", f"SUM_{SHAPE_AREA_FIELD_NAME}"], 
            "DELETE_FIELDS")

        return

class MSC(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mean Shape Complexity"
        self.description = "Based on https://info.undp.org/docs/pdc/Documents/ECU/MetricasFragstats-English.pdf (C19)"
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
       
       
        # Get the parameter definitions with custom data types for this class
        return get_parameter_info(
            input_datatype=["DEShapeFile","GPFeatureLayer"],  # Input is a Feature Layer
            output_datatype=["GPFeatureLayer", "DEShapeFile"],  # Output is a Feature Layer
            column_datatype="Field"  # Column is a field
        )    
        

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
        
        arcpy.management.AddField(OUTPUT_LAYER, 
                                  MSC_FIELD_NAME, 
                                  "DOUBLE")
        

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
    

    
class MSC_MW:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Mean Shape Complexity (Moving Window)"
        self.description = "Calculate mean Shape Complexity for Defined IDs \
            using 500 m focal window"

    def getParameterInfo(self):
        """Define the parameters for the tool."""
        
        # Define the first parameter (input layer or table)
        param0 = arcpy.Parameter(
            displayName="Input Layer",  # Label for the input parameter
            name="input_layer",  # Name of the input parameter
            datatype=["GPFeatureLayer", "DEShapeFile"],  # Data type for input (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Define the second parameter (output layer)
        param1 = arcpy.Parameter(
            displayName="Output Layer",  # Label for the input parameter
            name="output_layer",  # Name of the input parameter
            datatype=["GPFeatureLayer", "DEShapeFile"],  # Data type for input (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Output"  # Direction of data flow (output)
        )

        # Define the third parameter (column to be normalized)
        # ID for each will the MSC be calculated. Around this ID will be crated 500 m focal window
        param2 = arcpy.Parameter(
            displayName="Column with ID",  # Label for the column parameter
            name="column_small",  # Name of the column parameter
            datatype="Field",  # Data type for the column (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Combine the parameters into a list
        params = [param0, param1, param2]

        # Set parameter dependencies (output depends on input)
        param1.parameterDependencies = [param0.name]
        param1.schema.clone = True  # Cloning schema to ensure output matches input format

        # Set a filter to accept only valid field types for the 'column' parameter
        params[2].filter.list = ['Short', 'Long', 'Float', 'Double']  # Field data types
        params[2].parameterDependencies = [params[0].name]  # Column depends on input layer


        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""  
        return

    def execute(self, parameters, messages):
        input_layer = parameters[0].valueAsText
        output_layer = parameters[1].valueAsText
        small_col = parameters[2].valueAsText

        arcpy.management.RepairGeometry(input_layer, 'DELETE_NULL', 'ESRI')

        unique_ids = set()
        with arcpy.da.SearchCursor(input_layer, [small_col]) as cursor:
            for row in cursor:
                unique_ids.add(row[0])

        merged_outputs = []

        for current_id in unique_ids:
            arcpy.AddMessage(f"\nProcessing ID: {current_id}")

            # Step 1: select current feature
            arcpy.management.MakeFeatureLayer(input_layer, "input_layer_lyr")
            clause = f"{small_col} = '{current_id}'" if isinstance(current_id, str) else f"{small_col} = {current_id}"
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", clause)

            # Step 2: first-order neighbors
            arcpy.management.MakeFeatureLayer(input_layer, "first_neighbors_lyr")
            arcpy.management.SelectLayerByLocation("first_neighbors_lyr", "INTERSECT", "input_layer_lyr", selection_type="NEW_SELECTION")
            first_neighbor_ids = {row[0] for row in arcpy.da.SearchCursor("first_neighbors_lyr", [small_col]) if row[0] != current_id}

            if not first_neighbor_ids:
                arcpy.AddMessage(f"No neighbors found for ID {current_id}")
                continue

            # Step 3: second-order neighbors
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", " OR ".join(
                [f"{small_col} = '{nid}'" if isinstance(nid, str) else f"{small_col} = {nid}" for nid in first_neighbor_ids]
            ))
            arcpy.management.MakeFeatureLayer(input_layer, "second_neighbors_lyr")
            arcpy.management.SelectLayerByLocation("second_neighbors_lyr", "INTERSECT", "input_layer_lyr", selection_type="NEW_SELECTION")

            second_neighbor_ids = {
                row[0] for row in arcpy.da.SearchCursor("second_neighbors_lyr", [small_col])
                if row[0] != current_id and row[0] not in first_neighbor_ids
            }

            all_ids = first_neighbor_ids.union(second_neighbor_ids)
            all_ids.add(current_id)

            where_clause = " OR ".join(
                [f"{small_col} = '{nid}'" if isinstance(nid, str) else f"{small_col} = {nid}" for nid in all_ids]
            )
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", where_clause)

            # Step 4: Calculate MSC for each selected feature
            perims = []
            areas = []

            with arcpy.da.SearchCursor("input_layer_lyr", ["SHAPE@"]) as cursor:
                for row in cursor:
                    geom = row[0]
                    perim = geom.length  # perimeter in map units
                    area = geom.area     # area in map units^2
                    perims.append(perim)
                    areas.append(area)

            # Calculate MSC values for each feature (MSC = perimeter / (2 * sqrt(pi * area)))
            msc_values = []
            for p, a in zip(perims, areas):
                if a > 0:
                    msc = p / (2 * math.sqrt(math.pi * a))
                    msc_values.append(msc)

            # Calculate average MSC over selected features
            if msc_values:
                avg_msc = sum(msc_values) / len(msc_values)
            else:
                avg_msc = None

            if avg_msc is None:
                arcpy.AddMessage(f"No MSC values to calculate average for ID {current_id}")
                continue

            # Step 5: Write avg MSC only to rows with current_id
            field_name = "MSC_MW"  # new field for average MSC

            existing_fields = [f.name for f in arcpy.ListFields("input_layer_lyr")]
            if field_name not in existing_fields:
                arcpy.management.AddField("input_layer_lyr", field_name, "DOUBLE")

            with arcpy.da.UpdateCursor("input_layer_lyr", [small_col, field_name]) as cursor:
                for row in cursor:
                    if row[0] == current_id:
                        row[1] = avg_msc
                        cursor.updateRow(row)

            # Select only features with current_id in the layer
            clause_current = f"{small_col} = '{current_id}'" if isinstance(current_id, str) else f"{small_col} = {current_id}"
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", clause_current)

            # Save only the features with current_id
            out_fc = os.path.join("in_memory", f"nei2_{current_id}")
            if arcpy.Exists(out_fc):
                arcpy.Delete_management(out_fc)

            arcpy.management.CopyFeatures("input_layer_lyr", out_fc)
            merged_outputs.append(out_fc)

        # Step 7: merge all outputs
        if merged_outputs:
            arcpy.management.Merge(merged_outputs, output_layer)
            arcpy.AddMessage(f"Merged {len(merged_outputs)} layers into {output_layer}")
        else:
            arcpy.AddMessage("No valid results to merge.")


class ED_MW:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Edge Density (Moving Window)"
        self.description = "Calculate mean Shape Complexity for Defined IDs \
            using 500 m focal window"

    def getParameterInfo(self):
        """Define the parameters for the tool."""
        
        # Define the first parameter (input layer or table)
        param0 = arcpy.Parameter(
            displayName="Input Layer",  # Label for the input parameter
            name="input_layer",  # Name of the input parameter
            datatype=["GPFeatureLayer", "DEShapeFile"],  # Data type for input (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Define the second parameter (output layer)
        param1 = arcpy.Parameter(
            displayName="Output Layer",  # Label for the input parameter
            name="output_layer",  # Name of the input parameter
            datatype=["GPFeatureLayer", "DEShapeFile"],  # Data type for input (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Output"  # Direction of data flow (output)
        )

        # Define the third parameter (column to be normalized)
        # ID for each will the MSC be calculated. Around this ID will be crated 500 m focal window
        param2 = arcpy.Parameter(
            displayName="Column with ID",  # Label for the column parameter
            name="column_small",  # Name of the column parameter
            datatype="Field",  # Data type for the column (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Define the third parameter (column to be normalized)
        param3 = arcpy.Parameter(
            displayName="Number of Grid Sides (4 -> square)",  # Label for the column parameter
            name="grid_sides",  # Name of the column parameter
            datatype="GPLong",  # Data type for the column (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Combine the parameters into a list
        params = [param0, param1, param2, param3]

        # Set parameter dependencies (output depends on input)
        param1.parameterDependencies = [param0.name]
        param1.schema.clone = True  # Cloning schema to ensure output matches input format

        # Set a filter to accept only valid field types for the 'column' parameter
        params[2].filter.list = ['Short', 'Long', 'Float', 'Double']  # Field data types
        params[2].parameterDependencies = [params[0].name]  # Column depends on input layer


        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""  
        return

    def execute(self, parameters, messages):


        input_layer = parameters[0].valueAsText
        output_layer = parameters[1].valueAsText
        small_col = parameters[2].valueAsText
        grid_sides = parameters[3].value  # you might want to use this if relevant later

        ED_MW_NAME = "ED_MW"

        arcpy.management.RepairGeometry(input_layer, 'DELETE_NULL', 'ESRI')

        unique_ids = set()
        with arcpy.da.SearchCursor(input_layer, [small_col]) as cursor:
            for row in cursor:
                unique_ids.add(row[0])

        merged_outputs = []

        for current_id in unique_ids:
            arcpy.AddMessage(f"\nProcessing ID: {current_id}")

            # Select focal feature
            arcpy.management.MakeFeatureLayer(input_layer, "input_layer_lyr")
            clause = f"{small_col} = '{current_id}'" if isinstance(current_id, str) else f"{small_col} = {current_id}"
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", clause)

            # First-order neighbors
            arcpy.management.MakeFeatureLayer(input_layer, "first_neighbors_lyr")
            arcpy.management.SelectLayerByLocation("first_neighbors_lyr", "INTERSECT", "input_layer_lyr", selection_type="NEW_SELECTION")
            first_neighbor_ids = {row[0] for row in arcpy.da.SearchCursor("first_neighbors_lyr", [small_col]) if row[0] != current_id}

            if not first_neighbor_ids:
                arcpy.AddMessage(f"No neighbors found for ID {current_id}")
                continue

            # Second-order neighbors
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", " OR ".join(
                [f"{small_col} = '{nid}'" if isinstance(nid, str) else f"{small_col} = {nid}" for nid in first_neighbor_ids]
            ))
            arcpy.management.MakeFeatureLayer(input_layer, "second_neighbors_lyr")
            arcpy.management.SelectLayerByLocation("second_neighbors_lyr", "INTERSECT", "input_layer_lyr", selection_type="NEW_SELECTION")
            second_neighbor_ids = {
                row[0] for row in arcpy.da.SearchCursor("second_neighbors_lyr", [small_col])
                if row[0] != current_id and row[0] not in first_neighbor_ids
            }

            all_ids = first_neighbor_ids.union(second_neighbor_ids)
            all_ids.add(current_id)

            where_clause = " OR ".join(
                [f"{small_col} = '{nid}'" if isinstance(nid, str) else f"{small_col} = {nid}" for nid in all_ids]
            )
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", where_clause)

            # Copy selected features into memory
            subset_fc = os.path.join("in_memory", f"subset_{current_id}")
            arcpy.management.CopyFeatures("input_layer_lyr", subset_fc)

            # Calculate shared borders
            polygon_neighb = os.path.join("in_memory", f"neighb_{current_id}")
            arcpy.analysis.PolygonNeighbors(
                subset_fc,
                polygon_neighb,
                small_col,
                area_overlap="NO_AREA_OVERLAP",
                both_sides="NO_BOTH_SIDES",
                out_area_units="SQUARE_METERS"
            )

            arcpy.management.AlterField(polygon_neighb, f"src_{small_col}", small_col, small_col)

            arcpy.management.JoinField(
                subset_fc,
                small_col,
                polygon_neighb,
                small_col,
                ["LENGTH"]
            )

            with arcpy.da.UpdateCursor(subset_fc, ["LENGTH"]) as cursor:
                for row in cursor:
                    if row[0] is None:
                        row[0] = 0
                    cursor.updateRow(row)

            # Geometry attributes
            arcpy.management.CalculateGeometryAttributes(
                subset_fc,
                [["L_LAND", "PERIMETER_LENGTH"], ["A_LAND", "AREA"]],
                length_unit="METERS",
                area_unit="SQUARE_METERS"
            )

            # Summary stats
            summary_table = os.path.join("in_memory", f"summary_{current_id}")
            arcpy.analysis.Statistics(
                subset_fc,
                summary_table,
                [["L_LAND", "SUM"], ["A_LAND", "SUM"]],
                small_col
            )

            arcpy.management.JoinField(
                subset_fc,
                small_col,
                summary_table,
                small_col,
                ["SUM_L_LAND", "SUM_A_LAND"]
            )

            # Calculate grid area and perimeter using exponentiation for sqrt
            arcpy.management.CalculateField(subset_fc, "A_GRID", "!SUM_A_LAND!", "PYTHON3")
            arcpy.management.CalculateField(subset_fc, "L_GRID", "!SUM_A_LAND! ** 0.5 * 4", "PYTHON3")

            expression = """(
            float(str(!SUM_L_LAND!).replace(',', '.')) -
            float(str(!L_GRID!).replace(',', '.')) -
            float(str(!LENGTH!).replace(',', '.'))
            ) / float(str(!A_GRID!).replace(',', '.'))"""

            arcpy.management.CalculateField(
                subset_fc,
                ED_MW_NAME,
                expression,
                expression_type="PYTHON3")

            # Get value for focal feature
            edge_density = None
            with arcpy.da.SearchCursor(subset_fc, [small_col, ED_MW_NAME]) as cursor:
                for row in cursor:
                    if row[0] == current_id:
                        edge_density = row[1]
                        break

            if edge_density is not None:
                existing_fields = [f.name for f in arcpy.ListFields("input_layer_lyr")]
                if ED_MW_NAME not in existing_fields:
                    arcpy.management.AddField("input_layer_lyr", ED_MW_NAME, "DOUBLE")

                with arcpy.da.UpdateCursor("input_layer_lyr", [small_col, ED_MW_NAME]) as cursor:
                    for row in cursor:
                        if row[0] == current_id:
                            row[1] = edge_density
                            cursor.updateRow(row)
                # no break here, updates all matching rows


                # Select focal feature only
                clause_current = f"{small_col} = '{current_id}'" if isinstance(current_id, str) else f"{small_col} = {current_id}"
                arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", clause_current)

                out_fc = os.path.join("in_memory", f"ed_{current_id}")
                if arcpy.Exists(out_fc):
                    arcpy.Delete_management(out_fc)

                arcpy.management.CopyFeatures("input_layer_lyr", out_fc)
                merged_outputs.append(out_fc)

        # Merge all output features
        if merged_outputs:
            arcpy.management.Merge(merged_outputs, output_layer)
            arcpy.AddMessage(f"Merged {len(merged_outputs)} layers into {output_layer}")
        else:
            arcpy.AddMessage("No valid results to merge.")


class SDI_MW:
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Shanon Diversity Index (Moving Window)"
        self.description = "Calculate mean the Index for Defined IDs \
            using 500 m focal window"

    def getParameterInfo(self):
        """Define the parameters for the tool."""
        
        # Define the first parameter (input layer or table)
        param0 = arcpy.Parameter(
            displayName="Input Layer",  # Label for the input parameter
            name="input_layer",  # Name of the input parameter
            datatype=["GPFeatureLayer", "DEShapeFile"],  # Data type for input (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Define the second parameter (output layer)
        param1 = arcpy.Parameter(
            displayName="Output Layer",  # Label for the input parameter
            name="output_layer",  # Name of the input parameter
            datatype=["GPFeatureLayer", "DEShapeFile"],  # Data type for input (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Output"  # Direction of data flow (output)
        )

        # Define the third parameter (column to be normalized)
        # ID for each will the MSC be calculated. Around this ID will be crated 500 m focal window
        param2 = arcpy.Parameter(
            displayName="Column with ID",  # Label for the column parameter
            name="column_small",  # Name of the column parameter
            datatype="Field",  # Data type for the column (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Define the third parameter (column to be normalized)
        param3 = arcpy.Parameter(
            displayName="Number of Grid Sides (4 -> square)",  # Label for the column parameter
            name="grid_sides",  # Name of the column parameter
            datatype="GPLong",  # Data type for the column (customizable for each class)
            parameterType="Required",  # This parameter is required
            direction="Input"  # Direction of data flow (input)
        )

        # Combine the parameters into a list
        params = [param0, param1, param2, param3]

        # Set parameter dependencies (output depends on input)
        param1.parameterDependencies = [param0.name]
        param1.schema.clone = True  # Cloning schema to ensure output matches input format

        # Set a filter to accept only valid field types for the 'column' parameter
        params[2].filter.list = ['Short', 'Long', 'Float', 'Double']  # Field data types
        params[2].parameterDependencies = [params[0].name]  # Column depends on input layer


        return params

    def isLicensed(self):
        """Set whether the tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter. This method is called after internal validation."""  
        return

    def execute(self, parameters, messages):


        input_layer = parameters[0].valueAsText
        output_layer = parameters[1].valueAsText
        small_col = parameters[2].valueAsText
        grid_sides = parameters[3].value  # Currently unused

        SDI_MW_NAME = "SDI_MW"

        arcpy.management.RepairGeometry(input_layer, 'DELETE_NULL', 'ESRI')

        unique_ids = set()
        with arcpy.da.SearchCursor(input_layer, [small_col]) as cursor:
            for row in cursor:
                unique_ids.add(row[0])

        merged_outputs = []

        for current_id in unique_ids:
            arcpy.AddMessage(f"\nProcessing ID: {current_id}")

            # Select focal feature
            arcpy.management.MakeFeatureLayer(input_layer, "input_layer_lyr")
            clause = f"{small_col} = '{current_id}'" if isinstance(current_id, str) else f"{small_col} = {current_id}"
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", clause)

            # First-order neighbors
            arcpy.management.MakeFeatureLayer(input_layer, "first_neighbors_lyr")
            arcpy.management.SelectLayerByLocation("first_neighbors_lyr", "INTERSECT", "input_layer_lyr", selection_type="NEW_SELECTION")
            first_neighbor_ids = {row[0] for row in arcpy.da.SearchCursor("first_neighbors_lyr", [small_col]) if row[0] != current_id}

            if not first_neighbor_ids:
                arcpy.AddMessage(f"No neighbors found for ID {current_id}")
                continue

            # Second-order neighbors
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", " OR ".join(
                [f"{small_col} = '{nid}'" if isinstance(nid, str) else f"{small_col} = {nid}" for nid in first_neighbor_ids]
            ))
            arcpy.management.MakeFeatureLayer(input_layer, "second_neighbors_lyr")
            arcpy.management.SelectLayerByLocation("second_neighbors_lyr", "INTERSECT", "input_layer_lyr", selection_type="NEW_SELECTION")
            second_neighbor_ids = {
                row[0] for row in arcpy.da.SearchCursor("second_neighbors_lyr", [small_col])
                if row[0] != current_id and row[0] not in first_neighbor_ids
            }

            all_ids = first_neighbor_ids.union(second_neighbor_ids)
            all_ids.add(current_id)

            where_clause = " OR ".join(
                [f"{small_col} = '{nid}'" if isinstance(nid, str) else f"{small_col} = {nid}" for nid in all_ids]
            )
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", where_clause)

            # Copy selected features into memory
            subset_fc = os.path.join("in_memory", f"subset_{current_id}")
            arcpy.management.CopyFeatures("input_layer_lyr", subset_fc)

            SHAPE_AREA_FIELD_NAME = "AREA_M"
            LANDSCAPE_PROPORTION_FIELD_NAME = "pi"
            LN_MULTIPLY_PROPORTION = "pilnpi"
            SHANON_INDEX = "SDI"
            FIELD_TYPE = "DOUBLE"

            # Calculate area of each patch
            arcpy.management.CalculateGeometryAttributes(
                subset_fc,
                [[SHAPE_AREA_FIELD_NAME, "AREA"]],
                area_unit="SQUARE_METERS"
            )

            # Summarize total area per ID
            summary_table_sdi = os.path.join("in_memory", f"summary_area_{current_id}")
            arcpy.analysis.Statistics(
                subset_fc,
                summary_table_sdi,
                [[SHAPE_AREA_FIELD_NAME, "SUM"]],
                small_col
            )

            # Join SUM_AREA back to subset_fc
            arcpy.management.JoinField(
                subset_fc,
                small_col,
                summary_table_sdi,
                small_col,
                [f"SUM_{SHAPE_AREA_FIELD_NAME}"]
            )

            # Calculate pi = patch area / total area
            arcpy.management.CalculateField(
                subset_fc,
                LANDSCAPE_PROPORTION_FIELD_NAME,
                expression=f"!{SHAPE_AREA_FIELD_NAME}! / !SUM_{SHAPE_AREA_FIELD_NAME}!",
                expression_type="PYTHON3",
                field_type=FIELD_TYPE
            )

            # Calculate pi * ln(pi)
            arcpy.management.CalculateField(
                subset_fc,
                LN_MULTIPLY_PROPORTION,
                expression=f"!{LANDSCAPE_PROPORTION_FIELD_NAME}! * math.log(!{LANDSCAPE_PROPORTION_FIELD_NAME}!)",
                expression_type="PYTHON3",
                field_type=FIELD_TYPE
            )

            # Sum pilnpi per group
            summary_table_pilnpi = os.path.join("in_memory", f"summary_pilnpi_{current_id}")
            arcpy.analysis.Statistics(
                subset_fc,
                summary_table_pilnpi,
                [[LN_MULTIPLY_PROPORTION, "SUM"]],
                small_col
            )

            # Join pilnpi back
            arcpy.management.JoinField(
                subset_fc,
                small_col,
                summary_table_pilnpi,
                small_col,
                [f"SUM_{LN_MULTIPLY_PROPORTION}"]
            )

            # Calculate final SDI = -1 * SUM_pilnpi
            arcpy.management.CalculateField(
                subset_fc,
                SHANON_INDEX,
                expression=f"!SUM_{LN_MULTIPLY_PROPORTION}! * -1",
                expression_type="PYTHON3",
                field_type=FIELD_TYPE
            )

            # Get SDI value for current_id
            sdi_value = None
            with arcpy.da.SearchCursor(subset_fc, [small_col, SHANON_INDEX]) as cursor:
                for row in cursor:
                    if row[0] == current_id:
                        sdi_value = row[1]
                        break

            # Optionally clean up intermediate fields
            arcpy.management.DeleteField(
                subset_fc,
                [SHAPE_AREA_FIELD_NAME, LANDSCAPE_PROPORTION_FIELD_NAME, LN_MULTIPLY_PROPORTION,
                f"SUM_{SHAPE_AREA_FIELD_NAME}", f"SUM_{LN_MULTIPLY_PROPORTION}"],
                "DELETE_FIELDS"
            )

            # Select focal feature only again
            clause_current = f"{small_col} = '{current_id}'" if isinstance(current_id, str) else f"{small_col} = {current_id}"
            arcpy.management.SelectLayerByAttribute("input_layer_lyr", "NEW_SELECTION", clause_current)

            out_fc = os.path.join("in_memory", f"ed_{current_id}")
            if arcpy.Exists(out_fc):
                arcpy.Delete_management(out_fc)

            arcpy.management.CopyFeatures("input_layer_lyr", out_fc)

            # Add SDI_MW field if not exists
            field_names = [f.name for f in arcpy.ListFields(out_fc)]
            if SDI_MW_NAME not in field_names:
                arcpy.management.AddField(out_fc, SDI_MW_NAME, "DOUBLE")

            # Write SDI value to SDI_MW field
            if sdi_value is not None:
                with arcpy.da.UpdateCursor(out_fc, [SDI_MW_NAME]) as cursor:
                    for row in cursor:
                        row[0] = sdi_value
                        cursor.updateRow(row)

            merged_outputs.append(out_fc)

        # Merge all output features
        if merged_outputs:
            arcpy.management.Merge(merged_outputs, output_layer)
            arcpy.AddMessage(f"Merged {len(merged_outputs)} layers into {output_layer}")
        else:
            arcpy.AddMessage("No valid results to merge.")
