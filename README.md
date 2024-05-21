# Edge Density (ED) for ArcGIS Pro

This repository contains a Python script for calculating Edge Density based on the Fragstats metrics (https://fragstats.org/index.php/fragstats-metrics/patch-based-metrics/area-and-edge-metrics/l4-edge-density).


# Requirements
Originally created for **ArcGIS Pro v. 3.3.0.**

*	Python (tested on 3.9.9)
*	arcpy (tested on 2.9)

# Input Data

* Vector layers
	* Shapefile layer with (multi)polygons for which ED is detected
  * the (multi)polygons must be cut using a grid with a regular number of sides (e.g. square, octagon)
	* each (multi)polygon must be assigned a numeric ID
	* the (multi)polygons must be cut using a grid with a regular number of sides (e.g. square, octagon)
   
https://github.com/terezapohankova/edge_density/assets/60270092/82909cd1-7c7a-4c63-bf76-c980758776d2


# Output Data
* Output table is aggregated by Grid ID (each row corresponds to one ID number)

|Paramater              |DataType            |Note        |
|-----------------------|----------------|----------------|
|Output Table          	|String          | in .xlsx format

Structure of Output Table in .xlsx format is
|FID              |GRID ID            |Edge Density        |
|-----------------------|----------------|----------------|
