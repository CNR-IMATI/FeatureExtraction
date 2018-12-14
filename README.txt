/********************************************
 * README.txt
 * Copyright (c) CNR-IMATI Ge - 2018
 * Author: Maria-Laura Torrente
 * All rights reserved
********************************************/


MATLAB routine to extract and localise the features on a given model.
This is an implementation of the feature curve recognition method described in the paper:
M.-L. Torrente, S. Biasotti, and B. Falcidieno. (2018) Recognition of feature curves on 3D shapes using an algebraic approach to Hough transforms. Pattern recognition ISSN: 0031-3203 Pergamon Press. 73


Tested on macOS Sierra version 10.12.6 and Ubuntu 14.04.

To properly work, the software needs (in the same folder of the function main.m) the following external software:
toolbox graph (add(genpath('toolbox_graph')): freely available at https://it.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph
cocoa5 package: freely available at http://cocoa.dima.unige.it/download/install5-mac.shtml
DBSCAN, for instance the software: http://yarpiz.com/255/ypml110-dbscan-clustering
RGB2Lab, for instance the function: https://it.mathworks.com/matlabcentral/fileexchange/24009-rgb2lab

/////////////////////////////////////////////////////////////////////////////////////////

Usage: OUT = main(modelName, curve, parN, minPts, eps, facet, nameFct, erosion, fileout)

The input is contained in the file ‘modelName’ and it is supposed to be a 3D model in .ply format.

The parameters are:
- ‘curve’: used to select the family of curves: 
    'CC’ is the D-convexities curve
    'CF’ is the circumference
    'C' is the citrus curve               
    'E' is the ellipse 
    'GP' is the geometric petal curve

-   'D' is the number of convexities of the D-convexities curve (default value is 8) 
or the degree of 'eccentricity' of the geometric petal curve (default value is 40); 
for the other curves this value is not taken into account;

- ‘parN’: positive integer number (less or equal to 60) used to estimate the threshold for DBSCAN (default=50);

- ‘minPts’: minimum number of points required from DBSCAN to form a dense region (default=10);

- ‘eps’: maximum admissible distance of the selected points from the recognized curve (default=0.85);

- ‘facet’: boolean indicating whether to use the faceting results or not: ‘1’ use facets or '0' do not use facets (default=0); 

- ‘nameFct’: name of the facet (it should be a fct file); it must be included if facet=1;

- ‘erosion’: number of layers to be deleted from the border of the facet (default=3).


The output, contained in `fileout’, is a sorted list of structure arrays, each of which contains the following fields:

- ‘curveId’: the chosen family of curves 

- ‘curvePars’: matrix containing the values of the parameters of the recognized curve;

- ‘axis1’: cell array containing the horizontal axis of the feature (struct with two points);

- ‘axis2’: cell array containing the vertical axis of the feature (struct with two points);

- ‘vertices’: list of the recognized feature points, that is, points close to the recognized curve.    


Example:
out=main('D_292_100K.ply','D_292_100K.fct', 'GP', 40, 55, 30, 0.5, 1,  7, 'D_292_100K_Features.json');

The file D_292_100K.ply is the 3D model in .ply format.
The file D_292_100K.fct contains the model faceting, in .fct format.
The file D_292_100K_Features.txt contains the output obtained in the previous example.

