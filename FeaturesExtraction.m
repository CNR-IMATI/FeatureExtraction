function [curveId, curvePars, axis1, axis2, indPts] = FeaturesExtraction(modelName, nameFct, curve, D, parN, minPts, eps, facet,  erosion)
    % FeaturesExtraction: 
    % Performs the whole extraction pipeline using the curve selected among the atlas of available curves
    % 
    % Input:    
    %  'modelName'  is the name of the input model (it should be a ply file)
    %  'nameFct'    is the name of the facet (it should be a fct file)
    %  'curve'      is the selected family of curves: 
    %               'CC' is the D-convexities curve 
    %               'C' is the citrus curve
    %               'CF' is the cirsumference
    %               'E' is the ellipse 
    %               'GP' is the geometric petal curve 
    %               'S' is the spiral curve 
    %  'D'          is the number of convexities of the D-convexities curve (default value is 8) 
    %               or the degree of 'eccentricity' of the geometric petal
    %               curve (default value is 40); for the other curves this value is not taken into account
    %  'parN'       is the positive integer number (less or equal to 60) used to estimate the threshold 
    %               for DBSCAN (default=50)
    %  'minPts'     minimum number of points required from DBSCAN to form a dense region (default=10)
    %  'eps'        is the maximum admissible distance of the selected points from the recognized curve 
    %               (default=0.85)
    %  'facet'      is a boolean indicating whether to use the faceting results or not. 
    %               '1' use facets or '0' do not use facets (default=0)
    %  'erosion'    is the number of layers to be deleted from the border of the facet (default=3)
    % 
    % Output: 
    %   'curveId'   the chosen family of curves
    %   'curvePars' matrix containing the values of the parameters of the recognized curve
    %               - first parameter
    %               - second parameter
    %   'axis1'     is a cell array containing the horizontal axis of the feature (struct with two points)
    %   'axis2'     is a cell array containing the vertical axis of the feature (struct with two points)
    %   'indPts'    is a cell array; each entry contains the list of the recognized feature points 
    %                   (that is, points close to the recognized curve of
    %                   parameters a and c)
    % 
    % 
    
if nargin < 3
        disp(sprintf('Error: Model, facet and family of curves are required.'));
    end
    if nargin ==3
        if curve == 'CC' 
            D=8; 
        else
            D=40;
        end
        parN=50;
        minPts=10;
        Eps =0.85;
        facet=0; 
        erosion=3;
    end
    if nargin ==4
        parN=50;
        minPts=10;
        Eps =0.85;
        facet=0; 
        erosion=3;
    end
    if nargin ==5
        minPts=10;
        Eps = 0.85;
        facet=0; 
        erosion=3;
    end
    if nargin ==6
        Eps=0.85;
        facet=0; 
        erosion=3;
    end
    if nargin ==7
        facet=0; 
        erosion=3;
    end
    if nargin ==8
        erosion=3;
    end

    if parN>60 
        disp(sprintf('Error: The parameter parN is out of range.'));
    end
   
    % set the paths
    path(path, 'toolbox_graph/');
    path(path, 'toolbox_graph/toolbox/');
    clear options;
    
    switch curve
        case 'C'
            disp('Use the Citrus curve')
            [curveId, curvePars, axis1, axis2, indPts] = extract_citrus(modelName, nameFct, parN, minPts, eps, facet,  erosion);
       case 'CC'
            disp('Use the convex curve')
            [curveId, curvePars, axis1, axis2, indPts] = extract_convexCurve(modelName, nameFct, D, parN, minPts, eps, facet, erosion);
        case 'CF'
            disp('Use the circumference')
            [curveId, curvePars, axis1, axis2, indPts] = extract_circumf(modelName, nameFct, parN, minPts, eps, facet, erosion);
        case 'E'
            disp('Use the ellipse curve')
            [curveId, curvePars, axis1, axis2, indPts] = extract_ellipse(modelName, nameFct, parN, minPts, eps, facet, erosion);
        case 'GP' 
            disp('Use the Geometric Petal curve')
            [curveId, curvePars, axis1, axis2, indPts] = extract_geometricPetalPolar_2(modelName, nameFct, D, parN, minPts, eps, facet, erosion);
        case 'S' 
            disp('Use the spiral curve')
            [curveId, curvePars, axis1, axis2, indPts] = extract_spiral(modelName, nameFct, parN, minPts, eps, facet, erosion);    
        otherwise
            disp('The given family of curves is not included in the atlas')
    end
        
    
end
    