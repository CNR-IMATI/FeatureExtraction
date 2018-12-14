function OUT = main(modelName, nameFct, curve, D, parN, minPts, eps, facet, erosion, fileout)

    %call the function for feature extraction FeaturesExtraction
    [curveId, curvePars, axis1, axis2, indPts] = FeaturesExtraction(modelName, nameFct, curve, D, parN, minPts, eps, facet,  erosion);
    
    s = size(curvePars,1);
    myCell = cell(1,s);
    for i=1:s
        res.curveId = curveId; 
        res.param = curvePars(i,:);
        res.axis1 = axis1{1,i};
        res.axis2 = axis2{1,i};
        VV = indPts{1,i};
        VV(1,:) = VV(1,:)-1;
        res.vertices = VV;
        myCell{1,i} = res;
    end
    myOut.features = myCell;
    OUT = jsonencode(myOut);
    fileID = fopen(fileout,'w');
    fprintf(fileID, '%s\n', OUT);
    fclose(fileID);
    %save(fileout,'-ascii','OUT');
end
        
