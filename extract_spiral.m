function [curveId, curvePars, axis1, axis2, indPts] = extract_spiral(modelName, nameFct, parN, minPts, eps, facet, erosion)

    % Perform the whole extraction pipeline using spiral curve in polar form
    % 
    % Input:    
    %  'modelName'  is the name of the input model (it should be a ply file)
    %  'nameFct'    is the name of the facet (it should be a fct file)
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
    %               - first entry is the parameter a in the geometric petal curve
    %               - second entry is the stretching coefficient c
    %   'indPts' is a cell array; each entry contains the list of the recognized feature points 
    %               (that is, points close to the recognized curve)
    % 
    
    % declare the chosen family of curves
    curveId = 'spiral';
    
    % read the mesh
    clear options;
    options.name = modelName; % useful for displaying
    [vertices, faces] = read_mesh(modelName);
    numVertices = size(vertices, 1);

    % compute the curvatures
    options.curvature_smoothing = 10;
    options.verb = 0;
    [Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertices, faces);

    %histogram of mean curvature
    hmean=histogram(Cmean, 'visible', 'off');
    set(gcf,'Visible','off');
    numBins = hmean.NumBins;
    i=0;
    Values = hmean.Values;
    V=0;
    perc = 80/100;
    while (i<=numBins & V < perc*numVertices)
     i = i+1;
     V = V+Values(i);
    end
    BL = hmean.BinLimits;
    valueCmean = BL(1)+i*hmean.BinWidth;

    % if facet=1, read the external facet
    if facet==0
        Indici = [1:numVertices]';
        IndiciBorder = [];
    else
        [n_facets,filename,facets]=read_faceting_file(nameFct);
        for i=1:n_facets
            if (facets{i}.type == 'external')
                Indici = facets{i}.vertices;
                FacetPts = vertices(Indici,:);
                IndiciBorder = facets{i}.border;
                BorderPts = vertices(IndiciBorder,:);
                [IDX,DDD] = knnsearch(FacetPts, BorderPts,'K',erosion,'Distance','euclidean');
                break
            end
        end
        ElimIndici = [];
        for i=1:size(DDD,1)
            ElimIndici = [ElimIndici; Indici(IDX(i,:))];
        end
        Indici = setdiff(Indici, ElimIndici);
    end
    
    % select feature points with high (over the theshold) mean curvature
    numIndici = size(Indici,1);
    select=[];
    for i=1:numIndici
        if  (Cmean(Indici(i)) > valueCmean)
            select=[select' Indici(i)]';
        end
    end
    punti = zeros(size(select,1),3);
    for i=1:size(select,1)
        punti(i,:)=vertices(select(i),:);
    end

    % Find the minimum distance among the selected points and epsilon
    [IDX,DDD] = knnsearch(punti,punti,'k', 60);
    epsilon = mean(DDD(:, parN));
    
    % Run DBSCAN Clustering Algorithm
     IDX=DBSCAN(punti, epsilon, minPts);
    
%     %Plot Results
%     figure
%     PlotClusterinResult(punti, IDX);
%     title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', minPts = ' num2str(minPts) ')']);
    numCluster=max(IDX);
    Colors=hsv(numCluster);
    Style = 'x';
    MarkerSize = 8;
    Legends = {};

    UsedClusters = [];
    MM =cell(1,numCluster);
    RR =cell(1,numCluster);
    %MMx = cell(1,numCluster);
    %MMy = cell(1,numCluster);
    newRR =cell(1,numCluster);
    %newRR2 =cell(1,numCluster);
    newRR3 =cell(1,numCluster);
    SC = cell(1,numCluster);
    clusterCardinality = cell(1,numCluster);
    RelMaxValue = cell(1, numCluster);
    MaxCoords = cell(1, numCluster);
    
    % iteration on each cluster 
    for i=1:numCluster
        Xi=punti(IDX==i,:);
        Color = Colors(i,:);
        Legends{end+1} = ['Cluster #' num2str(i)];
        xyz = Xi;
        Ni = size(xyz,1);
        for ii=2:Ni
            dd = zeros(1,ii-1);
            for jj=1:(ii-1)
                dd(jj) = norm(Xi(ii,:)-Xi(jj,:));
            end
            if (max(dd)/2 > 2*epsilon)
                break;
            end
        end
        % we do not consider the very small clusters (whose radius < 2*epsilon)
        clusterCardinality{1,i}=Ni;
        if (max(dd)/2 < 2*epsilon)
            continue;
        else  
            % we sample the big clusters (the ones containing too many points)
            sogliaCluster = 500;
            if (Ni > sogliaCluster)
                perc = sogliaCluster./Ni;
                ptCloud = pointCloud(xyz);
                ptCloud = pcdownsample(ptCloud,'random',perc);
                xyz = ptCloud.Location;
                clusterCardinality{1,i}=sogliaCluster;
            end    
            s=size(UsedClusters,2);
            UsedClusters(1,s+1) =i;
            % make the mean (0,0,0)
            MM{1,i} = mean(xyz); 
            xyz = bsxfun(@minus, xyz, MM{1,i}); 
            %compute the plane of linear regression
            x=xyz(:,1); y=xyz(:,2); z=xyz(:,3);
            X=[ones(size(x)) x y];
            b=regress(z,X);
            
            n=[b(2) b(3) -1];
            v1=[1 0 b(2)];
            v2=cross(v1,n);
            R=[v1/norm(v1); v2/norm(v2); n/norm(n)]';
            RR{1,i} = R;
            xyznew=xyz*R;
            
            % compute minimal bounding box, and rotate
            xy2=[xyznew(:,1) xyznew(:,2)];
            mBB = minBoundingBox(xy2');
            
            if (norm(mBB(:,1)-mBB(:,2)) < norm(mBB(:,1)-mBB(:,4)) )
                m = (mBB(2,4)-mBB(2,1))/(mBB(1,4)-mBB(1,1));
                b = -m*mBB(1,4)+mBB(2,4);
                alpha = atan(m);
            else
                m = (mBB(2,2)-mBB(2,1))/(mBB(1,2)-mBB(1,1));
                b = -m*mBB(1,2)+mBB(2,2);
                alpha = atan(m);
            end
            R = [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
            newRR{1,i} = R;
            xy2new=xy2*R;
            mBB=R'*mBB;

            SP=zeros(4,2);
            [Mx,I]=sort(mBB(1,:));
            if mBB(2,I(1))<mBB(2,I(2))
                K=I(1);
            else
                K=I(2);
            end
            for T=1:4 
                I1=K+T-1;
                if I1>4
                    I1=I1-4;
                end
                I2=K+T;
                if I2>4
                    I2=I2-4;
                end
                SP(T,:) = mean([mBB(:,I1)';mBB(:,I2)']);
            end         
            
%             figure
%             axis equal
%             hold on
%             plot(xy2new(:,1), xy2new(:,2), Style,'MarkerSize',MarkerSize,'Color',[0 0 0]);
%             plot(SP(1,1), SP(1,2), '*', 'MarkerSize', 15,'Color',[1 0 0]);
%             plot(SP(3,1), SP(3,2), '*', 'MarkerSize', 15,'Color',[1 0 0]);
%             plot(SP(2,1), SP(2,2), '*', 'MarkerSize', 15,'Color',[0 0 1]);
%             plot(SP(4,1), SP(4,2), '*', 'MarkerSize', 15,'Color',[0 0 1]);
%             plot(mBB(1,[1:end 1]),mBB(2,[1:end 1]),'r');
%             title({'The bounding box is shown in red: to VALIDATE it press RETURN';'otherwise SELECT the bounding box starting from the bottom-left and clockwise'});
            m1 = (SP(1,2)-SP(3,2))/(SP(1,1)-SP(3,1));
            b1 = -m1*SP(3,1)+SP(3,2);
            xretta1=min(xy2new(:,1)): (max(xy2new(:,1))-min(xy2new(:,1)))/100 : max(xy2new(:,1)); 
            yretta1=m1*xretta1+b1*ones(1,size(xretta1,2));
%             plot(xretta1, yretta1, 'r');
            m2 = (SP(2,2)-SP(4,2))/(SP(2,1)-SP(4,1));
            b2 = -m2*SP(4,1)+SP(4,2);
            yretta2=min(xy2new(:,2)): (max(xy2new(:,2))-min(xy2new(:,2)))/100 : max(xy2new(:,2)); 
            xretta2=(1/m2)*(yretta2-b2*ones(1,size(yretta2,2))); %%(y-b)/m=x
%             plot(xretta2, yretta2, 'b');
           
%             [xx, yy] = ginput(4);
%             if (size(xx, 1)==4 && size(yy, 1)==4)
%                 newBB = [xx, yy];
%                 for T=1:3
%                     SP(T,:) = (1/2)*[newBB(T,1)+newBB(T+1,1) newBB(T,2)+newBB(T+1,2)];
%                 end
%                 SP(4,:) = (1/2)*[newBB(4,1)+newBB(1,1) newBB(4,2)+newBB(1,2)];
%                 m = (SP(1,2)-SP(3,2))/(SP(1,1)-SP(3,1));
%                 b = -m*SP(3,1)+SP(3,2);
%                 alpha = atan(m);
%                 if ((SP(3,1)-SP(1,1))<0)
%                     alpha=alpha+pi;
%                 end
%                 R = [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
%                 newRR2{1,i} = R;
%                 xy2new=xy2new*R;
%                 SP=SP*R;
%                 xretta=min(xy2new(:,1)): (max(xy2new(:,1))-min(xy2new(:,1)))/100 : max(xy2new(:,1)); 
%                 yretta=b*ones(1,size(xretta,2));
%          
%                 figure
%                 hold on; axis equal;
%                 plot(xretta, yretta, 'r');
%                 plot(xy2new(:,1), xy2new(:,2), Style,'MarkerSize',MarkerSize,'Color',[0 0 1]);
%                 plot(SP(:,1), SP(:,2), '*', 'MarkerSize', 15,'Color',[1 0 0]);
%             else
%                 newRR2{1,i} = eye(2);
%             end

            %riallineo
%             MMx{1,i} = SP(4,1);
%             MMy{1,i} = SP(4,2);
%             xy2new(:,1)=xy2new(:,1)-SP(4,1);
%             xy2new(:,2)=xy2new(:,2)-SP(4,2);
%             SP(:,1) = SP(:,1)-SP(4,1);
%             SP(:,2) = SP(:,2)-SP(4,2);
            
            % stretch and transform to polar coordinates
%             D1 = norm(SP(2,:)-SP(4,:));
%             meanY = (1/2)*SP(1,2);
%             stretchCoeff = sqrt((D1^2-meanY^2)/SP(1,1)^2);
%             SC{1,i} = stretchCoeff;
%             xy2new(:,1) = stretchCoeff*xy2new(:,1);
            %
            thetaRho = zeros(size(xy2new,1),2);
            for k=1:size(xy2new,1)
                rho = sqrt(xy2new(k,1)^2+xy2new(k,2)^2);
                thetaRho(k,2) = rho;
                if xy2new(k,1)>0
                    theta = atan(xy2new(k,2)/xy2new(k,1));
                else
                    theta = atan(xy2new(k,2)/xy2new(k,1))+pi;
                end
                thetaRho(k,1) = theta;
            end;

            %rotate to default position
            [alpha,II] = min(thetaRho(:,1));
            MMM = max(thetaRho(:,2));
            M= [cos(-alpha) sin(-alpha); -sin(-alpha) cos(-alpha)];
            newRR3{1,i} = M;
            xy2new = xy2new*M;
            thetaRho(:,1) = thetaRho(:,1)-alpha;
          
            D1 = norm(SP(1,:)-SP(3,:));
            D2 = norm(SP(2,:)-SP(4,:));
            A = abs(2*D1-D2)/2;
            B = abs(D2-D1)/pi; 

            PtMin = [0.5*A 0.5*B]; 
            PtMax = [1.5*A 1.5*B];
            Eps = [max((1/20)*(PtMax(1)-PtMin(1)),10^(-4)) max((1/50)*(PtMax(2)-PtMin(2)),10^(-4))]; 
            
            % write parameters limits and epsilon on a file
            nameFile = strcat('parameters_', int2str(i) , '.cocoa5');
            fileID = fopen(nameFile,'w');
            fprintf(fileID, '%9s\n','PtMin:=');
            fprintf(fileID,'[%6.4f, %6.4f] \n', PtMin);
            fprintf(fileID, '%3s\n',';');
            fprintf(fileID, '%9s\n','PtMax:=');
            fprintf(fileID,'[%6.4f, %6.4f] \n', PtMax);
            fprintf(fileID, '%3s\n',';');
            fprintf(fileID, '%9s\n','Eps:=');
            fprintf(fileID,'[%6.4f, %6.4f] \n', Eps);
            fprintf(fileID, '%3s\n',';');
            fclose(fileID);
            
            % write cluster points on a file
            nameFile = strcat('clusterPoints_', int2str(i) , '.cocoa5');
            fileID = fopen(nameFile,'w');
            fprintf(fileID, '%9s\n','Points:=[');
            fprintf(fileID,'[%14.12f, %14.12f], \n', thetaRho(1:(size(xy2new,1)-1),:)');
            fprintf(fileID,'[%14.12f, %14.12f] \n', thetaRho(size(xy2new,1),:)');
            fprintf(fileID, '%3s\n','];');
            fclose(fileID);
        end
    end
    
    s = size(UsedClusters,2);
    indPts = cell(1, s);
    axis1 = cell(1,s);
    axis2 = cell(1,s);
    currentFolder = pwd;
    nameFile = strcat('currentFolder.cocoa5');
    fileID = fopen(nameFile,'w');
    stringa = strcat('str:=',currentFolder,";");
    fprintf(fileID, '%s\n', stringa);
    fclose(fileID);
    for k=1:s
        i=UsedClusters(k);
        nameFile = strcat('input4CoCoA.cocoa5');
        fileID = fopen(nameFile,'w');
        stringa1 = strcat('Source "', currentFolder, '/clusterPoints_', int2str(i), '.cocoa5";');
        stringa2 = strcat('Source "', currentFolder, '/parameters_', int2str(i), '.cocoa5";');
        fprintf(fileID, '%s\n', stringa1);
        fprintf(fileID, '%s\n', stringa2);
        fclose(fileID);
       
        % run Hough transform (recognitionAlgorithm) on file
        command = './CoCoA-5/cocoa5 recognitionAlgorithm_spiral.cocoa5';
        [status,cmdout] = system(command);
        M = importdata('RelMaxValue');
        MaxCoords{1,i} = importdata('MaxCoords');
        RelMaxValue{1,i} = M;
        
        % delete the files
        command = strcat('rm clusterPoints_', int2str(i), '.cocoa5;');
        [status,cmdout] = system(command);
        command = strcat('rm parameters_', int2str(i), '.cocoa5;');
        [status,cmdout] = system(command);
        command = 'rm input4CoCoA.cocoa5';
        [status,cmdout] = system(command);
        command = 'rm RelMaxValue';
        [status,cmdout] = system(command);
        command = 'rm MaxCoords';
        [status,cmdout] = system(command);
    end  
    
    %delete the file currentFolder
    command = 'rm currentFolder.cocoa5';
    [status,cmdout] = system(command);
        
    valori = zeros(numCluster,1);
    for jj=1:numCluster
        if ~(isempty(RelMaxValue{1,jj}))
            valori(jj,1)=RelMaxValue{1,jj};            
        end
    end
    [Val, II] = sort(valori);
  
    
    % plot tutti i risultati - punti del cluster scelto con geometric petal
    curvePars = [];
    kk=numCluster;
    while (Val(kk,1)>0 & kk>0)
        i = II(kk,1);
        xyzOrig=punti(IDX==i,:);
        xyz = bsxfun(@minus, xyzOrig, MM{1,i});
        xyz=xyz*RR{1,i};
        xy=[xyz(:,1) xyz(:,2)];
        xy = xy*newRR{1,i};
        %xy = xy*newRR2{1,i};
        %xy(:,1) = xy(:,1) - MMx{1,i};   
        %xy(:,2) = xy(:,2) - MMy{1,i};
        xy = xy*newRR3{1,i};
        jj=1;
        MMCC = MaxCoords{1,i};
        MC = MMCC(jj,:);
        a = MC(1);  
        b=MC(2);
        curvePars = [curvePars; [a b]];
        % select the points close to the recognized curve......
        % build sample points on the recognized curve
        puntiCurva = zeros(501,2);  
        for j=0:500
            theta=(j/500)*2*pi;
            rho=a+b*theta;
            xx = rho*cos(theta);
            yy = rho*sin(theta);
            puntiCurva(j+1,1)=xx;
            puntiCurva(j+1,2)=yy;
        end
       
        %figure
        %axis equal
        %hold on
        [I,dist] = knnsearch(puntiCurva,xy,'k', 1);
        outCross = find(dist<eps);
        outNotCross = find(dist>eps);
        %plot(xy(outCross,1),xy(outCross,2),'*k')
        %plot(xy(outNotCross,1),xy(outNotCross,2),'*r')
        %h = ezplot(@(x,y)geometricPetal(x,y,N,a^2,c^2), [min(xy(:,1)), max(xy(:,1)), 0, max(xy(:,2))])
        %set(h, 'Color', 'g','LineWidth',2)
        %title([])
        %plot(puntiCurva(:,1),puntiCurva(:,2),'.b')
        
        selectedPoints=[xyzOrig(outCross(1),:)];
        for k=2:size(outCross,1)    
            selectedPoints = [selectedPoints; xyzOrig(outCross(k,1),:)];                       
        end

%         % scrivo su file
%         nameFile = strcat('selectedPointsMax_', int2str(numCluster-kk+1), '.off');
%         fileID = fopen(nameFile,'w');
%         fprintf(fileID, '%-4s\n','COFF');
%         fprintf(fileID, '%-5d 0 0\n',size(selectedPoints,1));
%         fprintf(fileID,'%6.4f %6.4f %6.4f 0 0 255 255\n', selectedPoints');
%         fclose(fileID);
        
        resIndici = zeros(1,size(outCross,1));
        IndiciCluster = find(IDX==i);
        for k=1:size(outCross,1)
            resIndici(k) = select(IndiciCluster(outCross(k,1),1));
        end
        indPts{1, numCluster-kk+1} = resIndici;
        
        % salient points 
        theta = pi;
        rho = a+b*theta;
        XX1 = -rho;
        YY1 = 0;
        theta = 0;
        rho = a;
        XX2 = rho;
        YY2 = 0;
        theta = (3/2)*pi;
        rho = a+b*theta;
        XX3 = rho*cos(theta);
        YY3 = rho*sin(theta);
        theta = (1/2)*pi;
        rho = a+b*theta;;
        XX4 = rho*cos(theta);
        YY4 = rho*sin(theta);
        salientPts=[XX1 YY1; XX2 YY2; XX3 YY3; XX4 YY4];
        salientPts = salientPts*(newRR3{1,i})^(-1);
        salientPts = salientPts*(newRR{1,i})^(-1);
        salientPts3= [salientPts mean(xyz(:,3))*ones(4,1)];  %%better 0??
        salientPts3 = salientPts3*(RR{1,i})^(-1);       
        Asse1.Pt1 = salientPts3(1,:);
        Asse1.Pt2 = salientPts3(2,:);
        axis1{1, numCluster-kk+1} = Asse1;
        Asse2.Pt1 = salientPts3(3,:);
        Asse2.Pt2 = salientPts3(4,:);
        axis2{1, numCluster-kk+1} = Asse2;
        
%         nameFile = strcat('salientPoints_', int2str(numCluster-kk+1), '.off');
%         fileID = fopen(nameFile,'w');
%         fprintf(fileID, '%-4s\n','COFF');
%         fprintf(fileID, '%-5d 0 0\n', 4);
%         fprintf(fileID,'%6.4f %6.4f %6.4f 0 0 255 255\n', salientPts3');
%         fclose(fileID);

        kk= kk-1;
    end   
    %delete all the created files
  
end
 
 

