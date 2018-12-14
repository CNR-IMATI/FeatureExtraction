function [n_facets,filename,facets]=read_faceting_file(fct_path)

% read_faceting_file:  Reads .fct file that contains the results of the
%                      faceting.
%
% [n_facets,filename,facets]=read_faceting_file(fct_path)
% 
%
% Input:    
%  'fct_path'   is the full path of the .fct file containing the
%               faceting results
% 
% Output: 
%   'n_facets'           is the number of Facets as mentioned in the .fct
%                        file's header.                        
%   'filename'           is the fragment's name as mentioned in the .fct 
%                        file's header.
%   'facets'             is the the faceting results where facet{i} is a 
%                        struct including: id, type, border and vertices of
%                        facet i.
%   'facets{i}.id'       is the id of facet i (ex: 'A', 'B', ...)
%   'facets{i}.type'     is the type of facet i:
%                        external/internal or fracture.
%   'facets{i}.border'   is the indices of vertices reprsenting the border
%                        of facet i.
%   'facets{i}.vertices' is the indices of vertices reprsenting inner
%                        vertices of facet i.


% 
%  Copyright (c) 2017 UvA


fileID = fopen(fct_path,'r');

%%Reading the file header
Intro = textscan(fileID,'%s',3,'Delimiter','\n');
temp=strsplit(char(Intro{1}{1}),':');
filename=char(temp{end}(2:end));
temp=strsplit(char(Intro{1}{2}),':');
n_facets=str2num(temp{end});


facets=[];
formatSpec2= '%d ';
i=1;
while (~feof(fileID)) 
    %%Reading Facet ID
    InputText = textscan(fileID,'%s',1,'delimiter','\n');
    temp=strsplit(char(InputText{1}));
    facets{i}.id=char(temp{end});
    
    %%Reading Facet Type
    InputText= textscan(fileID,'%s',1,'delimiter','\n');
    temp=strsplit(char(InputText{1}));
    facets{i}.type=char(temp{end});

    %%Reading Facet Border
    InputText = textscan(fileID,'%s',1,'delimiter','\n');
    facets{i}.border = fscanf(fileID,formatSpec2);
    
    %%Reading Facet Vertices
    InputText = textscan(fileID,'%s',1,'delimiter','\n');
    facets{i}.vertices = fscanf(fileID,formatSpec2);
    i=i+1;
    
end
facets=facets';
fclose(fileID);
disp(sprintf('%s is successfully read.',fct_path));
