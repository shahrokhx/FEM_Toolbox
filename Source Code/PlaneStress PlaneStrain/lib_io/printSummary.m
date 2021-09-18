function printSummary(Model,fid)
% PURPOSE: reads input files and generates a structure including all data    
%
% INPUT(S): 
%   - Model: a structure including the following fields
%   - fid  : file id to store
%
% OUTPUT(S):
%  
%
% USAGE:
%   >> printSummary(Model)       print input details on command window
%   >> printSummary(Model,fid)   print input details on a file with fid
%                                identifier

% DEVELOPMENT HISTORY:
%   [2018, Oct, 29] Shahrokh Shahi -- initial development
%   [2018, XXX, XX] Shahrokh Shahi -- 

    if nargin < 2 
        fid = 1;
    end
    BULLET_CHAR = char(15);
    
    printDLine(fid)
    tText = 'A    S U M M A R Y    O F    T H E    I N P U T    M O D E L';
    printCTitle(fid,tText)
    printDLine(fid)


    if isfield(Model,'info')
       fprintf(fid,'%c About: \n',BULLET_CHAR);
       fprintf(fid,'%s\n', Model.info);
    end
    if isfield(Model,'analysisType')
       fprintf(fid,'\n%c Structure Type: \n',BULLET_CHAR);
       fprintf(fid,'%s\n', Model.analysisType);
    end
    printLine(fid)

    
    fprintf(fid,'%c Control Variables: \n',BULLET_CHAR);
    fprintf(fid,'   - Dimension:             %dD\n' ,Model.nDim);
    fprintf(fid,'   - Number of DOF(s):      %d \n' ,Model.nDof);
    fprintf(fid,'\n');
    fprintf(fid,'   - Number of Nodes:    %4d \n' ,Model.nNode);
    fprintf(fid,'   - Number of Elements: %4d \n' ,Model.nElem);
    fprintf(fid,'   - Nodes per Elements:    %d\n',Model.nElemNode);
    fprintf(fid,'\n');
    fprintf(fid,'   - Number of Sections:    %d \n',Model.nSection);
    printDLine(fid)

    
    printCTitle(fid,'N O D A L    I N F O R M A T I O N')
    printLine(fid)
    headerTxt = ['Node     Coordinates          Nodal Loads       ',...
                 '         Nodal Restraints    \n'];
    fprintf(fid,headerTxt);
    printSLine(fid)
    switch Model.nDim
        case 1 
            headerTxt1 =' ID          X               ';
            headerTxt2 =' Fx                      ';
            headerTxt3 =' ux          \n';
        
        case 2
            headerTxt1 =' ID       X     Y            ';
            headerTxt2 =' Fx       Fy            ';
            headerTxt3 =' ux   uy     \n';
            
        case 3
            headerTxt1 =' ID     X    Y    Z          ';
            headerTxt2 =' Fx       Fy       Fz   ';
            headerTxt3 =' ux   uy   uz \n';
    end
    fprintf(fid,[headerTxt1,headerTxt2,headerTxt3]);
    printLine(fid)
    
    % nodal forces to displaying
    F = zeros(Model.nNode,Model.nDim);
    for i = 1 : Model.nNodalForce
        node = Model.loading.nodalForces(i,1);
        F(node, 1:Model.nDim) = Model.loading.nodalForces(i,2:end);
    end
        
    % printing table contents
    for i = 1 : Model.nNode
        coord = Model.geometry.coordinates(i,:);
        fprintf(fid,'%2d  ',i);
        for j = 1 : length(coord)
            fprintf(fid,'%7.2f  ',coord(j));
        end
        fprintf(fid,'\t');
        for j = 1 : size(F,2)
            fprintf(fid,'%10.2f  ',F(i,j));
        end
        
        fprintf(fid,'\t');
        index = find(Model.boundary(:,1)==i);
        dx='-';dy='-';dz='-';bc='';
     
        if ~isempty(index)
            for k = 1 : length(index)
                if Model.boundary(index(k),2)==1
                    dx=num2str(Model.boundary(index(k),3));
                elseif Model.boundary(index(k),2)==2
                    dy=num2str(Model.boundary(index(k),3));
                elseif Model.boundary(index(k),2)==3
                    dz=num2str(Model.boundary(index(k),3));
                end
             end
             if length(index)==1
                 bc = 'Roller ';
             elseif length(index)==2
                 bc = 'Pin';
             elseif length(index)==3
                 bc = 'Fixed';
             end
        end
        presc = [dx, dy, dz]; %prescribed values
        for j = 1 : Model.nDof
            fprintf(fid,'   %s',presc(j));
        end
        fprintf(fid,'       %s\n',bc);
    end
    printDLine(fid);
    
    printCTitle(fid,'E L E M E N T S     I N F O R M A T I O N')
    printLine(fid)
    fprintf(fid,' ID\t ');
    for i = 1 : Model.nElemNode
        fprintf(fid,'Node%d\t',i);
    end
    fprintf(fid,'     Section No.     Material Properties\n');
    printLine(fid);
    
    for i = 1 : Model.nElem
        fprintf(fid,'%2d\t ',i);
        for j = 1 : Model.nElemNode
            fprintf(fid,'%2d\t\t',Model.geometry.elements(i,j));
        end
        fprintf(fid,'\t');
        fprintf(fid,'\t%2d\t\t\t[',Model.elemSectId(i));
        for j = 1 : Model.nMaterialProp
            fprintf(fid,'%g\t\t',Model.sections(Model.elemSectId(i),j));
        end
        fprintf(fid,']\n');
    end
    
    printDLine(fid);
    
end


% Helper Functions
function printLine(fid,n)
    if nargin < 2
        n = 80; 
    end
    if nargin < 1
        fid = 1;
    end
    fprintf(fid,'%s\n',repmat('_',1,n));
end

function printDLine(fid,n)
    if nargin < 2
        n = 80; 
    end
    if nargin < 1
        fid = 1;
    end
    fprintf(fid,'%s\n',repmat('=',1,n));
end

function printSLine(fid,n)
    if nargin < 2
        n = 80; 
    end
    if nargin < 1
        fid = 1;
    end
    fprintf(fid,'%s\n',repmat('-',1,n));
end
% printing centered text title
function printCTitle(fid,text,n)
    if nargin < 3
        n = 80;
    end
    if nargin < 2
        text =''; 
    end
    if nargin < 1
        fid = 1;
    end
    tLen = length(text);
    space = floor((n-tLen)/2);
    fprintf(fid,'%s\n',[repmat(' ',1,space),text]);
end