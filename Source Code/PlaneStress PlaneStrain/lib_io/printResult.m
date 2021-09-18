function printResult(Solution,fid)
    
    if nargin < 2
        fid = 1;
    end
    BULLET_CHAR = char(15);
    
    fprintf(fid,'\n\n');
    printDLine(fid)
    printCTitle(fid,'O U T P U T     S U M M A R Y')
    printDLine(fid)
    
    if isfield(Solution,'info')
       fprintf(fid,'%c About: \n',BULLET_CHAR);
       fprintf(fid,'%s\n', Solution.info);
    end
    printLine(fid)
    
    
    if isfield(Solution, 'nodalDisplacements')
        u = Solution.nodalDisplacements;
        [nNode,nDof] = size(u);
        printCTitle(fid,'Nodal Displacements');
        printSLine(fid);
        for i = 1 : nNode
            formatStr = repmat('%10.5f\t',1,nDof);
            fprintf(fid,['[%3d]   ',formatStr,'\n'],i,u(i,:));
        end
    end
    printLine(fid)
    if isfield(Solution, 'reactionForces')
        rf = Solution.reactionForces;
        [nNode,nDof] = size(rf);
        printCTitle(fid,'Reaction Forces');
        printSLine(fid);
        for i = 1 : nNode
            formatStr = repmat('%10.5f\t',1,nDof);
            fprintf(fid,['[%3d]   ',formatStr,'\n'],i,rf(i,:));
        end
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
