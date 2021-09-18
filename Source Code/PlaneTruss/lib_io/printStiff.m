function printStiff(Solution,fid)

    if nargin < 2
        fid = 1;
    end
    
    if isfield(Solution,'data')
        % element stiffness
        if isfield(Solution.data,'elementStiff')
            printCTitle(fid,'S T I F F N E S S     M A T R I C E S')
            printLine(fid);
    
            [n,m,nElem] = size(Solution.data.elementStiff);
            for iElem = 1 : nElem
                stif = Solution.data.elementStiff(:,:,iElem);
                fprintf(fid, 'Element No. [%2d]\n',iElem);
                printSLine(fid);
                for i = 1 : n
                    for j = 1 : m
                        fprintf(fid,'%14.3e\t',stif(i,j));
                    end
                    fprintf(fid,'\n');
                end
            end
        end
        printLine(fid);
        
        % global stiffness
        if isfield(Solution.data,'globalStiff')
            fprintf(fid, 'Global Stiffness Matrix\n');
            stif = Solution.data.globalStiff;
            [n,m] = size(stif);
            for i = 1 : n
                for j = 1 : m
                    fprintf(fid,'%14.3e\t',stif(i,j));
                end
                fprintf(fid,'\n');
            end
        end
        printLine(fid);
        
        % global forces
        if isfield(Solution.data,'globalForce')
            fprintf(fid, 'Global Nodal Forces\n');
            F = Solution.data.globalForce;
            [n,m] = size(F);
            for i = 1 : n
                fprintf(fid,'[%3d] \t',i);
                for j = 1 : m
                    fprintf(fid,'%14.3e\t',F(i,j));
                end
                fprintf(fid,'\n');
            end
        end
        
        printDLine(fid);
    end
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
