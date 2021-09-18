function printStiff(Solution,fid)

    if isfield(Solution,'data')
        printCTitle(fid,'S T I F F N E S S    M A T R I C E S',96);

        nElem = size(Solution.data.ke_local,3);
        for i = 1 : nElem
            printLine(fid,96);
            fprintf(fid,' lement No = %d \n',i);
            fprintf (fid,' - Local Stiffness Matrix: \n');
            disp_mat(Solution.data.ke_local,1,fid)
            printSLine(fid,96);

            fprintf (fid,' - Transformation Matrix: \n');
            disp_mat(Solution.data.T6,2,fid)        
            printSLine(fid,96);

            fprintf (fid,' - GlobalStiffness Matrix: K=T`kT = \n');    
            disp_mat(Solution.data.ke_global,3,fid)
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


function disp_mat(mat,mode,fid)
format compact

if mode==1 
    format shortEng
    fprintf(fid,'  +----                           ----+   +----                                               ----+ \n');
    fprintf(fid,'  | k11   k12   k13   k14   k15   k16 |   | EA/L   0         0        -EA/L     0         0       | \n');
    fprintf(fid,'  | k21   k22   k23   k24   k25   k26 |   | 0      12EI/L^3  6EI/L^2   0       -12EI/L^3  6EI/L^2 | \n');
    fprintf(fid,'  | k31   k32   k33   k34   k35   k36 |   | 0      6EI/L^2   4EI/L     0       -6EI /L^2  2EI/L   | \n');
    fprintf(fid,'  | k41   k42   k43   k44   k45   k46 | = |-EA/L   0         0         EA/L     0         0       |=\n');
    fprintf(fid,'  | k51   k52   k53   k54   k55   k56 |   | 0     -12EI/L^3 -6EI/L^2   0        12EI/L^3 -6EI/L2  | \n');
    fprintf(fid,'  | k61   k62   k63   k64   k65   k66 |   | 0      6EI/L^2   2EI/L     0       -6EI /L^2  4EI/L   | \n');
    fprintf(fid,'  +----                           ----+   +----                                               ----+ \n');
    for i = 1 : 6
        for j = 1 : 6
            fprintf(fid,'%14.3e\t',mat(i,j));
        end
        fprintf(fid,'\n');
    end

elseif mode==2
    format short g
    c=mat(1,1);
    s=mat(1,2);
    fprintf(fid,'  +---                  ---+   +---                               ---+\n');
    fprintf(fid,'  |  c   s   0   0   0   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',c ,s ,0,0,0,0);
    fprintf(fid,'  | -s   c   0   0   0   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',-s,c ,0,0,0,0);
    fprintf(fid,'  |  0   0   1   0   0   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,1,0,0,0);
    fprintf(fid,'  |  0   0   0   c   s   0 | = | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,0,c,s,0);
    fprintf(fid,'  |  0   0   0  -s   c   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,0,-s,c,0);
    fprintf(fid,'  |  0   0   0   0   0   1 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,0,0,0,1);
    fprintf(fid,'  +---                  ---+   +---                               ---+\n');

elseif mode==3
    for i = 1 : 6
        for j = 1 : 6
            fprintf(fid,'%14.3e\t',mat(i,j));
        end
        fprintf(fid,'\n');
    end
end


fprintf('\n')
format long 
end