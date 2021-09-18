function print_results(outfile, ...
    nprops,materialprops,ncoord,ndof,nnode,coords, ...
    nelem,maxnodes,connect,nelnodes,elident, ...
    nfix,fixnodes,ndload,dloads,dofs)
%      Print nodal displacements, element strains && stresses:a file
    fprintf(outfile,'\n\n');
    printDLine(outfile)
    printCTitle(outfile,'O U T P U T     S U M M A R Y')
    printDLine(outfile)
    
    printCTitle(outfile,'N O D A L    D I S P L A C E M E N T S');
    printLine(outfile)
   if (ndof == 2) 
     fprintf(outfile,' Node      Coords         u1       u2 \n');
     for i = 1:nnode
      fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f\n', ...
                               i,coords(1,i),coords(2,i),dofs(2*i-1),dofs(2*i));
     end
   elseif (ndof == 3) 
     fprintf(outfile,' Node            Coords            u1       u2       u3 \n');
     for i = 1:nnode
      fprintf(outfile,'%3d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n', ...
                    i,coords(1,i),coords(2,i),coords(3,i),dofs(3*i-2),dofs(3*i-1),dofs(3*i));
     end
   end
   
   printDLine(outfile)
   printCTitle(outfile,'S T R E S S E S    A N D   S T R A I N S') 

   lmncoord = zeros(ncoord,maxnodes);
   displacements = zeros(ndof,maxnodes);

%
%   Loop over all the elements
%
   for lmn = 1:nelem
    printLine(outfile)
    fprintf(outfile,' \n Element: %d ',lmn);
    if (ncoord == 2)   
    fprintf(outfile,'  \n int pt    Coords          e_11      e_22     e_12      s_11       s_22      s_12 \n');

    elseif (ncoord == 3) 
    fprintf(outfile,'\n int pt         Coords            e_11      e_22     e_33      e_12       e_13      e_23      s_11      s_22      s_33      s_12      s_13      s_23 \n');
    end
%
%   Extract coords of nodes, DOF for the current element
%
      for a = 1:nelnodes(lmn)
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(a,lmn));
        end
        for i = 1:ndof
          displacements(i,a) = dofs(ndof*(connect(a,lmn)-1)+i);
        end
      end
      n = nelnodes(lmn);
      ident = elident(lmn);
 
      npoints = numberofintegrationpoints(ncoord,n);
      dNdx = zeros(n,ncoord);
      dxdxi = zeros(ncoord,ncoord);
      strain = zeros(ndof,ncoord);
      xi = zeros(ncoord,1);
      x = zeros(ncoord,1);
%
%  Set up integration points 
%
      xilist = integrationpoints(ncoord,n,npoints);
%
%  Loop over the integration points
%
     for intpt = 1:npoints

%     Compute shape functions && derivatives wrt local coords
%
       for i = 1:ncoord
         xi(i) = xilist(i,intpt);
       end
       N = shapefunctions(n,ncoord,ident,xi);      
       dNdxi = shapefunctionderivs(n,ncoord,ident,xi);
%
%     Compute the coords of the integration point
%
      for i = 1:ncoord
        x(i) = 0.;
        for a = 1:n
          x(i) = x(i) + lmncoord(i,a)*N(a);
        end
      end
%
%     Compute the jacobian matrix && its determinant
%
      for i = 1:ncoord
        for j = 1:ncoord
          dxdxi(i,j) = 0.;
          for a = 1:n
            dxdxi(i,j) = dxdxi(i,j) + lmncoord(i,a)*dNdxi(a,j);
          end
        end
      end

      dxidx = inv(dxdxi);
%
%     Convert shape function derivatives:derivatives wrt global coords
%
      for a = 1:n
        for i = 1:ncoord
          dNdx(a,i) = 0.;
          for j = 1:ncoord
            dNdx(a,i) = dNdx(a,i) + dNdxi(a,j)*dxidx(j,i);
          end
        end
      end
%
%     Compute the (infinitesimal) strain by differentiating displacements
%
      for i = 1:ncoord
         for j = 1:ncoord
            strain(i,j) = 0.;
            for a = 1:n
              strain(i,j) = strain(i,j) + 0.5*(displacements(i,a)*dNdx(a,j)+displacements(j,a)*dNdx(a,i));
            end
         end
      end

      stress = materialstress(ndof,ncoord,strain,materialprops);

      if (ncoord == 2) 

      fprintf(outfile,'%5d %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n', ...
        intpt,x(1),x(2),strain(1,1),strain(2,2),strain(1,2),stress(1,1),stress(2,2),stress(1,2));


      elseif (ncoord == 3) 

      fprintf(outfile,'%5d %7.4f %7.4f %7.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f \n',...
              intpt,x(1),x(2),x(3), ...
              strain(1,1),strain(2,2),strain(3,3),strain(1,2),strain(1,3),strain(2,3), ...
              stress(1,1),stress(2,2),stress(3,3),stress(1,2),stress(1,3),stress(2,3));
      end
     end
     
   end
   printDLine(outfile)
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