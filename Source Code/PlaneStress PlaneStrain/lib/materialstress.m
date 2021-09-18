%================= Material Stress ==================================
%
%   Computes stress sigma_{ij} given strain epsilon_{ij}
%
function stress = materialstress(ndof,ncoord,strain,materialprops)

   C = materialstiffness(ndof,ncoord,strain,materialprops);
   stress = zeros(ndof,ncoord);
   for i = 1 : ndof
     for j = 1 : ncoord
        for k = 1 : ndof
          for l = 1: ncoord
            stress(i,j) = stress(i,j) + C(i,j,k,l)*strain(k,l);
          end
        end
     end
   end
end