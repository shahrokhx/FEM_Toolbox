%===================== Assemble the global traction vector =============
%
function r = globaltraction(ncoord,ndof,nnodes,ndload,coords,nelnodes,elident,connect,dloads,dofs)

   r = zeros(ndof*nnodes,1);
%    traction = zeros(ndof,1);
   nTraction = length(dloads);
   for load = 1  :  nTraction  % ndload
%
%     Extract the coords of the nodes on the appropriate element face
%
      traction = dloads(load);
      
      lmn  = traction.element; % dloads(1,load);
      face = traction.face;    % dloads(2,load);
      n = nelnodes(lmn);
      ident = elident(lmn);
      nfnodes = nfacenodes(ncoord,n,ident,face); 
      nodelist = facenodes(ncoord,n,ident,face);     
      lmncoord = zeros(ncoord,nfnodes);
      for a = 1:nfnodes
        for i = 1:ncoord
          lmncoord(i,a) = coords(i,connect(nodelist(a),lmn));  %,dloads(1,load)));
        end
        for i = 1:ndof
%         lmndof(i,a) = dofs(ndof*(connect(nodelist(a),dloads(1,load))-1)+i);
          lmndof(i,a) = dofs(ndof*(connect(nodelist(a),lmn)-1)+i);
        end
      end
%
%    Compute the element load vector
%

%      for i = 1:ndof
%        traction(i) = dloads(i+2,load);
%      end

     loadEq  = traction.eq;
     
     rel = eldload(ncoord,ndof,nfnodes,ident,lmncoord,loadEq); %traction);
%
%    Assemble the element load vector into global vector
%
     for a = 1:nfnodes
       for i = 1:ndof
%          rw = (connect(nodelist(a),dloads(1,load))-1)*ndof+i;
         rw = (connect(nodelist(a),lmn)-1)*ndof+i;
         r(rw) = r(rw) + rel((a-1)*ndof+i);
       end
     end

   end
end  