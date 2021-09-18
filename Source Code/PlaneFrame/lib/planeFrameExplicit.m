function Solution = planeFrameExplicit(Model)

Model = ModelRevise(Model);

Solution.info = ['Solution obtained on ', date];

FirstNode  = Model.geometry.elements(:,1);
SecondNode = Model.geometry.elements(:,2);
ifix       = Model.bound.ifix;

% processing element load information
ndof = Model.nNode * Model.nDof;
K    = zeros(ndof,ndof); %initialize the global stiffness matrix of size ndof
Feq  = zeros(ndof,1); %nodal force vector due to element loads
ke_local (:,:,Model.nElem) = zeros(6);
ke_global(:,:,Model.nElem) = zeros(6);
T6(:,:,Model.nElem)        = zeros(6);

for e = 1 : Model.nElem
    c  = Model.geometry.Cx(e);
    s  = Model.geometry.Cy(e);
    le = Model.geometry.L(e);

    ae = Model.properties.A(e) * Model.properties.E(e); %axial rigidity 
    ei = Model.properties.E(e) * Model.properties.Iz(e); %flexural rigidity
    %Degrees of freedom for element e with nodes m and n are 3m-2, 3m-1 ,3m, 3n-2, 3n-1 ,3n
    %These degrees of freedom are stored in the vector {dof} as shown below
    dof = [3*FirstNode(e)-2, 3*FirstNode(e)-1, 3*FirstNode(e), ...
           3*SecondNode(e)-2, 3*SecondNode(e)-1, 3*SecondNode(e)];
       
    % rotation transformation matrix for element e
    T=[c    s    0    0     0     0;
      -s    c    0    0     0     0;
       0    0    1    0     0     0;
       0    0    0    c     s     0;
       0    0    0   -s     c     0;
       0    0    0    0     0     1;];
   
    % element stiffness matrix w.r.t local axes
    estif = [ ae/le     0             0          -ae/le     0            0; 
               0       12*ei/le^3    6*ei/le^2       0    -12*ei/le^3   6*ei/le^2;
               0       6*ei/le^2     4*ei/le         0    -6*ei/le^2    2*ei/le ;
            -ae/le      0             0           ae/le      0           0; 
               0      -12*ei/le^3   -6*ei/le^2       0     12*ei/le^3  -6*ei/le^2;
               0       6*ei/le^2     2*ei/le         0    -6*ei/le^2    4*ei/le ];
               
    if Model.geometry.hingIndex(SecondNode(e)) == 1   % if the right node is a hinge
        % Modifications due to the hing on the right hand
        estif( 6   , :   ) = 0;
        estif( :   , 6   ) = 0;
        estif([2,5],[2,5]) = (1/4) * estif([2,5],[2,5]);
        estif([2,5], 3   ) = (1/2) * estif([2,5], 3   );
        estif( 3   ,[2,5]) = (1/2) * estif( 3   ,[2,5]);
        estif( 3   , 3   ) = (3/4) * estif( 3   , 3   );
    end
     
     
     ke=(T'*estif*T) ;  %element stiffness matrix with reference to global axes
     
     ke_local(:,:,e) = estif;
     T6(:,:,e)= T;
     ke_global(:,:,e) = ke;
     
     
     K(dof,dof)=K(dof,dof)+ke;  %assembly of global stiffness matrix K
     %{Feq} is vector of equivalent nodal forces due to element loads
     Feq(dof) = Feq(dof) + T'*Model.loading.eforce(:,e);     
end

Fc = Model.loading.F+Feq; %combined force vector {Fc} is obtained by adding vector of nodal loads {F} to {Feq}

Solution.data.ke_local = ke_local;
Solution.data.ke_global = ke_global;
Solution.data.T6 = T6;
Solution.data.F_global = Fc;

%-------------------------------------------------------------------------
%applying boundary conditions (specified by restraints along X and Y axes)
%on global stiffness matrix and global force vector
K(ifix,:)=0; %Elements in rows with d.o.f corresponding to ifix are made 0                 
K(:,ifix)=0; %Elements in columns with d.o.f corresponding to ifix are made 0                    
K(ifix,ifix)=eye(length(ifix));%Diagonal elements with d.o.f corresponding to ifix are made 1     
Fc(ifix) = 0; %corresponding columns of Global force vector are also made zero
%---------------------------------------------------------------------
%obtaining vector of nodal displacements along global axes
u=inv(K)*Fc;

Solution.data.u = u;



%---------------------------------------------------------------------
%P O S T P R O C E S S I N G
%---------------------------------------------------------------------
%Computation of internal forces and moments in elements of Plane Frame
for e = 1 : Model.nElem
    %degrees of freedom corresponding to element e
    dof=[3*FirstNode(e)-2 3*FirstNode(e)-1 3*FirstNode(e) 3*SecondNode(e)-2 3*SecondNode(e)-1 3*SecondNode(e)];

    % Vector of internal forces of element w.r.t. local axes   
    % lforce(:,e) = estif*T*u(dof)-Model.loading.eforce(:,e);   
    lforce(:,e) = ke_local(:,:,e) * T6(:,:,e) * u(dof)-Model.loading.eforce(:,e);
end

Solution.data.lforce = lforce;
Solution.nodalDisplacements = reshape(u,Model.nDof,Model.nNode)';


% Calculating the Reactions
reactionForces = [];
for i = 1 : Model.nNode
    if(Model.bound.resnode(i)==0) %restrained nodes (at support) are identified 
       Rx = 0.0; %components of support reactions along global X and Y directions are initialized
       Ry = 0.0; 
       Mz = 0.0;%support moment is initialized
       for e = 1 : Model.nElem %elements with "node i" as one of the nodes are identified
          %components of force in element e along global X and Y
          %directions are extracted using direction cosines and are
          %then its contribution to Rx and Ry is computed
           if(i == FirstNode(e))  %"node i" is the first node of element e
              Rx = Rx + Model.geometry.Cx(e)*lforce(1,e) - Model.geometry.Cy(e)*lforce(2,e);
              Ry = Ry + Model.geometry.Cy(e)*lforce(1,e) + Model.geometry.Cx(e)*lforce(2,e);
              Mz = Mz + lforce(3,e);
           end
           if(i == SecondNode(e)) %"node i" is the second node of element e
              Rx = Rx + Model.geometry.Cx(e)*lforce(4,e) - Model.geometry.Cy(e)*lforce(5,e);
              Ry = Ry + Model.geometry.Cy(e)*lforce(4,e) + Model.geometry.Cx(e)*lforce(5,e);
              Mz = Mz + lforce(6,e);
           end
       end
       %support reactions are set zero for any unrestrained degree of freedom
       %associated with a support
       if(Model.bound.resx(i) == 1)  
           Rx = 0.0;
       end
       if(Model.bound.resy(i) == 1)
           Ry = 0.0;
       end
       if(Model.bound.resz(i) == 1)
           Mz = 0.0;
       end
       reactionForces = [reactionForces; i,Rx,Ry,Mz];
    end
end
Solution.reactionForces = reactionForces;






%--------------------------------------------------------------------------
% this part will be moved in the future
if isfield(Model,'analysis')
    if Model.analysis.showResults
        printResult(Solution,1)
    end
    if Model.analysis.showDetails
        % print stiffness matrices
        printStiff(Solution,1)
    end
    if Model.analysis.saveToFile
        outFile = fopen(Model.analysis.outputFileName, 'w');
        printSummary(Model,outFile);
        printResult(Solution,outFile);

        fclose(outFile);
    end
end
%--------------------------------------------------------------------------
end



%--------------------------------------------------------------------------

function Model = ModelRevise(Model)

%--------------------------------------------------------------------------
%                 G E N E R A L      I N F O R M A T I O N
%--------------------------------------------------------------------------
nNode = Model.nNode;              
nElem = Model.nElem;             
nDof  = Model.nDof;           

%--------------------------------------------------------------------------
%         R E A D I N G      N O D A L      I N F O R M A T I O N
%--------------------------------------------------------------------------

Xval  = Model.geometry.coordinates(:,1); % x coordinates of nodes 
Yval  = Model.geometry.coordinates(:,2); % y coordinates of nodes 

% restraints of a node
res = ones(nNode,nDof);
for i = 1 : Model.nBoundary
    row = Model.boundary(i,1);
    col = Model.boundary(i,2);
    res(row,col) = 0;
end
resx  = res(:,1); % restraint along x-direction
resy  = res(:,2); % restraint along y-direction
resz  = res(:,3); % restraint about z-direction

%--------------------------------------------------------------------------
%  Generation of  nodal force vector (F) (of size ndof) along global axes
%--------------------------------------------------------------------------

hingeIndex = zeros(nNode,1);
if isfield(Model,'hinge')
    if ~isempty(Model.hinge)
        hingeIndex(Model.hinge)=1;
    end
end
% hinge = 1    regular node = 0

% nodal forces
F = zeros(nNode,nDof);
for i = 1 : Model.nNodalForce
    node = Model.loading.nodalForces(i,1);
    F(node, 1:Model.nDof) = Model.loading.nodalForces(i,2:end);
end
% reshape to vector form
F = reshape(F',1,nDof*nNode)';

%--------------------------------------------------------------------------
%      B O U N D A R Y    C O N D I T I O N    I N F O R M A T I O N
%--------------------------------------------------------------------------
% bn is a counter for the total number of restraints
bn = 0; % initialized to zero
%--------------------------------------------------------------------------
% notation used to indicate restraint associated with a given dof
% 0 indicates RESTRAINED   degree of freedom
% 1 indicates UNRESTRAINED degree of freedom
%--------------------------------------------------------------------------
% A new vector "{resnode}" of size nNode is initialized with ones 
% {resnode} is used to indicate the status of restraint for nodes:
%  - Restrained   nodes (at supports): 0  
%                       (with atleast one RESTRAINED degree of freedom) 
%  - Unrestrained nodes : 1resnode = ones(nNode,1);

resnode = ones(nNode,1); % nodes are initialized as UNRESTRAINED by default

for i = 1 : nNode
   if(resx(i) == 0) % if this dof along x-direction is restrained
        bn = bn+1;
        ifix(bn) = 3*i-2;  % add to the list of restrained dofs
        resnode(i) = 0; 
   end
   if(resy(i) == 0) % if this dof along x-direction is restrained
        bn = bn+1;
        ifix(bn) = 3*i-1;  % add to the list of restrained dofs
        resnode(i) = 0;
   end
   if(resz(i) == 0) % if this dof ABOUT z-direction is restrained
        bn = bn+1;
        ifix(bn) = 3*i;    % add to the list of restrained dofs
        resnode(i) = 0;
   end
end
% NOTE: bn is incremented by 3  for a fixed support 
%       bn is incremented by 2  for a pin support 
%       bn is incremented by 1  for a roller support in the above for-loop
% at the end of the above loop, 
%       bn   is total number of restrained degrees of freedom
%       ifix is list of restrained degrees of freedom for the Plane Frame

%--------------------------------------------------------------------------
%  R E A D I N G     E L E M E N T     I N F O R M A T I O N
%--------------------------------------------------------------------------

Cx      = zeros(nElem,1);    % direction cosines of elements Cx and Cy
Cy      = zeros(nElem,1);
L       = zeros(nElem,1);    % length of each element
eforce  = zeros(6,nElem);    % element force vector with six dofs 
                             % (3 at each node of element)
FirstNode  = Model.geometry.elements(:,1);
SecondNode = Model.geometry.elements(:,2);


sectProp = Model.sections(Model.elemSectId,:);
A        = sectProp(:,3);
E        = sectProp(:,1);
Iz       = sectProp(:,2);


for e = 1 : nElem
    x1 = Xval(FirstNode(e));    
    y1 = Yval(FirstNode(e));  %x and y coordinates of first node
    x2 = Xval(SecondNode(e));   
    y2 = Yval(SecondNode(e)); %x and y coordinates of second node
    le = sqrt((x2-x1)^2+(y2-y1)^2); %length of element e
    %direction cosines of element e
    c  = (x2-x1)/le;
    s  = (y2-y1)/le;
    %Storing values of direction cosines c and s in vectors {Cx} and {Cy} for subsequent use 
    Cx(e) = c;
    Cy(e) = s;
    L (e) = le;
end

distLoad = zeros(nElem,4);
for i = 1 : Model.nTractionForce
    e  = Model.loading.tractionForces(i).element;
    Ty = Model.loading.tractionForces(i).eq{2}; 
    distLoad(e,:) = [Ty(0), 0, Ty(L(e)), L(e)];
end
% reading element load information
% uniformly distributed load of intensity w1 acts from a distance of d1 to d2 
% as measured from left end of element
wr = distLoad(:,1);  % distributed load w1  acting on elements from column 7
d1 = distLoad(:,2);  % left limit d1 of distributed load w1 from columns 8
wl = distLoad(:,3);  % distributed load w2  acting on elements from column 9
d2 = distLoad(:,4); % right limit d2 of distributed load w1 from column 10

%-----------------------
w1 = -wr;
% w1 = wr;

%Computation of equivalent nodal loads owing to uniformly distributed load over some portion of the span
for e = 1 : nElem
    %computation of equivalent nodal loads owing to uniformly distributed load over some portion of the span 
    c1 = d2(e) - d1(e);
    a  = d1(e) + 0.5*c1;
    b  = 0.5*c1 + le - d2(e);
    p1 = (-w1(e)*c1 / (12*le^2)) * (12*a*b^2 + c1^2*(le-3*b));
    p2 = ( w1(e)*c1 / (12*le^2)) * (12*a^2*b+c1^2*(le-3*a));
    p3 = ( w1(e)*c1 * b/le) - (p1+p2)/le;
    p4 = ( w1(e)*c1 * a/le) + (p1+p2)/le;
    eforce(:,e) = eforce(:,e) + [0.0 -p3 p1 0.0 -p4 p2]';
    %Negative of vector of fixed end forces in added to the element force vector
end
%--------------------------------------------------------------------------
%            F O R M I N G     M O D E L     S T R U C T U R E 
%--------------------------------------------------------------------------
% Model.info         = 'Plane Beam/Frame Program';
% Model.analysisType = 'BEAM';
% Model.nDim         = 2;  % 2D (plane beam/frame)
% Model.nNode        = nNode;
% Model.nElem        = nElem;
% Model.nDof         = 3;
% Model.nElemNode    = 2;

% Model.geometry.coordinates = [Xval, Yval];
% Model.geometry.elements    = [FirstNode, SecondNode];
Model.geometry.Cx  = Cx;
Model.geometry.Cy  = Cy;
Model.geometry.L   = L;
Model.geometry.hingIndex = hingeIndex;

Model.properties.A = A;
Model.properties.E = E;
Model.properties.Iz= Iz;

Model.nBoundary    = bn;
Model.bound.resx = resx;
Model.bound.resy = resy;
Model.bound.resz = resz;
Model.bound.resnode = resnode;
Model.bound.ifix = ifix;

Model.loading.F = F;
Model.loading.eforce = eforce;
Model.loading.data = [d1,d2,wl,wr];
end