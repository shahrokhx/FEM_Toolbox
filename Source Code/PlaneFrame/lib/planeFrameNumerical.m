function Solution = planeFrameNumerical(Model)

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
    
    coord = [Model.geometry.coordinates(FirstNode(e),:) ; ...
             Model.geometry.coordinates(SecondNode(e),:)]';
    % element stiffness matrix w.r.t local axes
    estif = elstif(Model.nDof,coord,...
                   Model.properties.E(e), ...
                   Model.properties.A(e), ...
                   Model.properties.Iz(e));
               
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
     T6(:,:,e) = T;
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
K(ifix,:) = 0; %Elements in rows with d.o.f corresponding to ifix are made 0                 
K(:,ifix) = 0; %Elements in columns with d.o.f corresponding to ifix are made 0                    
K(ifix,ifix) = eye(length(ifix));%Diagonal elements with d.o.f corresponding to ifix are made 1     
Fc(ifix) = 0; %corresponding columns of Global force vector are also made zero
%---------------------------------------------------------------------
%obtaining vector of nodal displacements along global axes
u = inv(K) * Fc;

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


function k = elstif(nDof,coordinates,E,A,I)

    xiList = [-0.5773502692,  0.5773502692];
    wList  = [1.0          ,  1.0         ];
    nElemNode = 2;  %hard-coded for truss/beams
    nDim      = size(coordinates,1);
    
    dCoord = diff(coordinates');
    L      = norm(dCoord);
    coordinates = [0 L]; % in 1D view

    k2 = zeros(2);
    k4 = zeros(4);
    for intpt = 1 : length(xiList)

      xi      = xiList(intpt);  % local corrdinate
      w       = wList (intpt);
      % shape functions and derivatives at local integration points
      N       = shapefunctions     (nDim-1,nElemNode,xi);
      dNdxi   = shapefunctionderivs(nDim-1,nElemNode,xi);
      
      d2Hdxi2 = [(3/2)*xi, (-1 + 3*xi)/2 * (L/2), ...
                -(3/2)*xi, ( 1 + 3*xi)/2 * (L/2)];
      d2Hdxi2 = d2Hdxi2' * d2Hdxi2;

      
      % Jacobian matrix
      J      = coordinates * dNdxi;

      % inverse and determinent of the Jacobian matrix
      invJ   = inv(J);
      detJ   = abs(det(J));

      % derivatives of the shape functions in the global coordinate 
      dNdx   = dNdxi * invJ;

      % B matrix
      B      = dNdx';

      % numerical summation
      k2     = k2     +        B'*E*B  *  detJ * A * w;
      k4     = k4     +        (8*E*I/L^3) * d2Hdxi2  ;
      
    end
    
    k = zeros(nElemNode * nDof);
    k([1,4]    ,[1,4]    ) = k2;
    k([2:3,5:6],[2:3,5:6]) = k4;    
end



%--------------------------------------------------------------------------



function Model = ModelRevise(Model)

%--------------------------------------------------------------------------
%                 G E N E R A L      I N F O R M A T I O N
%--------------------------------------------------------------------------
nNode = Model.nNode;              
nElem = Model.nElem;             
nDof  = Model.nDof;           
nDim  = Model.nDim;
nElemNode = Model.nElemNode;
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
%--------------------------------------------------------------------------
if Model.nTractionForce > 0
    NPOINTS = 3;
    traction    = Model.loading.tractionForces;
    
    for iLoad = 1 : length(traction)
        element = traction(iLoad).element;
        face    = traction(iLoad).face;
        loadEq  = traction(iLoad).eq;

        nodeIdList = facenodes(nDim,nElemNode,face);
        nElemReduc = length(nodeIdList);
        
        elemId     = Model.geometry.elements(element,nodeIdList);
        faceCoord  = Model.geometry.coordinates(elemId,:)';   % here (in 2D) it is an "edge"
        Le         = norm(diff(faceCoord'));


        [xiList, wList] = intQuad(NPOINTS);

        fel = zeros(nElemNode*nDof , 1);
        for intpt = 1 : length(xiList)
            xi    = xiList(:,intpt)';
            w     = wList (intpt)   ;

            % shape functions and derivatives at local integration points
            [N,dNdxi] = beamShapeFunctions(xi);

            N(2) = (Le/2)*N(2); 
            N(4) = (Le/2)*N(4);

            % Jacobian matrix
            detJ = Le/2;

            load = zeros(length(loadEq),1);
            % coordinate transform to calculate loads at integration point
            x = (1+xi)/2 * Le;
            for i = 1 : length(load)
              load(i,1) = loadEq{i}(x);
            end
            load;
            
            fel = fel + ...
            [N(1) * eye(2); ...
             N(2) , N(2)  ; ...
             N(3) * eye(2); ...
             N(4) , N(4)  ] * load      *     detJ * w ;
        end
        % fel = reshape(fel,nDof,length(nodeIdList))'
        % f(nodeIdList,:)=fel
        eforce(:,element) = eforce(:,element) + fel;
    end

end
% distLoad = zeros(nElem,4);
% for i = 1 : Model.nTractionForce
%     e  = Model.loading.tractionForces(i).element;
%     Ty = Model.loading.tractionForces(i).eq{2}; 
%     distLoad(e,:) = [Ty(0), 0, Ty(L(e)), L(e)];
% end
% 
% wr = distLoad(:,1);  % distributed load w1  acting on elements from column 7
% d1 = distLoad(:,2);  % left limit d1 of distributed load w1 from columns 8
% wl = distLoad(:,3);  % distributed load w2  acting on elements from column 9
% d2 = distLoad(:,4); % right limit d2 of distributed load w1 from column 10
% 
% 
% w1 = -wr;
% %Computation of equivalent nodal loads owing to uniformly distributed load over some portion of the span
% for e = 1 : nElem
%     %computation of equivalent nodal loads owing to uniformly distributed load over some portion of the span 
%     c1 = d2(e) - d1(e);
%     a  = d1(e) + 0.5*c1;
%     b  = 0.5*c1 + le - d2(e);
%     p1 = (-w1(e)*c1 / (12*le^2)) * (12*a*b^2 + c1^2*(le-3*b));
%     p2 = ( w1(e)*c1 / (12*le^2)) * (12*a^2*b+c1^2*(le-3*a));
%     p3 = ( w1(e)*c1 * b/le) - (p1+p2)/le;
%     p4 = ( w1(e)*c1 * a/le) + (p1+p2)/le;
%     eforce(:,e) = eforce(:,e) + [0.0 -p3 p1 0.0 -p4 p2]';
%     %Negative of vector of fixed end forces in added to the element force vector
% end

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%            F O R M I N G     M O D E L     S T R U C T U R E 
%--------------------------------------------------------------------------
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
% Model.loading.data = [d1,d2,wl,wr];
end

%--------------------------------------------------------------------------
% Gauss-Quadrature Numerical Integration Points
function [xi, w] = intQuad(npoints)
    % xi: integration points (xi)
    % w : associated weights
    switch npoints
        case 1
            xi = 0.0;
            w  = 2.0;
        case 2
            xi = [-1/sqrt(3), 1/sqrt(3)];
            w  = [1.0 ,1.0];
        case 3
            xi = [-sqrt(3/5), 0.0, sqrt(3/5)];
            w  = [5/9, 8/9, 5/9];
        case 4
            xi1= sqrt(3/7 - (2/7)*sqrt(6/5));
            xi2= sqrt(3/7 + (2/7)*sqrt(6/5));
            xi = [-xi2, -xi1, xi1, xi2];
            w1 = (18 + sqrt(30)) / 36; 
            w2 = (18 - sqrt(30)) / 36;
            w  = [w2, w1, w1, w2];
        case 5
            xi1 = (1/3)*sqrt(5-2*sqrt(10/7));
            xi2 = (1/3)*sqrt(5+2*sqrt(10/7));
            xi  = [-xi2, -xi1, 0, xi1, xi2];
            w1  = (322 + 13*sqrt(70))/900;
            w2  = (322 - 13*sqrt(70))/900;
            w   = [w2, w1, 128/225, w1, w2];
    end
end