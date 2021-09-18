function Solution = trussAnalysisNumerical(Model)
% PURPOSE: Analysing trusses with numerical integration stiffness matrices    
%
% INPUT(S): 
%   - Model: a structure including the model data
%
% OUTPUT(S):
%   - Solution: a structure including the solution of the analysis

% DEVELOPMENT HISTORY:
%   [2018, Oct, 23] Shahrokh Shahi -- development (function mode)
%   [2018, Nov, 13] Shahrokh Shahi -- Handlig traction (distributed axial) 

time = clock;
Solution.info = sprintf(['Solution obtained on %s %d:%d:%d',...
                         '\nNumerical Integration Stiffness Matrices' ],...
                         date,time(4),time(5),floor(time(6)));
% Pre-processing
%--------------------------------------------------------------------------
%                        I N P U T     D A T A
%--------------------------------------------------------------------------

% control variables

nNode     = Model.nNode;
nElem     = Model.nElem;
nDof      = Model.nDof;
nElemNode = Model.nElemNode;
nDim      = Model.nDim;

%nodal coordinates
coordinates = Model.geometry.coordinates;

% nodal forces
F = zeros(size(coordinates));
for i = 1 : Model.nNodalForce
    node = Model.loading.nodalForces(i,1);
    F(node, 1:Model.nDim) = Model.loading.nodalForces(i,2:end);
end

% reshape to vector form
Feq=reshape(F',1,nDof*nNode)';  
 
% elements connectivity
elements = Model.geometry.elements;

% section/materaial properties
A = zeros(1, Model.nElem);
E = zeros(1, Model.nElem);
for i = 1 : Model.nElem
    E(i) = Model.sections(Model.elemSectId(i), 1);
    A(i) = Model.sections(Model.elemSectId(i), 2);
end

% boundary condition
boundary = Model.boundary';


%--------------------------------------------------------------------------
%                         P R O C E S S I N G
%--------------------------------------------------------------------------
% initializing the loop variables
K = zeros(nDof*nNode);
k(:,:,nElem) = zeros(2 * nDof);
L = zeros(1,nElem);

% loop over elements
for i = 1 : nElem

    a = coordinates(elements(i,1),1:2); % start node coordinates
    b = coordinates(elements(i,2),1:2); % end node coordinates
   
    % calculating length
    L(i) = sqrt((a(1)-b(1)).^2 + (a(2)-b(2)).^2);  
    
    % sin and cos
    c = (b(1) - a(1))/L(i);  % cos(t) = (x2 - x1)/L
    s = (b(2) - a(2))/L(i);  % sin(t) = (y2 - y1)/L
    
    % TODO: this only works for 2D
    T = [[c s; -s c] zeros(2)    ; 
         zeros(2)    [c s; -s c]];
    T_sigma(:,:,i) = [c s 0 0; 0 0 c s];
    T_rot  (:,:,i) = T;
    
    k(:,:,i) = T' * elstif(nDof,[a;b]',E(i),A(i)) * T;
    
    % node i --------> dof ndof*i - ndof + 1 : ndof*i
    % e.g.  node 1 --> dof 1 : 2
    %       node 3 --> dof 5 : 6
    inode = elements(i,1);
    jnode = elements(i,2);
    
    index = [nDof*(inode-1)+1:nDof*inode,nDof*(jnode-1)+1:nDof*jnode];
    K(index,index) = K(index,index) + k(:,:,i);
end
Solution.data.elementStiff = k;
Solution.data.globalStiff = K;

% handling distributed axial force
F_dist = zeros(size(coordinates));
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
        faceCoord  = Model.geometry.coordinates(elemId,:)';
        L          = norm(diff(faceCoord'));
        
        [xiList, wList] = intQuad(NPOINTS);
        fel = zeros(nElemNode*nDof , 1);
        for intpt = 1 : length(xiList)
            xi    = xiList(:,intpt)';
            w     = wList (intpt)   ;

            % shape functions and derivatives at local integration points
            N     = shapefunctions     (nDim-1,nElemReduc,xi);
            dNdxi = shapefunctionderivs(nDim-1,nElemReduc,xi);

            % Jacobian matrix
            J     = faceCoord * dNdxi;

            % inverse and determinent of the Jacobian matrix
            if size(J) == [2,1]
              detJ = norm(J);
            else
              detJ  = abs(det(J));
            end
            detJ;     % detJ=L/2 

            load = zeros(length(loadEq),1);
            
            % coordinate transform to calculate loads at integration point
            x = (1+xi)/2 * L;
            for i = 1 : length(load)
              load(i,1) = loadEq{i}(x);
            end
            load;

            fel = fel + ...
            [N(1) * eye(2); ...
             N(2) * eye(2)] * load      *     detJ * w ;
        end
        fel_global = T_rot(:,:,element)' * fel;
        fel_global = reshape(fel_global,nDof,length(nodeIdList))';
        F_dist(elemId,:) = F_dist(elemId,:) + fel_global;
    end
end
Feq_dist = reshape(F_dist',1,nDof*nNode)';
Feq = Feq + Feq_dist;
Solution.data.globalForce = Feq;

% imposing displacement boundary conditions
debc = double(nDof) .* (boundary(1,:) - 1) + boundary(2,:);
ebcVals = boundary(3,:)';
[dofs, rf] = solution(K, Feq, debc, ebcVals);
Solution.nodalDisplacements = reshape(dofs,Model.nDim,Model.nNode)';

Solution.reactionForces = zeros(Model.nNode, Model.nDim);
for i = 1 : size(boundary,2)
    node = boundary(1,i);
    dof = boundary(2,i);
    Solution.reactionForces(node,dof) = rf(i);
end

%--------------------------------------------------------------------------
%                      P O S T - P R O C E S S I N G
%--------------------------------------------------------------------------    
% calculating element stresses

Solution.elementStresses.sigma11 = zeros(nElem,1);
Solution.internalForces.axial = zeros(nElem,1);

for i = 1 : nElem    
    
    inode = elements(i,1);
    jnode = elements(i,2);
    
    index = [nDof*(inode-1)+1:nDof*inode,nDof*(jnode-1)+1:nDof*jnode];
    
    sigma = E(i) * [-1/L(i), 1/L(i)] * T_sigma(:,:,i) * dofs(index);
    
    Solution.elementStresses.sigma11(i) = sigma;
    Solution.internalForces.axial(i,1) = sigma * A(i);
end
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
function [d, rf] = solution(K, R, debc, ebcVals)
    dof = length(R);
    df = setdiff(1:dof, debc);
    Kf = K(df, df);
    Rf = R(df) - K(df, debc)*ebcVals;
    dfVals = Kf\Rf;
    d = zeros(dof,1);
    d(debc) = ebcVals;
    d(df) = dfVals;
    rf = K(debc,:)*d - R(debc);
end

%--------------------------------------------------------------------------
function k = elstif(nDof,coordinates,E,A)

    xiList = [-0.5773502692,  0.5773502692];
    wList  = [1.0          ,  1.0         ];
    nElemNode = 2;  %hard-coded for trusses
    nDim      = size(coordinates,1);
    
    dCoord = diff(coordinates');
    L      = norm(dCoord);
    coordinates = [0 L]; % in 1D view

    k2 = zeros(nDof);
    
    for intpt = 1 : length(xiList)

      xi    = xiList(intpt);  % local corrdinate
      w     = wList (intpt);
      % shape functions and derivatives at local integration points
      N     = shapefunctions     (nDim-1,nElemNode,xi);
      dNdxi = shapefunctionderivs(nDim-1,nElemNode,xi);

      % Jacobian matrix
      J     = coordinates * dNdxi;

      % inverse and determinent of the Jacobian matrix
      invJ  = inv(J);
      detJ  = abs(det(J));

      % derivatives of the shape functions in the global coordinate 
      dNdx  = dNdxi * invJ;

      % B matrix
      B     = dNdx';

      % numerical summation
      k2     = k2     +        B'*E*B  *  detJ * A * w;

    end
    k = zeros(nElemNode * nDof);
    k([1,3],[1,3]) = k2;
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
%--------------------------------------------------------------------------