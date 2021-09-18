function Solution = trussAnalysisExplicit(Model)
% PURPOSE: Analysing trusses with explicit stiffness matrices    
%
% INPUT(S): 
%   - Model: a structure including the model data
%
% OUTPUT(S):
%   - Solution: a structure including the solution of the analysis

% storing information in Solution structs
time = clock;
Solution.info = sprintf(['Solution obtained on %s %d:%d:%d',...
                         '\nExplicit Stiffness Matrices' ],...
                         date,time(4),time(5),floor(time(6)));
% Pre-processing
%--------------------------------------------------------------------------
%                        I N P U T     D A T A
%--------------------------------------------------------------------------

% control variables

nNode   = Model.nNode;
nElem   = Model.nElem;
nDof    = Model.nDof;

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
    lambda = [c*c c*s; c*s s*s];
    T = [lambda -lambda; -lambda lambda];
    
    k(:,:,i) = E(i)*A(i)/L(i)   *    T;
    
    % TODO : this only works for 2D
    T_sigma(:,:,i) = [c s 0 0; 0 0 c s];
    
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