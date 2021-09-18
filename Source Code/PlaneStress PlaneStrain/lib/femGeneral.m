function Solution = femGeneral(Model)
% PURPOSE: finite element analysis of a given model
%
% INPUT(S): 
%   - Model: a structure including the following fields

% DEVELOPMENT HISTORY:
%   [2018, May, 25] Shahrokh Shahi -- transformation of the basic one
%   [2018, XXX, XX] Shahrokh Shahi -- 
%   [2018, XXX, XX] Shahrokh Shahi -- 

    %% 2D Plane Stress/Strain  and General 3D case
    
    nprops            = Model.nSection;
    materialprops     = Model.sections;               % NOTE:
    % in gui the first material property is "E"; however, in my
    % elasticity-based FEM code, the first material property is "mu" (shear
    % modulus), so, a transformation is inevitable here:
    % G or mu = E / 2(1+nu)
    materialprops(:,1) = materialprops(:,1) ./ ...
                        (2 .* (1 + materialprops(:,2)));
    
    ncoord            = Model.nDim;
    ndof              = Model.nDof;
    nnode             = Model.nNode;
    coords            = Model.geometry.coordinates';  % trn
    nelem             = Model.nElem;
    maxnodes          = double(Model.nElemNode);
    connect           = Model.geometry.elements';     % trn
    nelnodes          = ones(nelem, 1) * maxnodes; 
    elident           = Model.elemSectId;% use this reserved var. as a key
    nfix              = Model.nBoundary;
    fixnodes          = Model.boundary';               %trn
    ndload            = Model.nTractionForce;
    
%   dloads            = Model.loading.tractionForces';  %trn
    traction            = Model.loading.tractionForces; 

    Solution.info = ['Solution obtained on ', date];
    
    
    % fem linear procedure:
    dofs = zeros(ndof*nnode,1); 
    
    K = globalstiffness(ncoord,ndof,nnode,coords,nelem,maxnodes, ...
                        elident,nelnodes,connect,materialprops,dofs);
    
    % traction forces vector
    r = globaltraction(ncoord,ndof,nnode,ndload,coords,nelnodes, ...
                       elident,connect,traction,dofs);                
    
    R = r;
    
    % imposing displacement boundary conditions
    debc = (fixnodes(1,:) - 1) * ndof + fixnodes(2,:);
    ebcVals = fixnodes(3,:)';
    
    [dofs, rf] = solution(K, R, debc, ebcVals);
    
    
    if isfield(Model,'analysis')
        if Model.analysis.saveToFile
            fileSaveId = fopen(Model.analysis.outputFileName, 'w');
            printSummary(Model,fileSaveId);
            print_results(fileSaveId, ...
                          nprops,materialprops,ncoord,ndof,nnode,coords, ...
                          nelem,maxnodes,connect,nelnodes,elident, ...
                          nfix,fixnodes,ndload,traction,dofs);
            fclose(fileSaveId);
        end
        
        if Model.analysis.showDetails
            % print stiffness matrices
        end
    end
    
    
    % temporary displyin on command window (TODO: improve the output format)
    print_results(1, ...
                  nprops,materialprops,ncoord,ndof,nnode,coords, ...
                  nelem,maxnodes,connect,nelnodes,elident, ...
                  nfix,fixnodes,ndload,traction,dofs);

    Solution.nodalDisplacements = reshape(dofs,ncoord,nnode)';
end
















