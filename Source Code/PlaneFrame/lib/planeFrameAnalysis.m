function Solution = planeFrameAnalysis(Model)

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
    c = Model.geometry.Cx(e);%direction cosines
    s = Model.geometry.Cy(e);
    le = Model.geometry.L(e);
    ae = Model.properties.A(e)*Model.properties.E(e); %axial rigidity 
    ei = Model.properties.E(e)*Model.properties.Iz(e); %flexural rigidity
    %rotation transformation matrix for element e
    T=[c    s    0    0     0     0;
      -s    c    0    0     0     0;
       0    0    1    0     0     0;
       0    0    0    c     s     0;
       0    0    0   -s     c     0;
       0    0    0    0     0     1;];
    %element stiffness matrix w.r.t local axes
    estif = [ ae/le     0             0          -ae/le     0            0; 
           0       12*ei/le^3    6*ei/le^2       0    -12*ei/le^3   6*ei/le^2;
           0       6*ei/le^2     4*ei/le         0    -6*ei/le^2    2*ei/le ;
         -ae/le      0             0           ae/le      0           0; 
           0      -12*ei/le^3   -6*ei/le^2       0     12*ei/le^3  -6*ei/le^2;
           0       6*ei/le^2     2*ei/le         0    -6*ei/le^2    4*ei/le ];

    lforce(:,e) = estif*T*u(dof)-Model.loading.eforce(:,e);   %Vector of internal forces of element w.r.t. local axes
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
Solution.data.reactionForces = reactionForces;

end