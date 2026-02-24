function SIMP_2D_GNTO(DL,DH,nelx,nely,volfrac,rmin)
%--------------------------------------------------------------------------
% % Geometric Nonlinear Topology Optimization base on SIMP
% SIMP-based 2D geometrically nonliear topology optimization
% DL: length (x-direction)
% DH: height (y-direction)
% nelx = number of elements in x
% nely = number of elements in y
% volfrac = volume fraction (constraint)
% rmin = filter radius (in physic units)

p=1;      % penalization (can be increased gradually)
dp=0.1;   % The increment of the penal coefficient
E=3e6;    % Young's modulus
nu=0.4;   % Poisson ratio
xthreshold=0.01; % Low density threshold
force=2000 ;  % Magnitude of external force
IterMax=16;  % Maximum number of Newton-Raphson iterations
pmin=1e-4;   % Minumum density

parent_dir_name ='GNTO_2D_results';
if  exist(parent_dir_name,"dir") 
rmdir(parent_dir_name, 's') % Delete the existed file
end
mkdir(parent_dir_name);

% Calculate nodal coordinates of all elements
% Element size
EL=DL/nelx; EH=DH/nely;
EA=EL*EH; % Element area

% Nodal coordinates (2D grid)
[ey1,ex1] = meshgrid(EH * (0 : nely),EL * (0: nelx));
Nodes=[ex1(:) ey1(:)];   % [x  y]

% -------------------------------------------------------------------------
% Typical 2D cantilever beam boundary condition
% Left side (x = 0) completed fixed
% Load: vertical downward at the middle of the right end (or adjust)
% -------------------------------------------------------------------------

% Nodes on the right end (x = DL)
right_nodes = find(abs(Nodes(:,1) -DL) < 1e-6);
bottom_y = 0;
[~, idx] = min(abs(Nodes(right_nodes,2) - bottom_y));
load_node = right_nodes(idx);

loaddof = 2 * load_node;        % uy dof (2nd dof of that node)
% loaddof = 2*load_node - 1;    % horizontal load 

% Fixed dofs - nodes on left side (x=0)
left_nodes = find(abs(Nodes(:,1)) < 1e-6 );
Fixeddofs = [2*left_nodes-1; 2*left_nodes];      %  ux & uy


% -------------------------------------------------------------------------
% General setup
% -------------------------------------------------------------------------

% Prepare for geometrically nonlinear FEA
nele = nelx*nely;
ndof = 2*(nelx+1)*(nely+1);         % 2 dofs per node
Alldofs=1:ndof;
Freedofs=setdiff(Alldofs,Fixeddofs);

ForceVector = sparse(loaddof,1,-force,ndof,1);   % External force vector

% After setting ForceVector
fprintf(['Load node: %d, coordinates [%.4f, %.4f], ' ...
    'force = %.0f\n'], load_node, Nodes(load_node,:), force);


% Element → node mapping (connectivity)
% Ordering: bottom-left, bottom-right, top-right, top-left

nodegrd = reshape(1:(nelx+1)*(nely+1),nelx+1,nely+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);    % bottom-left node of each element

edofMat = zeros(nele, 8);

e = 0;
for ely = 1:nely
    for elx = 1:nelx
        e = e + 1;
        
        % Node number of bottom-left corner of current element (1-based)
        n_bl = (ely-1)*(nelx + 1) + elx;
        
        n1 = n_bl;               % bottom-left
        n2 = n_bl + 1;           % bottom-right
        n3 = n_bl + (nelx + 1) + 1;   % top-right
        n4 = n_bl + (nelx + 1);       % top-left
        
        % Assign dofs: ux uy for each node
        % Order: n1, n2, n3, n4 (counter-clockwise)
        edofMat(e, :) = [
            2*n1 - 1, 2*n1,     ...  % n1: ux uy
            2*n2 - 1, 2*n2,     ...  % n2
            2*n3 - 1, 2*n3,     ...  % n3
            2*n4 - 1, 2*n4      ...  % n4
        ];
    end
end

% Node numbers per element (for possible visualization or other uses)

Eles= repmat(nodeids(:),1,4)+repmat([0 1 nelx+[2 1]], nele,1); % Nodal number corresponding to each element


% ────────────────────────────────────────────────────────────────
% From here you continue with:
%   - density initialization: x = volfrac * ones(nele,1);
%   - filter setup
%   - main optimization loop
%   - nonlinear FE solver (total Lagrangian + Newton-Raphson)
%   - sensitivity analysis (adjoint or direct)
%   - update densities (e.g. OC or MMA)
% ────────────────────────────────────────────────────────────────

fprintf('2D setup ready. nele = %d, ndof = %d\n', nele, ndof);

% initial density values
%x=volfrac*ones(nele,1);
x = volfrac + 0.01*(rand(nele,1)-0.5);

[H,Hs] =prepareFilter(nele,nelx,nely,EL,rmin);% Prepare for the filter
[Nodes2Ele]=NodeSurroundingElement(Nodes,Eles); % Calculate the center and nodal coordinates of all elements
% Maximum and minimum values for the design variables
xmin    = 0.001*ones(nele,1);
xmax    = ones(nele,1);
% MMA parameters
m=1;                    % one constraint (volume)
nn=nele;                % number of design variables

low = xmin;             % xmin = 1e-3
upp = xmax;             % xmax = 1.0    

c=1000*ones(m,1);       % large → volume constraint almost = equality

a0=1;
[d,a]=deal(zeros(m,1));

%Initial design
xy00=x;                 % make sure column vector, length = nele
[xold1,xold2]=deal(xy00);
Loop=1; change=1.0; maxloop=2; 

% Adjoint vector (for compliance sensitivity)
Lamada=zeros(ndof,1);

% Show the figure for the optimized result
clf; showDensity(Nodes,Eles,xy00(:)); 
while change>1e-5 && Loop<maxloop
    % The update for the penal coefficient
    if(Loop~=0)
        p=p+dp;
        fprintf(['\t\t *** Successful loop converge ***   p= ',num2str(p,3),' \n']);
        if p>2.998
            p=3;
        end
    end
    % =============================================================
    % 1. Visualization
    % =============================================================

    % Show the figure for the optimized result
    clf;
    showDensity(Nodes,Eles,xy00(:));            % 2D visualization function

    
    FileName = sprintf('%s\\Fig_%04d.png', parent_dir_name, Loop);
    saveas(gcf,FileName);

    % Better quality alternative:
    % exportgraphics(gcf, FileName, 'Resolution', 300);

    % =============================================================
    % 2. Optional: identify nodes surrounded by void (if needed)
    % =============================================================
    [uselessNode]=NodeSurroundWithVoidElement(Nodes2Ele, xy00,xthreshold);
    % Geometrically nonlinear FEA + Residual force sensitivity with the element density
%     profile clear
%     profile on
    % =============================================================
    % 3. 2D Finite Element Analysis (linear or geometrically nonlinear)
    % =============================================================

    [U,GKF,dRdpe]=GNFEA_2D(ndof,nele,Freedofs,IterMax,ForceVector, ...
                        Nodes,Eles,edofMat,E,nu,xy00,p, ...
                        uselessNode,pmin,'plane_stress');
    % plot_disp_magnitude(Nodes, Eles, U, ...
    %     'title_str', sprintf('Displacement Magnitude - Iter %d', Loop), ...
    %     'cmap', 'jet', ...
    %     'scale', 1);   % exaggerate 5× if deformation is small

%     profile viewer
    Lamada(Freedofs,:)=GKF(Freedofs,Freedofs)\ForceVector(Freedofs,:);%K*lamada=F

% ────────────────────────────────────────────────────────────────
% Compare displacement U and adjoint Lamada (should be almost equal)
% ────────────────────────────────────────────────────────────────

% Make sure both are column vectors of length ndof
U_full = zeros(ndof, 1);
U_full(Freedofs) = U(Freedofs);   % U is already on Freedofs

Lamada_full = zeros(ndof, 1);
Lamada_full(Freedofs) = Lamada(Freedofs);

% 1. Norm difference
diff_norm = norm(U_full - Lamada_full);
rel_diff  = diff_norm / max(norm(U_full), norm(Lamada_full));

fprintf('\n=== Comparison U vs Lamada (Iter %d) ===\n', Loop);
fprintf('  ||U - Lamada||_2          = %.2e\n', diff_norm);
fprintf('  Relative difference        = %.2e\n', rel_diff);
fprintf(['  Max absolute difference   ' ...
    ' = %.2e\n'], max(abs(U_full - Lamada_full)));
fprintf(['  Mean absolute difference   ' ...
    '= %.2e\n'], mean(abs(U_full - Lamada_full)));

% 2. Check if they are "practically" equal (within solver tolerance)
tol_equal = 1e-8;   % adjust based on your solver tolerance
if rel_diff < tol_equal
    fprintf(['  U and Lamada are practically identical' ...
        ' (within %.1e)\n'], tol_equal);
else
    fprintf(['  U and Lamada differ noticeably ' ...
        '(rel diff = %.2e) — investigate!\n'], rel_diff);
end

% 3. Optional: plot difference field (visual check)
if exist('Nodes','var') && exist('Eles','var')
    figure('Name','U vs Lamada Difference','NumberTitle','off');
    
    % Difference in magnitude
    diff_vec = U_full - Lamada_full;
    diff_ux = diff_vec(1:2:end);
    diff_uy = diff_vec(2:2:end);
    diff_mag = sqrt(diff_ux.^2 + diff_uy.^2);
    
    patch('Vertices', Nodes, 'Faces', Eles, ...
          'FaceVertexCData', diff_mag, ...
          'FaceColor', 'interp', 'EdgeColor', 'none');
    
    colormap('jet'); colorbar;
    title('||U - Lamada|| (magnitude difference)');
    xlabel('x'); ylabel('y');
    axis equal tight;
end

fprintf('======================================\n');

%==========================================================================
    dcdpe=full(-Lamada'*dRdpe)';
    % fprintf('size of dRdpe = [%4d, %4d] \n', size(dRdpe));
    % fprintf('size of dcdpe = %4d \n', size(dcdpe,1));
    % objective function, the sensitivity and the filter
    Comp=ForceVector'*U; % objective
    df0dx= H*(xy00(:).*dcdpe(:))./Hs./max(1e-3,xy00(:));
    df0dx(:)=df0dx(:)./max(abs(df0dx(:)));
    fval=sum(xy00)-(nelx*nely*volfrac);
    dfdx=reshape(EA*ones(1,nelx*nely),1,nelx*nely);
    dfdx=H*(xy00(:).*dfdx(:))./Hs./xy00(:);
    %MMA
    xval=xy00;
%     xold=xy00;
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,nn,Loop,xval,xmin,xmax,xold1,xold2,...
        Comp,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    xold2 = xold1;
    xold1 = xy00;
    xy00 = xmma;
    change=max(max(abs( xy00-xold1)));
    disp([' It.: ' sprintf('%4i\t',Loop) ' Obj.: ' sprintf('%6.3f\t',Comp) ' Vol.: ' ...
        sprintf('%6.4f\t',mean(xy00)) 'ch.:' sprintf('%6.4f\t',change)]); 
    Loop = Loop + 1; 
end
end
%=== Geometrically nonlinear FEA + Residual force sensitivity with the element density===
function [U,GKF,dRdpe] = GNFEA_2D(ndof,nele,Freedofs,IterMax,...
                                ForceVector,Nodes,Elements,edofMat,...
                                E,nu,x,p,uselessNode,pmin, analysis_type)
    % analysis_type = 'plane_stress' or 'plane_strain' (default 'plane_stress')

    if nargin < 14 || isempty(analysis_type)
        analysis_type = 'plane_stress';
    end

    if strcmp(analysis_type, 'plane_strain')
        D0 = LinearElasticD_PlaneStrain(E, nu);
    else
        D0 = LinearElasticD_PlaneStress(E, nu);
    end

    U=zeros(ndof,1);
    loop = 0;           % Iteration number
    du=zeros(ndof,1);   % Displacement increment
    ResidualForceMax=max(abs(ForceVector(Freedofs)));

    TOL=1e-5;           % Newtonian iteration tolerance
    while loop<IterMax  && ResidualForceMax>TOL
        Dimension=2;
        Fint = zeros(ndof,1);% Initialize the internal force vector
        GaussCoordinate=[-1/sqrt(3), 1/sqrt(3)];
        GaussWeight=[1.0D0, 1.0D0];

        % For assembly: 8 dofs/element → max 64 terms per element 
        iK = zeros(64,nele);
        jK = zeros(64,nele);
        sK = zeros(64,nele);

        for i=1:nele
            D=(x(i)^(p)*(1-pmin)+pmin)*D0;

            ElementNodeCoord=Nodes(Elements(i,:),:);                        % 4 x 2
            edof=edofMat(i,:);                                              % 8 dofs
            ElementDisp0 = U(edof);
            ElementDisp  = reshape(ElementDisp0, Dimension, 4);             % 2 x 4
            % Calculate all Gaussian points in each element
            tmp=zeros(8,8);                                                 % element tagent stiffness accumulator 
           
            for LX=1:2
                for LY=1:2
                        xi=GaussCoordinate(LX); eta=GaussCoordinate(LY);
                        [dNdx, detJ] = ShapeFunction2D([xi; eta], ElementNodeCoord);

                        FAC = GaussWeight(LX)*GaussWeight(LY)*detJ;

                        F = ElementDisp*dNdx' + eye(2);                     % 2x2 deformation gradient
                        
                        C = F' * F;                                         % right Cauchy-Green

                        Egreen = 0.5*(C - eye(2));                          % Green-Lagrane strain
                        
                        % Voigt form for plane stress or plane strain: [E11 E22 2*E12]'
                        Strain=[Egreen(1,1); Egreen(2,2); 2*Egreen(1,2)];

                        Stress=D*Strain;        % [S11 S22 S12]' 

                        % B matrix (linear part in total Lagrangian) 3x8
                        
                        BN=zeros(3,8);
                        BG=zeros(4,8);
                        for I=1:4
                            COL=(I-1)*2 + 1:(I-1)*2 + 2;
                            BN(:,COL) = [
                                dNdx(1,I)*F(1,1), dNdx(1,I)*F(2,1);
                                dNdx(2,I)*F(1,2), dNdx(2,I)*F(2,2);
                                dNdx(1,I)*F(1,2) + dNdx(2,I)*F(1,1),...
                                dNdx(1,I)*F(2,2) + dNdx(2,I)*F(2,1)
                               ];
                            % Geometric part for tangent stiffness (SIG in
                            % voigt form)
                            % BG-like for geometric stiffness (4x8 in 2D)
                            BG(:,COL)=[
                                dNdx(1,I) 0;
                                dNdx(2,I) 0;
                                0         dNdx(1,I);
                                0         dNdx(2,I)
                                ];
                        end

                        % Internal force contribution
                        Fint(edof) = Fint(edof) + FAC*BN'*Stress;           % Internal force
                        SIG=[Stress(1) Stress(3); Stress(3) Stress(2)];     % 2x2 stress matrix 

                        % Tangent stiffness contribution
                        EKF = BN'*D*BN + BG'*kron(SIG,eye(2))*BG;             % adjust kron for 2D
                        tmp = tmp+FAC*EKF;
                end
            end
            % Assembly indices
            [dof1, dof2] = meshgrid(edof);
            iK(:,i) = dof1(:);
            jK(:,i) = dof2(:);
            sK(:,i) = tmp(:);
        end

        GKF=sparse(iK,jK,sK,ndof,ndof);% Assemble global stiffness matrix       

        ResidualForce = ForceVector-Fint;
        du(Freedofs) = GKF(Freedofs,Freedofs)\ResidualForce(Freedofs);
        U = U + du;

        if  loop>0
            % Update free dofs if useless nodes exist
            if ~isempty(uselessNode)
                uselessDOF = [2*uselessNode-1; 2*uselessNode];  % ux, uy for each node
                Freedofs = setdiff(Freedofs, uselessDOF);
            end
            ResidualForceMax = max(abs(ResidualForce(Freedofs)));
            fprintf(' N-R iter %4d : Residual = %10.4e\n', loop,...
                full(ResidualForceMax));
        end
        loop = loop + 1;
    end

% =============================================================
% Sensitivity dR/dρe (for adjoint sensitivity in TO)
% =============================================================    
    iK=zeros(8, nele);
    jK=zeros(8, nele);
    sK=zeros(8, nele);

    for i=1:nele
        dDdx=(p*x(i)^(p-1)*(1-pmin))*D0;
    
        ElementNodeCoord = Nodes(Elements(i,:),:);
        edof=edofMat(i,:);
        ElementDisp0 = U(edof);
        ElementDisp  = reshape(ElementDisp0, Dimension, 4);
        tmp=zeros(8,1);
        for LX=1:2
            for LY=1:2
                xi=GaussCoordinate(LX);
                eta=GaussCoordinate(LY);

                [dNdx, detJ]=ShapeFunction2D([xi; eta],ElementNodeCoord);
                FAC = GaussWeight(LX)*GaussWeight(LY)*detJ;

                F = ElementDisp*dNdx' + eye(2);                     % 2x2 deformation gradient                        
                C = F' * F;                                         % right Cauchy-Green
                Egreen = 0.5*(C - eye(2));                          % Green-Lagrane strain;                                
                
                % Voigt form for plane stress or plane strain: [E11 E22 2*E12]'
                Strain=[Egreen(1,1); Egreen(2,2); 2*Egreen(1,2)];

                BN=zeros(3,8);
                for I=1:4
                    COL=(I-1)*2 + 1:(I-1)*2 + 2;
                    BN(:,COL) = [
                        dNdx(1,I)*F(1,1), dNdx(1,I)*F(2,1);
                        dNdx(2,I)*F(1,2), dNdx(2,I)*F(2,2);
                        dNdx(1,I)*F(1,2) + dNdx(2,I)*F(1,1),...
                        dNdx(1,I)*F(2,2) + dNdx(2,I)*F(2,1)
                       ];
                end
                tmp=tmp+ FAC*BN'*(dDdx*Strain);
%                 dRdpe(edof,i) = dRdpe(edof,i) + FAC*BN'*dDdx*Strain;
            end
        end
        iK(:,i)=edof;
        jK(:,i)=i;
        sK(:,i)=tmp;
    end

    dRdpe=sparse(iK,jK,sK,ndof,size(Elements,1));
end
% %===Obtain elements and nodes during deformation large displacement===

% =========================================================================
% Utility function
% =========================================================================

function [Nodes2Ele]=NodeSurroundingElement(Nodes,Eles)
    Nodes2Ele=sparse(size(Nodes,1),size(Eles,1));
    for i=1:size(Nodes,1)
        [hang,~]=find(Eles==i);
        Nodes2Ele(i,hang)=ones(1,size(hang,1));        
    end
end

function [uselessNode]=NodeSurroundWithVoidElement(Nodes2Ele, xPhy,xmin)
    xPhyMatrix=repmat(xPhy',size(Nodes2Ele,1),1);
    a=Nodes2Ele.*xPhyMatrix;
    [max_a,~]=max(a,[],2);
    [uselessNode]=find(max_a<xmin); % Find node numbers surrounded by low densities
end
% Calculate the derivation of shape function and Jacobi matrix

function [dNdx, JacobiDET] = ShapeFunction2D(GaussPoint, ElementNode)
    % GaussPoint = [xi, eta]
    % ElementNode = 4 x 2 (x,y coordinates of 4 nodes)

    % Node ordering: 1: (-1,-1), 2(1,-1), 3:(1,1), 4:(-1,1)
    ParentNodes = [-1  1  1  -1;
                   -1 -1  1   1];

    ParentNDerivative = zeros(2,4);

    for I = 1:4
        XPoint = ParentNodes(1,I);
        YPoint = ParentNodes(2,I);
        ShapePart = [1 + GaussPoint(1)*XPoint, 1 + GaussPoint(2)*YPoint];

        ParentNDerivative(1,I) = 0.25*XPoint*ShapePart(2);
        ParentNDerivative(2,I) = 0.25*YPoint*ShapePart(1); 
    end

    Jacobi = ParentNDerivative * ElementNode;       % 2x2
    JacobiDET = det(Jacobi);

    dNdx = Jacobi \ ParentNDerivative;              % 2x4

end


function [dNdx, JacobiDET] = ShapeFunction3D(GaussPoint, ElementNode)
ParentNodes=[-1  1  1 -1 -1  1  1 -1;
    -1 -1  1  1 -1 -1  1  1;
    -1 -1 -1 -1  1  1  1  1];
ParentNDerivative=zeros(3,8);
for I=1:8
    XPoint = ParentNodes(1,I);
    YPoint = ParentNodes(2,I);
    ZPoint = ParentNodes(3,I);
    ShapePart = [1+GaussPoint(1)*XPoint 1+GaussPoint(2)*YPoint 1+GaussPoint(3)*ZPoint];
    ParentNDerivative(1,I) = 0.125*XPoint*ShapePart(2)*ShapePart(3);
    ParentNDerivative(2,I) = 0.125*YPoint*ShapePart(1)*ShapePart(3);
    ParentNDerivative(3,I) = 0.125*ZPoint*ShapePart(1)*ShapePart(2);
end
Jacobi = ParentNDerivative*ElementNode; % Calculate Jacobi matrix
JacobiDET = det(Jacobi);
% JacobiINV=inv(Jacobi);
% dNdx=JacobiINV*ParentNDerivative;
dNdx=Jacobi\ParentNDerivative;
end

%  Linear elastic material stress-strain matrix
function [D0]=LinearElasticD(E,u) 
D0=(E*(1-u))/((1+u)*(1-2*u))*[1 u/(1-u) u/(1-u) 0 0 0
                             u/(1-u) 1 u/(1-u) 0 0 0
                             u/(1-u) u/(1-u) 1 0 0 0
                             0 0 0 (1-2*u)/(2*(1-u)) 0 0
                             0 0 0 0 (1-2*u)/(2*(1-u)) 0
                             0 0 0 0 0 (1-2*u)/(2*(1-u))];
end


function D = LinearElasticD_PlaneStress(E, nu)
    % Plane stress D matrix (Voigt: εxx, εyy, γxy)
    factor = E / (1 - nu^2);
    D = factor * [1    nu   0;
                  nu   1    0;
                  0    0  (1-nu)/2 ];
end


function D = LinearElasticD_PlaneStrain(E, nu)
    % Plane strain D matrix (Voigt: εxx, εyy, γxy)
    % This is the effective in-plane stiffness when εzz=0
    factor = E * (1 - nu) / ((1 + nu) * (1 - 2*nu));
    D = factor * [1          nu/(1-nu)    0;
                  nu/(1-nu)  1            0;
                  0          0          (1-2*nu)/(2*(1-nu)) ];
end

function plot_disp_magnitude(Nodes, Elements, U, varargin)
    % Optional inputs:
    %   'title_str'   - string for title
    %   'cmap'        - colormap name, e.g. 'jet', 'parula', 'hot'
    %   'scale'       - exaggerate displacement (default 1)

    p = inputParser;
    addParameter(p, 'title_str', 'Displacement Magnitude', @ischar);
    addParameter(p, 'cmap', 'parula', @ischar);
    addParameter(p, 'scale', 1, @isscalar);
    parse(p, varargin{:});

    % Extract ux, uy
    ux = U(1:2:end);
    uy = U(2:2:end);
    disp_mag = sqrt(ux.^2 + uy.^2);

    % Deformed nodes (optional exaggeration)
    Nodes_def = Nodes + p.Results.scale * [ux uy];

    figure('Name','Displacement Magnitude','NumberTitle','off');
    patch('Vertices', Nodes_def, ...
          'Faces', Elements, ...
          'FaceVertexCData', disp_mag, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    colormap(p.Results.cmap);
    colorbar('FontSize', 12);
    title(p.Results.title_str, 'FontSize', 14);
    xlabel('x'); ylabel('y');
    axis equal tight;
    grid on;
end