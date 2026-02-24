%  Linear elastic material stress-strain matrix
function [D0]=LinearElasticD(E,u) 
D0=(E*(1-u))/((1+u)*(1-2*u))*[1 u/(1-u) u/(1-u) 0 0 0
                             u/(1-u) 1 u/(1-u) 0 0 0
                             u/(1-u) u/(1-u) 1 0 0 0
                             0 0 0 (1-2*u)/(2*(1-u)) 0 0
                             0 0 0 0 (1-2*u)/(2*(1-u)) 0
                             0 0 0 0 0 (1-2*u)/(2*(1-u))];
end

%
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