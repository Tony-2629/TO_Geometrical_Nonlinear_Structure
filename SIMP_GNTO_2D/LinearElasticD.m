%  Linear elastic material stress-strain matrix
function [D0]=LinearElasticD(E,u) 
D0=(E*(1-u))/((1+u)*(1-2*u))*[1 u/(1-u) u/(1-u) 0 0 0
                             u/(1-u) 1 u/(1-u) 0 0 0
                             u/(1-u) u/(1-u) 1 0 0 0
                             0 0 0 (1-2*u)/(2*(1-u)) 0 0
                             0 0 0 0 (1-2*u)/(2*(1-u)) 0
                             0 0 0 0 0 (1-2*u)/(2*(1-u))];
end