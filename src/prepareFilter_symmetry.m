% Prepare Filter
function [H,Hs] =prepareFilter_symmetry(Nele,Nelx,Nely,El,rmin)
rmin=rmin/El;
iH = ones(Nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;

for i1 = 1:Nelx
    for j1 = 1:Nely
        e1 = (i1-1)*Nely+j1;        % 2D element index (column-major)

        % Search window in 2D
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),Nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),Nely)
                e2 = (i2-1)*Nely+j2;
                dist = sqrt( (i1-i2)^2 + (j1-j2)^2 );   % only x and y

                if dist <= rmin % important: keeps it circular, not square
                   k = k+1;
                   iH(k) = e1;
                   jH(k) = e2;
                   sH(k) = max(0,rmin-dist);
                end
            end
        end
    end
end

% Trim to actual number of entries
iH = iH(1:k);
jH = jH(1:k);
sH = sH(1:k);

H = sparse(iH,jH,sH);
Hs = sum(H,2);
end