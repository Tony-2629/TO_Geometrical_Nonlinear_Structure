function showDensity(Nodes, Elements, rho)
    % Show the optimized density distribution in 2D (quad elements)
    % Nodes     : [nNodes x 2] node coordinates (x,y)
    % Elements  : [nElem x 4] element connectivity (node numbers)
    % rho       : [nElem x 1] density value (0-1)

    % Threshold to decide which elements are "visible" (common: 0.1 or 0.3)
    threshold = 0.1;
    % Count how many elements have significant density
    highDensityElements = find(rho > threshold);
    nVisible = length(highDensityElements);

    if nVisible == 0
        warning ('No elements above density threshold %2f', threshold);
        return;
    end

    % Preallocate faces (one face = the quad itself) and colors
    AllFaces = zeros(nVisible, 4);
    FaceColor = zeros(nVisible, 3);
    % Fill faces and grayscale colors
    for k = 1:nVisible
        e = highDensityElements(k);         % element index
        AllFaces(k,:) = Elements(e,:);      % [n1 n2 n3 n4]

        % Grayscale: while(1) = voide → black (0) = solid
        % Most common mapping in 2D topo opt papers
        gray = 0.2 + 0.8*(1 - rho(e));      % rho=1 → ~0.2 (dark), rho=0 → 1.0 (white)
        % Alternative popular choice: gray = 1 - rho(e); rho=1 black, rho=0 white
        
        FaceColor(k,:) = [gray, gray, gray]; % RGB = same value → grayscale
    end
    
    % Plot using patch
    patch('Vertices', Nodes, ...
          'Faces', AllFaces,...
          'FaceVertexCData', FaceColor,...
          'FaceColor', 'flat', ...
          'EdgeColor', [0.7 0.7 0.7], ...
          'LineWidth', 0.3, ...
          'FaceAlpha', 1);
    % View setting for 2D
    view(2);                % force 2D view
    axis equal;
    axis tight;

    hold on

    % Optional: add a very light mesh outline (helps see structure)
    % patch('Vertices', Nodes, 'Faces', Elements, ...
    %       'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.3);

    drawnow;   % make sure it updates immediately    

end
    
   