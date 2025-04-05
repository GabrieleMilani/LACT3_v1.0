function [P, polyL, int_i] = Load_wall_data_2(filename, pathname, ax1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE MATRIX OF THE NODES P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%number of blocks and polylines in the walls
[lwpolylayers,lwpolylines] = dxf2coord(filename, pathname);
n_blocks = length(lwpolylayers);
n_polylines = length(lwpolylines);

P =[];

% Define the tolerance for floating-point comparison
tolerance = 1e-3;

% Loop through the original coordinates matrix
index = 1;
for i = 1:size(lwpolylines, 1)
    % Get the current coordinate (consider only the 2nd, 3rd, and 4th columns)
    current_coordinate = lwpolylines(i, 2:4);
    
    % Check if the current coordinate is already in the unique_coordinates matrix
    is_unique = true;
    for j = 1:size(P, 1)
        if all(abs(current_coordinate - P(j, 1:3)) < tolerance)
            is_unique = false;
            break;
        end
    end
    
    % If the coordinate is unique, add it to the matrix P
    if is_unique
        P(index, 1:3) = current_coordinate;
        index = index + 1;
    end
end

% Add the node number in the first coloumn
P = [(1:size(P, 1))', P];
P = round(P,3);

% plot (P(:,2), P(:,3), "o"), hold on % oppure plot (P(:,2),P(:,3),"o")

%initialise the matrix to store the centroids of the blocks and the area

coord_G = zeros(n_blocks,4);
A = zeros(n_blocks,1);
m=1;  %counter to avoid the empty row 
for n = 1 : n_blocks
    temp = [lwpolylines(n == lwpolylines(:,1),2)';...
            lwpolylines(n == lwpolylines(:,1),3)'];
    temp(:,end) = [] ;   % removes last column

    %plot all the elements
    polyin = polyshape(temp(1,:), temp(2,:));

    %create the matrix of the coordinate of the centroids
    coord_G(m,1) = m; 
    [coord_G(m,2), coord_G(m,3)] = centroid(polyin);
    coord_G(m,4) = 0;
        
    %Compute the area
    A(m,1) = area(polyin);
        
    %Write the block numbers in the plot
    m = m+1;
    clear temp;
end


if ~isempty(ax1)
    fig=ancestor(ax1,"figure","toplevel");
    d=uiprogressdlg(fig,'Title','Import model');
    d.Message = 'Get geometry from .dxf file...';
    d.Value = .33;
    d.Message = 'Plot model...';
    axis(ax1,"auto");
    for n = 1 : n_blocks
        i = 1;
        for j = 1 : length(lwpolylines)
            if n == lwpolylines(j,1)
                temp(1,i) = lwpolylines(j,2);
                temp(2,i) = lwpolylines(j,3);
                i = i+1;
            end
        end
        i = i-1;
        temp(:,i) = [] ;   % removes last column
    
        %plot all the elements
        polyin = polyshape(temp(1,:), temp(2,:));
        % figure(1),box on, axis equal, hold on
        % xlabel('x [mm]'),ylabel('y [mm]')
        plot(ax1,polyin), hold(ax1,"on")
        clear temp;
    end
    %plot(coord_G(:,2), coord_G(:,3), 'ok'), hold on
    text(ax1,coord_G(:,2), coord_G(:,3), [num2str(coord_G(:,1))]);
    % plot ruler
    model_aabb=[min(P(:,2)) min(P(:,3)) max(P(:,2)) max(P(:,3))];
    width=10^fix(log10(model_aabb(3)-model_aabb(1)));
    if width<(model_aabb(3)-model_aabb(1))
        width=width*fix((model_aabb(3)-model_aabb(1))/2/width);
    else
        width=width/2;
    end
    height=width/10;
    fill(ax1,[0,width/4,width/4,0], [-height*3,-height*3,-height*2,-height*2], [0 0 0])
    fill(ax1,[width/4,width/2,width/2,width/4], [-height*3,-height*3,-height*2,-height*2], [1 1 1])
    fill(ax1,[width/2,3*width/4,3*width/4,width/2], [-height*3,-height*3,-height*2,-height*2], [0 0 0])
    fill(ax1,[3*width/4,width,width,3*width/4], [-height*3,-height*3,-height*2,-height*2], [1 1 1])
    text(ax1,0,-height*4,"0.00 m");
    text(ax1,width,-height*4,sprintf("%.2f m",width/1000));
    
    % axis(ax1,"padded")
    % - rewrite for adaption of old MATLAB version
    set(ax1,"PlotBoxAspectRatioMode","manual");
    set(ax1,"DataAspectRatio",[1 1 1]);
    xl=get(ax1,"XLim");
    yl=get(ax1,"YLim");
    xm=((xl(2)-xl(1))-(model_aabb(3)-model_aabb(1)))/2;
    ym=((yl(2)-yl(1))-(model_aabb(4)-model_aabb(2)))/2;
    if xm <= ym
        set(ax1,"XLim",[mean(xl)-(model_aabb(3)-model_aabb(1))/2/0.7 mean(xl)+(model_aabb(3)-model_aabb(1))/2/0.7]);
    else
        set(ax1,"YLim",[mean(yl)-(model_aabb(4)-model_aabb(2))/2/0.7 mean(yl)+(model_aabb(4)-model_aabb(2))/2/0.7]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create topology matrix of the polyline polyL(k,:)=[N1 N2 N3 ... x_g y_g z_g A s rho]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Count the maximum number of node within all blocks
n_max_nodi = 0;

for n= 1 : n_blocks
    count = 0;
    for i = 1:length(lwpolylines)
        if lwpolylines(i,1) == n
            count = count+1;
        end
        if count > n_max_nodi 
            n_max_nodi = count;
        end
    end
end
n_max_nodi = n_max_nodi - 1;
polyL = zeros(n_blocks,n_max_nodi+6);

% Cycle within each element in lwpolylines
for element = 1:n_polylines
    % Find the rows that correspond to the current element in lwpolylines
    element_rows = lwpolylines(lwpolylines(:, 1) == element, 2:4);
    
    % Initialization of a index variable for polyL
    polyL_index = 1;
    
    % Cycle within each row of the furrent element
    for i = 1:size(element_rows, 1)
        current_coordinate = element_rows(i, :);
        
        % Compare the current coordinate with the coordinates in P
        for j = 1:size(P, 1)
            if all(abs(current_coordinate - P(j, 2:4)) < tolerance)
                % Check if the index is already present in polyL for the current element
                if ~ismember(j, polyL(element, :))
                    % Add index in the matrix polyL
                    polyL(element, polyL_index) = j;
                    polyL_index = polyL_index + 1;
                end
                break;
            end
        end
    end
end

% Add supplementary columns to polyL

polyL(:, n_max_nodi+1:n_max_nodi+3) = coord_G(:,2:4);
polyL(:, n_max_nodi+4) = A;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IDENTIFY THE INTERFACES int_i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(ax1)
    d.Value = .67;
    d.Message = 'Identify the interfaces...';
end
% Initialize the interface matrix
interfaces = [];

% Get the number of blocks
num_blocks = size(polyL,1);

% Iterate over each pair of blocks to find common edges (interfaces)
for i = 1:num_blocks
    for j = i+1:num_blocks
        % Get the nodes of block i and block j
        nodes_i = nonzeros(polyL(i, 1:n_max_nodi));
        nodes_j = nonzeros(polyL(j, 1:n_max_nodi));
        
        % Find common nodes (edges) between block i and block j
        common_nodes = intersect(nodes_i, nodes_j);
        
        % If two common nodes are found, it means there is an interface
        if length(common_nodes) == 2
            interfaces = [interfaces; common_nodes', i, j];
        end
    end
end

% Calculate direction vectors t_i and n_i
num_interfaces = size(interfaces, 1);
int_i = zeros(num_interfaces, 13);

k = [0, 0, 1]; %versor z

for i = 1:num_interfaces
    % if i==58
    %     a=0
    % end
    N1 = interfaces(i, 1);
    N2 = interfaces(i, 2);
    E1 = interfaces(i, 3);
    E2 = interfaces(i, 4);
    
    % Get coordinates of N1 and N2
    P_N1 = P(N1, 2:4);
    P_N2 = P(N2, 2:4);
    L = norm(P_N2 - P_N1);
    % Calculate t_i
    t_i = (P_N2 - P_N1) / norm(P_N2 - P_N1);
    t_i(3) = 0;
    % Calculate n_i (perpendicular to t_i)         
    n_i = cross(k,t_i);  
    % Calculate d2 (distance between point P_N2 and barycentre of element 2)
    d2=(P_N2-polyL(E2,(size(polyL,2)-5):(size(polyL,2)-3)));
    chek=d2*n_i';
    if chek<0
        E1 = interfaces(i, 4);
        E2 = interfaces(i, 3);
    end
   
    % Store results
    int_i(i, 1:11) = [N1, N2, E1, E2, t_i, n_i, L];
    
end
if ~isempty(ax1)
    % Close dialog box
    d.Value=1;
    close(d)
end
end
