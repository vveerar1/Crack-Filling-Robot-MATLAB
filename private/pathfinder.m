%% Visibility-graph shortest-path planner between two points in a bounded arena
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [waypoint_coordinates,weight] = pathfinder(start_point, end_point, external_boundaries)
% pathfinder  Obstacle-avoiding shortest path via a visibility graph.
%   [waypoint_coordinates, weight] = pathfinder(start_point, end_point, external_boundaries)
%   Finds the shortest path from start_point to end_point within the polygon
%   defined by external_boundaries, routing around boundary vertices as needed.
%
%   OUTPUTS:
%     waypoint_coordinates  Nx2 matrix of (x,y) waypoints from start to end
%     weight                total path length
%
%   INPUTS:
%     start_point          1x2 or 2x1 [x, y] of the start position
%     end_point            1x2 or 2x1 [x, y] of the end position
%     external_boundaries  Mx2 matrix of boundary polygon vertices (x, y per row);
%                          the polygon is treated as closed

%% Main body of function
spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);
spdist2 = @(Ps1,Ps2) sqrt((Ps1(:,1)-Ps2(:,1)).^2 + (Ps1(:,2)-Ps2(:,2)).^2);
total_length = @(Ps) sum(sqrt(sum(diff(Ps).*diff(Ps),2)));

%Initialize empty arrays
initial_combined_nodes = zeros(size(external_boundaries,1)+2,3);

%Create initial node list: start, boundary vertices, end
initial_combined_nodes(1,1:2) = [start_point(1), start_point(2)];
initial_combined_nodes(2:size(external_boundaries,1)+1,1:2) = external_boundaries;
initial_combined_nodes(size(external_boundaries,1)+2,1:2) = [end_point(1), end_point(2)];

%Assign placeholder/tentative distances to all nodes as infinity
initial_combined_nodes(:,3) = Inf*ones(size(initial_combined_nodes,1),1);

%Identify nodes visible from starting point
visible_nodes_ID = zeros(1,size(initial_combined_nodes,1));

%% Build visibility graph: enumerate all node pairs, test line of sight, query shortest path
visible_neighbours_library =  zeros(size(initial_combined_nodes,1),size(initial_combined_nodes,2),size(initial_combined_nodes,1));

visible_index = 0;
ii=nchoosek(1:size(initial_combined_nodes,1),2);ii(:,3)=spdist2(initial_combined_nodes(ii(:,1),1:2),initial_combined_nodes(ii(:,2),1:2));
visibility_h=line_of_sight(initial_combined_nodes(ii(:,1),1:2),initial_combined_nodes(ii(:,2),1:2), external_boundaries);

visible_edges = gather(ii(gather(visibility_h)>0,:));
G=graph(visible_edges(:,1),visible_edges(:,2),visible_edges(:,3));
[path,weight,edge] = shortestpath(G,1,size(initial_combined_nodes,1));
waypoint_coordinates = gather(initial_combined_nodes(path,[1,2]));
end
