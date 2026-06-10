%% Polygon-based segment visibility check
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [visibility,in,out] = line_of_sight_poly(observer_node, target_node, external_boundaries)
% line_of_sight_poly  Polygon-based segment visibility check.
%   [visibility,in,out] = line_of_sight_poly(observer_node, target_node, external_boundaries)
%   A segment is visible (visibility=1) iff it stays entirely inside the
%   boundary region; also returns the clipped inside (in) and outside (out) parts.

on = observer_node;
tn = target_node;
ed = polyshape(external_boundaries(1,:),external_boundaries(2,:));

visibility = ones(size(on,1),1);
for i = 1: size(on,1)
        [in,out]=intersect(ed,[on(i,:);tn(i,:)]);
        if ~isempty(out)
            visibility(i)=0;    % The Target node is not visible from the observing node
        end
end

end
