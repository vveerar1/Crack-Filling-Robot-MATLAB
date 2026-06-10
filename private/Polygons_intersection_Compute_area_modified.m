function [S S_area] = Polygons_intersection_Compute_area_modified(S)
% Modified by: Vishnu Veeraraghavan, Automated Control Systems and Robotics Lab (vveerar1@binghamton.edu).
% Guillaume JACQUENOT
% guillaume at jacquenot at gmail dot com
% 2007_10_08
% 2009_06_16
% Compute area of each polygon of in S.
% Results are stored as a field in S

S_area = struct('A', {});
for i=1:numel(S)
    S(i).area = 0;
    S_area(i).A = zeros(1,numel(S(i).P));
    for j=1:numel(S(i).P)
        S_area(i).A(j) = S(i).P(j).area;
        S(i).area      = S(i).area + S(i).P(j).area;        
    end
end

