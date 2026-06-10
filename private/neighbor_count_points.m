% neighbor_count_points  Classify binary-skeleton pixels by 8-neighbour count.
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [BW_out,Points] = neighbor_count_points(BW,operation)
% neighbor_count_points  Classify binary-skeleton pixels by 8-neighbour count.
%   [BW_out,Points] = neighbor_count_points(BW,operation)
%   operation 0: isolated pixels (0 neighbours), 1: endpoints (1 neighbour),
%   2: branch/intersection points (>2 neighbours).
%   Returns BW_out (binary mask) and Points ([row col] list) for the selected class.

% Operation:
% 0: Single pixels
% 1: End Points
% 2: Intersection Points

[ll,kk]=find(BW);
if operation < 2
    test=find(sum(cell2mat(arrayfun(@(row,col) [BW(row-1,col-1:col+1) ...
                    BW(row, col-1) BW(row, col+1)...
                    BW(row+1,col-1:col+1)],ll,kk,'UniformOutput',false)),2)==operation);
else
    test=find(sum(cell2mat(arrayfun(@(row,col) [BW(row-1,col-1:col+1) ...
                    BW(row, col-1) BW(row, col+1)...
                    BW(row+1,col-1:col+1)],ll,kk,'UniformOutput',false)),2)>operation);
end

Points=[ll(test) kk(test)];

BW_out=full(sparse([ll(test);size(BW,1)],[kk(test);size(BW,2)],1));BW_out(end,end)=0;
