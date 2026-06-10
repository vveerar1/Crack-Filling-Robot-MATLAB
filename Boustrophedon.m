%% Boustrophedon Path Planner
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function PathEdge = Boustrophedon(Path,splitReg,see,seP,init,wall_fol,known,sim)
%
% Generates a complete-coverage boustrophedon (zig-zag) path for a sequence
% of monotone cells following the traversal order given by Path.
%
% INPUTS:
%   Path     - Cell traversal sequence (from reeb_traversal); each row is [from, to, cellIndex].
%   splitReg - Decomposed cells (polyshape array).
%   see      - Per-cell sweep direction flags.
%   seP      - Start/end coordinates for each Reeb edge.
%   init     - Current robot position; 0 if no prior position.
%   wall_fol - Logical vector; 1 for cells that need a wall-follow pass.
%   known    - Flag: 1 for SCC (known environment), 0 for OnlineSCC.
%   sim      - Flag: 1 for simulation (no boundary clamping).
%
% OUTPUTS:
%   PathEdge - Complete path as sequential (x, y) coordinates.

    global reebEdge reebCell crackEdge crackRaw s allNode spdist a total_length
    o=0;fl_see=false;
    if init~=0;PathEdge=init;else;PathEdge=[];end;sensor=s/sqrt(2);                                          % Square inscribed in sensor circle

    for i = 1:length(Path(:,1))
        EE = Path(i,1:2);
        if ismember(EE,[reebEdge;reebEdge(:,2),reebEdge(:,1)],'rows')
            o=o+1;
            ind=find(ismember(reebEdge,[EE;fliplr(EE)],'rows'));
            if ~isempty(ind);cell=splitReg(reebCell(ind(1)));reebEdge(ind(1),:)=[0,0];end

            orgcell=cell;

                [subXY,flag] = BoustrophedonPath(cell,orgcell,EE,sensor,[],init,wall_fol(i),fl_see);

                if ~sim
                    subXY(subXY(:,1)<sensor,1)=sensor;                  %%% Add the sensor gap
                    subXY(subXY(:,1)>3048-sensor,1)=3048-sensor;
                    subXY(subXY(:,2)<sensor,2)=sensor;                  %%% Add the sensor gap
                    subXY(subXY(:,2)>2898-sensor,2)=2898-sensor;
                end

        else
            if (i~=size(Path,1)) && ismember(Path(i+1,1:2),[reebEdge;reebEdge(:,2),reebEdge(:,1)],'rows') && ~isempty(seP)
                [marker_x1,marker_y1] = addPtsLin([PathEdge(end,1),seP(o+1,1)],[PathEdge(end,2),seP(o+1,2)],s);
                subXY =[PathEdge(end,:);[marker_x1;marker_y1]';seP(o+1,:)];
            else
                if known;[subXY] = [PathEdge(end,:);fliplr(allNode(EE(2),:))];else;[subXY] = PathEdge(end,:);end
            end

            if ~sim
                subXY(subXY(:,1)<sensor,1)=sensor;                  %%% Add the sensor gap
                subXY(subXY(:,1)>3048-sensor,1)=3048-sensor;
                subXY(subXY(:,2)<sensor,2)=sensor;                  %%% Add the sensor gap
                subXY(subXY(:,2)>2898-sensor,2)=2898-sensor;
            end
        end

        if ~isempty(subXY)
            if ~isempty(PathEdge)
                [marker_x1,marker_y1] = addPtsLin([PathEdge(end,1),subXY(1,1)],[PathEdge(end,2),subXY(1,2)],s);
                PathEdge=[PathEdge;[marker_x1;marker_y1]';subXY];init=subXY(end,:);
            else
                PathEdge=[PathEdge;subXY];init=subXY(end,:);
            end
        end
    end
end

function [marker_x,marker_y] = addPtsLin(x,y,marker_dist)
    % Inserts equidistant waypoints along a straight segment.
    dist_from_start = cumsum( [0, sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2)] );marker_x=[];marker_y=[];
    marker_locs = marker_dist : marker_dist : dist_from_start(end);
    if ~isempty(marker_locs) && ~(length(marker_locs)==1 && marker_locs==marker_dist)
        marker_indices = interp1( dist_from_start, 1 : length(dist_from_start), marker_locs);
        marker_base_pos = floor(marker_indices);
        weight_second = marker_indices - marker_base_pos;
        marker_x = [marker_x, x(marker_base_pos) .* (1-weight_second) + x(marker_base_pos+1) .* weight_second];
        marker_y = [marker_y, y(marker_base_pos) .* (1-weight_second) + y(marker_base_pos+1) .* weight_second];
    end
end

function polyout=polyclean(polyin)
    polyout=polyshape();
    if isscalar(polyin)
        poly = regions(polyin);
        poly=poly(fix(poly.area*1e-03)>5);
        polyout = regJoin(poly);
    else
        poly = polyin;
        poly=poly(fix(poly.area*1e-03)>5);
        polyout = poly;
    end
end

function polyout = regJoin(polyin)
    polyout=polyshape();
    for i=1:length(polyin)
        polyout=addboundary(polyout,polyin(i).Vertices);
    end
end

function polyout = regCombine(polyin)
    polyout = polybuffer(union(polybuffer(polyin,1)),-1);
end
