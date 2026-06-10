%% Boustrophedon Path Planner with Cell Connection
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

% Optional args crack_reach (robot reach radius for SCC) and n_crack (number of crack-graph nodes);
% OnlineSCC omits both, which defaults wall clamping to the sensor inset.
function PathEdge = Boustrophedon_CellCon(Path,splitReg,see,seP,init,wall_fol,known,sim,crack_reach,n_crack)
%
% Generates a complete-coverage boustrophedon path with global cell-connection
% optimisation. When the decomposed workspace contains disjoint cells, a
% shortest-path digraph selects the lowest-cost start/end orientation for
% each cell before the zig-zag sweeps are concatenated.
%
% INPUTS:
%   Path        - Cell traversal sequence (from reeb_traversal).
%   splitReg    - Decomposed cells (polyshape array).
%   see         - Per-cell sweep direction flags.
%   seP         - Start/end coordinates for each Reeb edge.
%   init        - Current robot position; 0 if no prior position.
%   wall_fol    - Logical vector; 1 for cells needing a wall-follow pass.
%   known       - Flag: 1 for SCC (known environment), 0 for OnlineSCC.
%   sim         - Flag: 1 for simulation (no boundary clamping).
%   crack_reach - (optional) Robot reach radius for near-wall crack fills (SCC only).
%   n_crack     - (optional) Number of crack-graph nodes; edges to crack nodes use crack_reach.
%
% OUTPUTS:
%   PathEdge    - Complete path as sequential (x, y) coordinates.

    global reebEdge reebCell crackEdge crackRaw s allNode spdist a total_length spdist2
    o=0;fl_see=false;
    if init~=0;PathEdge=init;else;PathEdge=[];end;sensor=s/sqrt(2);                                          % Square inscribed in sensor circle
    if nargin<9;  crack_reach=[]; end                                                                         % SCC passes robot reach radius; OnlineSCC omits -> sensor inset used everywhere
    if nargin<10; n_crack=0;      end                                                                         % Number of crack-graph nodes; only edges to a crack node use crack_reach

    %%% Cell Connection
    reebT=reebEdge;ccNode=[];ccW=[];
    for i = 1:length(Path(:,1))
        EE = Path(i,1:2);
        if ismember(EE,[reebT;reebT(:,2),reebT(:,1)],'rows')
            ind=find(ismember(reebT,[EE;fliplr(EE)],'rows'));
            if ~isempty(ind);cell=splitReg(reebCell(ind(1)));reebT(ind(1),:)=[0,0];end
            orgcell=cell;

            if ~wall_fol(i)

                [subXY,~] = BoustrophedonPath(cell,orgcell,EE,sensor,0,[],wall_fol(i),known,fl_see);
                if ~isempty(subXY);ccNode=[ccNode;subXY(1,:);subXY(end,:)];ccW=[ccW;total_length(rmmissing(subXY))];clear subXY;end

                [subXY,~] = BoustrophedonPath(cell,orgcell,EE,sensor,1,[],wall_fol(i),known,fl_see);
                if ~isempty(subXY);ccNode=[ccNode;subXY(1,:);subXY(end,:)];ccW=[ccW;total_length(rmmissing(subXY))];clear subXY;end

            else
                [subXY,~] = BoustrophedonPath(cell,orgcell,EE,sensor,0,fliplr(allNode(EE(1),:)),wall_fol(i),known,fl_see);
                if ~isempty(subXY);ccNode=[ccNode;subXY(1,:);subXY(end,:)];ccW=[ccW;total_length(rmmissing(subXY))];clear subXY;end

                [subXY,~] = BoustrophedonPath(cell,orgcell,EE,sensor,1,fliplr(allNode(EE(1),:)),wall_fol(i),known,fl_see);
                if ~isempty(subXY);ccNode=[ccNode;subXY(1,:);subXY(end,:)];ccW=[ccW;total_length(rmmissing(subXY))];clear subXY;end
            end

        else
            subXY = [fliplr(allNode(EE(1),:));fliplr(allNode(EE(2),:))];
            if ~isempty(subXY);ccNode=[ccNode;subXY(1,:);subXY(end,:);subXY(1,:);subXY(end,:)];ccW=[ccW;total_length(rmmissing(subXY));total_length(rmmissing(subXY))];clear subXY;end
        end
    end
    i=1:2:length(ccNode(:,1));ccEdge=[i',(i+1)'];
    sDir=repmat([0,1],1,length(i));
    i=2:4:length(ccNode(:,1))-3;ccEdge = [ccEdge; i',(i+3)';i',(i+5)';(i+2)',(i+3)';(i+2)',(i+5)'];

    ccWeight = spdist2(ccNode(ccEdge(:,1),:),ccNode(ccEdge(:,2),:));ccWeight(1:length(ccW))=ccW;
    ccG=digraph(ccEdge(:,1),ccEdge(:,2),ccWeight);
    ccDE = [ccEdge; length(ccNode(:,1)),1;];
    ccDG = sparse(ccEdge(:,1),ccEdge(:,2),ccWeight);
    [path1(1,:),d1(1)] = shortestpath(ccG,1,length(ccNode(:,1)));if ~known;d1(1)=d1(1)+spdist2(init,ccNode(1,:));end
    [path1(2,:),d1(2)] = shortestpath(ccG,1,length(ccNode(:,1))-2);if ~known;d1(2)=d1(2)+spdist2(init,ccNode(1,:));end
    [path1(3,:),d1(3)] = shortestpath(ccG,3,length(ccNode(:,1)));if ~known;d1(3)=d1(3)+spdist2(init,ccNode(3,:));end
    [path1(4,:),d1(4)] = shortestpath(ccG,3,length(ccNode(:,1))-2);if ~known;d1(4)=d1(4)+spdist2(init,ccNode(3,:));end
    [~,l]=min(d1);

    i=1:2:length(ccNode(:,1));
    ccPath = path1(l,1:2:length(path1(l,:)));clear path1 d1
    se=sDir(ismember(i,ccPath)); se=se(sum(reshape(ismember([Path(:,1:2);fliplr(Path(:,1:2))],reebEdge,'rows'),[],2),2)>0);
    seP=ccNode(ccPath,:);seP=seP(sum(reshape(ismember([Path(:,1:2);fliplr(Path(:,1:2))],reebEdge,'rows'),[],2),2)>0,:);
    %%%


    for i = 1:length(Path(:,1))
        EE = Path(i,1:2);
        if ismember(EE,[reebEdge;reebEdge(:,2),reebEdge(:,1)],'rows')
            o=o+1;
            ind=find(ismember(reebEdge,[EE;fliplr(EE)],'rows'));
            if ~isempty(ind);cell=splitReg(reebCell(ind(1)));reebEdge(ind(1),:)=[0,0];end

            orgcell=cell;

            [subXY,~] = BoustrophedonPath(cell,orgcell,EE,sensor,se(o),init,wall_fol(i),known,fl_see);

            if ~sim
                subXY(subXY(:,1)<sensor,1)=sensor;                  %%% Add the sensor gap
                subXY(subXY(:,1)>3048-sensor,1)=3048-sensor;
                subXY(subXY(:,2)<sensor,2)=sensor;                  %%% Add the sensor gap
                subXY(subXY(:,2)>2898-sensor,2)=2898-sensor;
            end

        else
            if isempty(PathEdge); PathEdge=fliplr(allNode(EE(1),:));end
            if (i~=size(Path,1)) && ismember(Path(i+1,1:2),[reebEdge;reebEdge(:,2),reebEdge(:,1)],'rows') && ~isempty(seP)
                [marker_x1,marker_y1] = addPtsLin([PathEdge(end,1),seP(o+1,1)],[PathEdge(end,2),seP(o+1,2)],s);
                subXY =[PathEdge(end,:);[marker_x1;marker_y1]';seP(o+1,:)];
            else
                if known;[subXY] = [PathEdge(end,:);fliplr(allNode(EE(2),:))];else;[subXY] = PathEdge(end,:);end
            end

            if ~sim
                % Use crack_reach when driving to a crack node (near-wall fill); transit/reconnect edges keep the sensor inset.
                if ~isempty(crack_reach) && EE(2)<=n_crack; bound=crack_reach; else; bound=sensor; end
                subXY(subXY(:,1)<bound,1)=bound;
                subXY(subXY(:,1)>3048-bound,1)=3048-bound;
                subXY(subXY(:,2)<bound,2)=bound;
                subXY(subXY(:,2)>2898-bound,2)=2898-bound;
            end
        end

        if ~isempty(subXY)
            if ~isempty(PathEdge) && ~known
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
