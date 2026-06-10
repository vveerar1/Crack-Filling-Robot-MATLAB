%% Boustrophedon Path Generation for a Single Cell
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [subXY,flag] = BoustrophedonPath(cell,orgcell,reebEdge,bp_gap,dir,init,wall_fol,known,fl_see)
%
% Generates a zig-zag boustrophedon path within a single monotone cell.
% Optionally prepends a wall-follow pass when the cell topology requires it.
%
% INPUTS:
%   cell     - Cell polygon eroded by the robot footprint (polyshape).
%   orgcell  - Original (un-eroded) cell polygon.
%   reebEdge - [start, end] critical-point indices bounding this cell.
%   bp_gap   - Sweep-line spacing (sensor range inscribed in square).
%   dir      - Sweep direction (0 or 1); [] to auto-select from init.
%   init     - Current robot position, or 0 if none.
%   wall_fol - 1 if this cell requires a wall-follow pass first.
%   known    - Flag: 1 for SCC (known environment).
%   fl_see   - Simulation direction flag.
%
% OUTPUTS:
%   subXY    - Sequential (x, y) waypoints for this cell.
%   flag     - True when a wall-follow pass was NOT performed.

    global spdist s allNode a total_length
    subXY=[];rowBW=max(cell.Vertices(:,2));flag=~wall_fol;cell_loj=false;
    if wall_fol
        subXY = [subXY;wall_follow(cell,orgcell,bp_gap,fliplr(allNode(reebEdge(1),:)))];
        if ~isempty(subXY)
            [~,lo]=min(abs(subXY(1,1)-fliplr(allNode([reebEdge(1),reebEdge(2)],2))));
            if lo==1
                if ~isempty(subXY)
                    obj = polybuffer([init;subXY],'line',bp_gap*sqrt(2));     % Keep in mind
                    cell = subtract(cell,obj);
                    aa=regions(cell);
                    if length(aa)>1
                        cell=aa(aa.area==max(aa.area));
                    end
                    init=subXY(end,:);
                    dir=[];
                end
                Start = reebEdge(2); End=reebEdge(1);
            else
                subXY=[];flag=wall_fol;
                Start = reebEdge(2); End=reebEdge(1);
            end
        else
            Start = reebEdge(1); End=reebEdge(2);
        end
    else
        Start = reebEdge(1); End=reebEdge(2);
    end
    if ~cell.area==0 && cell.area > area(polybuffer([0,0],'points',a)) % Cell area to footprint area check
        if allNode(Start,2)>allNode(End,2)
           if max(cell.Vertices(:,1))-min(cell.Vertices(:,1))<=2*bp_gap
               j = cell.centroid;
           else
               j=max(cell.Vertices(:,1))-bp_gap:-2*bp_gap:min(cell.Vertices(:,1))+bp_gap;
               if isempty(j);j=mean([max(cell.Vertices(:,1)),min(cell.Vertices(:,1))]);
               else
                if (j(end)- min(cell.Vertices(:,1))>bp_gap*sqrt(2))
                    j = [j,min(cell.Vertices(:,1))+bp_gap];
                end
               end
           end
        else
           if max(cell.Vertices(:,1))-(min(cell.Vertices(:,1)))<=2*bp_gap
                j = cell.centroid;
           else
               j=min(cell.Vertices(:,1))+bp_gap:+2*bp_gap:max(cell.Vertices(:,1))-bp_gap;
               if isempty(j);j=mean([max(cell.Vertices(:,1)),min(cell.Vertices(:,1))]);
               else
                if (max(cell.Vertices(:,1)) - j(end)>bp_gap*sqrt(2))
                    j = [j,max(cell.Vertices(:,1))-bp_gap];
                end
               end
           end
        end

        jx=j;
        jy=[min(cell.Vertices(:,2));max(cell.Vertices(:,2))];

        ins=[];dis=[];
        for j=jx
            %%New
            in_p=intersect(cell,polybuffer([[j;j],jy],'line',s));
            in=[j,min(in_p.Vertices(:,2));j,max(in_p.Vertices(:,2))];
            ins=[ins,in(:,2)];
            dis=[dis,total_length(in)];
        end

        ins(:,dis>=2*bp_gap)=ins(:,dis>=2*bp_gap)+[+bp_gap ;-bp_gap];
        ins(:,~(dis>=2*bp_gap))=[mean(ins(:,~(dis>=2*bp_gap)));mean(ins(:,~(dis>=2*bp_gap)))];

        if isempty(dir)
            if init==0
                dir=0;dir=dir==2;
            else
                [~,dir] = min(spdist(init,[[jx(1);jx(1)],ins(:,1)]));dir=dir==2;
            end
        end

        if dir
            ins(:,mod(1:length(jx),2)~=0)=flipud(ins(:,mod(1:length(jx),2)~=0));
        else
            ins(:,mod(1:length(jx),2)==0)=flipud(ins(:,mod(1:length(jx),2)==0));
        end

        for k=1:length(jx)
            if ins(1,k)<ins(2,k) && dis(k)>bp_gap
                in_space=ins(1,k):bp_gap:ins(2,k);if ~(in_space(end)==ins(2,k));in_space=[in_space,ins(2,k)];end
            elseif ins(1,k)>ins(2,k) && dis(k)>bp_gap
                in_space=ins(1,k):-bp_gap:ins(2,k);if ~(in_space(end)==ins(2,k));in_space=[in_space,ins(2,k)];end
            else
                in_space=ins(:,k)';
            end
            subXY=[subXY;[repmat(jx(k),length(in_space),1),in_space']];
        end
    else
        [subXY(1),subXY(2)] = cell.centroid;
    end
end

function subXY = wall_follow(polyin,orgpolyin,sensor,start)                                 % Wall Follow
    global vertical spdist s a
    if isequal(orgpolyin.Vertices,polyin.Vertices)
        working = polybuffer(polyin,-sensor);ss=sensor;
        if ~area(working)==0
            while 1
                if area(polyclean(working))==0 || size(regions(working),1)>1
                    ss=ss-10;
                    working = polybuffer(polyin,-ss);
                else
                    break
                end
            end
        end
    else
        working = polyin;
    end
    subXY=[];
    if working.area~=0
        [corPtx,idx] = polycorner(working);
        [~,vertexid]=min(spdist(start,corPtx));vertexid=idx(vertexid);
        temp=abs(working.Vertices(:,1)-working.Vertices(vertexid,1))<2;
        temp=working.Vertices(temp,:);
        [~,t(1)]=min(temp(:,2));
        [~,t(2)]=max(temp(:,2));

        [~,t(3)]=min(spdist(working.Vertices(vertexid,:),temp(t,:)));
        [~,vertexid]=min(spdist(temp(t(t(3)),:),corPtx));vertexid=idx(vertexid);
        ind = [vertexid,mod(vertexid,length(working.Vertices))+1];
        if ~vertical(working.Vertices(ind,:))
            subXY=[subXY;working.Vertices(ind(1),:)];
            while true
                vid=ind(2);
                subXY = [subXY;working.Vertices(vid,:)];
                ind2 = [ind,mod(ind(2),length(working.Vertices))+1];
                ind = [ind(2),mod(ind(2),length(working.Vertices))+1];
                if sum(ind(1)==idx)>0 && (fix(working.Vertices(ind(1),1))==fix(min(working.Vertices(:,1))) || fix(max(working.Vertices(:,1)))==fix(working.Vertices(ind(1),1)))
                    if size(subXY,1)==2
                       if fix(subXY(1,1))==fix(subXY(2,1));subXY=[];
                       elseif spdist(subXY(1,:),subXY(2,:))<a;subXY=[];end
                    end
                    break
                end
            end
        else
            ind = [vertexid,mod(vertexid-2,length(working.Vertices))+1];
            if ~vertical(working.Vertices(ind,:))
                subXY=[subXY;working.Vertices(ind(1),:)];
                while true
                    vid=ind(2);
                    subXY = [subXY;working.Vertices(vid,:)];
                    ind2 = [ind,mod(ind(2)-2,length(working.Vertices))+1];
                    ind = [ind(2),mod(ind(2)-2,length(working.Vertices))+1];
                    if sum(ind(1)==idx)>0 && (fix(working.Vertices(ind(1),1))==fix(min(working.Vertices(:,1))) || fix(max(working.Vertices(:,1)))==fix(working.Vertices(ind(1),1)))
                        if size(subXY,1)==2
                           if fix(subXY(1,1))==fix(subXY(2,1));subXY=[];
                           elseif spdist(subXY(1,:),subXY(2,:))<a;subXY=[];end
                        end
                        break
                    end
                end
            end
        end
    end
    if ~isempty(subXY)
        if fix(subXY(end,1))==fix(subXY(end-1,1))
            subXY(end,:)=[];
        end
    end

    marker_x2=[];marker_y2=[];
    for h=1:size(subXY,1)-1
        [m,n] = addPtsLin(subXY([h h+1],1)',subXY([h h+1],2)',s);
        marker_x2=[marker_x2,[subXY(h,1) m]];marker_y2=[marker_y2,[subXY(h,2) n]];
    end
    if ~isempty(subXY)
        subXY=[[marker_x2;marker_y2]';subXY(end,:)];
    end

    % Reduces number of points in the wall follow
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

function [corPtx,idx] = polycorner(polyin)
    corPtx=[];idx=[];global vertical
    i=convhull(rmmissing(polyin.Vertices));
    corPtx=[corPtx;polyin.Vertices(i,:)];
	idx=[idx;i];
end

function [marker_x,marker_y] = addPtsLin(x,y,marker_dist)
    % Inserts equidistant waypoints along a straight segment.
    dist_from_start = cumsum( [0, sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2)] );marker_x=[];marker_y=[];
    marker_locs = marker_dist : marker_dist : dist_from_start(end);
    if ~isempty(marker_locs) && (length(rmmissing(x))>1)
        marker_indices = interp1( dist_from_start, 1 : length(dist_from_start), marker_locs);
        marker_base_pos = floor(marker_indices);
        weight_second = marker_indices - marker_base_pos;
        marker_x = [marker_x, x(marker_base_pos) .* (1-weight_second) + x(marker_base_pos+1) .* weight_second];
        marker_y = [marker_y, y(marker_base_pos) .* (1-weight_second) + y(marker_base_pos+1) .* weight_second];
    end
end
