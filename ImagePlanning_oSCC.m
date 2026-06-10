%% Image Planning Function for OnlineSCC
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [waypoint_coords,flag,crackR,ttt] = ImagePlanning_oSCC(BW3,a,s,cp,acp,endP,realEndP,contEndP,ppath,preCrack)
% ImagePlanning_oSCC  Plans the crack-filling waypoints for one detected crack cluster.
%
% Given the current sensor binary image (BW3) and the pre-identified skeleton
% endpoints, this function extracts crack polylines, builds a crack graph via
% Minkowski-sum footprint buffering and a visibility graph, and computes a
% Chinese-Postman tour over that graph.  Waypoints are returned in traversal
% order, densified at footprint-radius spacing, with a fill flag indicating
% whether the robot nozzle should be active.
%
% INPUTS:
%   BW3        - Binary image of the sensor-range region (rows x cols).
%   a          - Nozzle footprint radius (pixels).
%   s          - Sensor detection radius (pixels).
%   cp         - Simulated current robot position [row, col].
%   acp        - Actual current robot position [row, col].
%   endP       - All skeleton endpoints of the extracted cracks [row, col].
%   realEndP   - Subset of endP that are true crack endpoints [row, col].
%   contEndP   - Subset of endP that are continuing (boundary) endpoints [row, col].
%   ppath      - Robot path from start to the current iteration [row, col].
%   preCrack   - Crack polylines from the previous iteration (cell array).
%
% OUTPUTS:
%   waypoint_coords - Crack-fill waypoints [row, col, fill_flag], densified at
%                     footprint-radius spacing; fill_flag=1 where the nozzle
%                     should dispense.
%   flag            - 1 if the robot is actively filling a crack, 0 if traversing.
%   crackR          - Crack polylines extracted in this iteration (cell array).
%   ttt             - Computation time (seconds).

    if sum(sum(BW3))>a/4
        spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);
        midP= @(P1,P2) (P1(:)+ P2(:)).'/2;
        total_length = @(Ps) sum(sqrt(sum(diff(Ps).*diff(Ps),2)));
        aP=polybuffer([0,0],'points',a);
        BW=BW3;
        BW3=padarray(BW3,[0,1],'post');
        BW3=padarray(BW3,[1,0],'post');
        worki=false;
        
        I=BW3;I_te=I;
        
        dir_map=[-1 -1;-1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0 ; 1 1];
        dir_can=[1 2 4 3 6;2 1 3 4 5; 3 2 5 1 8; 4 1 6 2 7 ; 5 3 8 2 7; 6 4 7 1 8; 7 6 8 4 5; 8 5 7 3 6];
        colors={'y','m','c','r','g','b'};
        
        [crackRaw,line,pointX,pointY] = compCrack(I,endP,dir_map,colors);

        tic
        % Discard crack segments shorter than the footprint radius (spur removal).
        c_chq=cellfun(total_length,crackRaw);c_chq=c_chq<a;
        line(c_chq,:)=[];
        pointX(c_chq,:)=[];pointY(c_chq,:)=[];
        crackRaw(c_chq)=[];

        %% Endpoint validation
        
        if ~isempty(realEndP)
            endlogi = reshape(sum(cell2mat(cellfun(@(s) spdist(s,[line(:,1) line(:,2);line(:,3) line(:,4)]),num2cell(realEndP, 2),'un',0)')<5,2)>0,[],2);
            endlogi_a = reshape(spdist(cp,[line(:,1) line(:,2);line(:,3) line(:,4)])<a,[],2);
        else
            endlogi = zeros(size(line,1),2);endlogi_a = zeros(size(line,1),2);
        end
        line(~(endlogi(:,1)|endlogi(:,2)),:)=[];crackRaw(~(endlogi(:,1)|endlogi(:,2)))=[];
        pointX(~(endlogi(:,1)|endlogi(:,2)),:)=[];pointY(~(endlogi(:,1)|endlogi(:,2)),:)=[];
        endlogi_a(~logical(sum(endlogi,2)),:)=[];endlogi(~logical(sum(endlogi,2)),:)=[];
        link=struct('x',cell(1,size(line,1)),'y',cell(1,size(line,1)));

        %% Minkowski-sum footprint offset — acute-angle crack adjustment
        % For cracks with a sharp apex angle, reduce the buffer radius to
        % tan(half-angle)*a so the offset polygon fits inside the crack wedge.
        fP=repmat(a,[size(line,1),1]);
        [~,yy]=size(BW3);det=[];
        for l=1:length(crackRaw)
            if total_length(crackRaw{l})>2*a*sqrt(2)
                aRan=polybuffer(crackRaw{l}([1,end],:),'points',a);%plot(polybuffer(fliplr(crackRaw{l}([1,end],:)),'points',a))
                crackW=crackRaw{l}(~ismember(crackRaw{l},rmmissing(intersect(aRan,rmmissing(crackRaw{l}))),'rows'),:);
                if total_length(crackW)>a % New ~isempty(crackW)
                    mov=crackW-crackW(1,:);
                else
                    mov=crackRaw{l}-crackRaw{l}(1,:);
                end
            else
                mov=crackRaw{l}-crackRaw{l}(1,:);%plot(mov(:,2),mov(:,1));plot(mov(end,2),mov(end,1),'r*')
            end
            v1=[0,yy];
            v2=mov(end,:);
            ang=ab2v(v1,v2);
            rotv=[cos(deg2rad(ang)) sin(deg2rad(ang)); -sin(deg2rad(ang)) cos(deg2rad(ang))];
            rv2=(rotv*[mov]')';%plot(rv2(:,2),rv2(:,1))
            [~,d]=min(rv2(:,1)); peakPts_idx=d;
            [~,d]=max(rv2(:,1)); peakPts_idx=[peakPts_idx;d];
            [~,d]=max(abs(rv2(peakPts_idx,1)));peakPts_idx=peakPts_idx(d);%plot(rv2(peakPts_idx,2),rv2(peakPts_idx,1),'r*')
            v1=-rv2(peakPts_idx,:);v2=rv2(end,:)-rv2(peakPts_idx,:);
            ang=ab2v(v1,v2)/2; % Half apex angle; if <45° the crack is acute
            if ang<45 || ang >360-45
                fP(l)=tan(deg2rad(ang))*a;det=[det;l];
            end
        end

        % Retain only acute-crack edges that have at least one true endpoint.
        det=det(sum(endlogi(det,:))>0);h_line=[]; % Check at least one node is a real endpoint
        if ~isempty(det)
            h_line = [h_line;line(det,:)];
            h_pointX=pointX(det,:);h_pointY=pointY(det,:);
            h_crackRaw=crackRaw(det);h_fP=fP(det);

            line(det,:)=[];link(det)=[];
            pointX(det,:)=[];pointY(det,:)=[];
            crackRaw(det)=[];fP(det)=[]; endlogi(det,:)=[]; endlogi_a(det,:)=[];
        end

        %% Connecting the crack graph

        if ~isempty(line) || ~isempty(h_line)
            if ~isempty(line)
                detind=[];
                templine=line;

                endlogi_t = reshape(spdist(cp,[line(:,1) line(:,2);line(:,3) line(:,4)])<s*sqrt(2),[],2) & endlogi;

                ll= sum(endlogi,2)>0;

                tt=reshape([templine(ll,[1,3]),templine(ll,[2,4])],[],2);

                if isequal(cp,acp)
                ll= sum(endlogi_t,2)==max(sum(endlogi_t,2));
                end

                endlogi_t=endlogi;endlogi_t(~(ll),:)=0;ee=endlogi_t;
                tt2=tt(logical(ee(:)),:);  % True end points
                [cd,ee1] = min(spdist(cp,tt2));[~,ee]=ismember(tt2(ee1,:),tt,'rows');

                % Resolve intersection points: skip any endpoint shared by multiple edges.
                while(1)
                    if sum(ismember(tt,tt(ee,:),'rows'))>1
                        tt2(ismember(tt2,tt2(ee1,:),'rows'),:)=[];
                        [cd,ee1] = min(spdist(cp,tt2));[~,ee]=ismember(tt2(ee1,:),tt,'rows');
                    else
                        break
                    end
                end

                if cd<=s*sqrt(2)
                    ee=spdist(tt(ee,:),realEndP);ee=ee==min(ee);
                    [~,d]=dsearchn(realEndP(ee,:),[templine(:,1) templine(:,2);templine(:,3) templine(:,4)]);d=d<6;d=reshape(d,[],2);
                    d(~(sum(d,2)&ll),:)=0;ll=find(sum(d,2));pp1=[];

                    if length(ll)>1
                        for lll=ll'
                            [~,pp]=dsearchn(crackRaw{lll},acp);pp1=[pp1;pp];
                        end
                        d(ll(~(min(pp1)==pp1)),:)=0;
                    end

                    if any(d,'all')%any(d(:,1))
                        intt= find(d(:,1));if isempty(intt); intt= find(d(:,2));else intt=intt(1);end 
                        crackRawt=crackRaw;
                        if find(d(intt,:))==1
                            link(intt).x=pointX(intt,:);
                            link(intt).y=pointY(intt,:);
                        else
                            for intt_l=reshape(intt,1,[])
                                link(intt_l).x=fliplr(pointX(intt_l,:));
                                link(intt_l).y=fliplr(pointY(intt_l,:));
                                endlogi(intt_l,:)= fliplr(endlogi(intt_l,:));
                                crackRaw{intt_l}=flipud(crackRaw{intt_l});
                                line(intt_l,:)=[line(intt_l,(3:4)),line(intt_l,(1:2))];
                            end
                        end
                        templine(intt,:)=0;indt=intt;

                        if endlogi(intt,2)~=0 && isequal(cp,acp)   % Both endpoints are real — chain to the next adjacent segment
                            endlogi = reshape(ismember([templine(:,1) templine(:,2);templine(:,3) templine(:,4)],realEndP,'rows'),[],2);
                            ss=sum(endlogi,2);
                            if any(ss(1:size(ss,1),:)==2)             % Another segment with two real endpoints exists
                                intt2=find(ss(1:size(ss,1),:)==2)';
                                detind=[detind,find(~ismember(1:size(line,1),intt2) & ~ismember(1:size(line,1),intt))];
                                templine(setdiff(1:end,intt2),:)=0;

                                while(1)
                                    tempdist= spdist([link(intt).x(end),link(intt).y(end)],templine(:,1:2));
                                    inst=find(tempdist<=a);
                                    if ~isempty(inst)
                                        inst=find(tempdist==min(tempdist(tempdist<=a)));inst=inst(1);
                                        v1=crackRawt{indt}(1,:)-crackRawt{indt}(end,:);
                                        v2=crackRaw{intt}(1,:)-crackRaw{intt}(end,:);
                                        v3=crackRaw{inst}(end,:)-crackRaw{inst}(1,:);
                                        ang=ab2v(v1,v3);indt=inst;ang2=ab2v(v2,v3);
                                    end
                                    if ~isempty(inst) && ~(ang<45 || ang >360-45) && ~(ang2<45 || ang2 >360-45)
                                        [~,d]=min(tempdist);
                                        link(intt).x=[link(intt).x, pointX(d,:)];
                                        link(intt).y=[link(intt).y, pointY(d,:)];
                                        line(intt,3:4)=[link(intt).x(end),link(intt).y(end)];
                                        crackRaw{intt}=[crackRaw{intt};crackRaw{d}];
                                        templine(d,:)=0;
                                        detind=[detind,d];intt2(intt2==d)=[];
                                    else
                                        tempdist= spdist([link(intt).x(end),link(intt).y(end)],templine(:,3:4));
                                        inst=find(tempdist<=a);
                                        if ~isempty(inst)
                                            inst=find(tempdist==min(tempdist(tempdist<=a)));inst=inst(1);
                                            v1=crackRawt{indt}(1,:)-crackRawt{indt}(end,:);
                                            v2=crackRaw{intt}(1,:)-crackRaw{intt}(end,:);
                                            v3=crackRaw{inst}(1,:)-crackRaw{inst}(end,:);
                                            ang=ab2v(v1,v3);indt=inst;ang2=ab2v(v2,v3);
                                        end
                                        if ~isempty(inst) && ~(ang<45 || ang >360-45) && ~(ang2<45 || ang2 >360-45)
                                            [~,d]=min(tempdist);
                                            link(intt).x=[link(intt).x, fliplr(pointX(d,:))];
                                            link(intt).y=[link(intt).y, fliplr(pointY(d,:))];
                                            line(intt,3:4)=[link(intt).x(end),link(intt).y(end)];
                                            crackRaw{intt}=[crackRaw{intt};flipud(crackRaw{d})];
                                            crackRawt{d}=flipud(crackRawt{d});
                                            templine(d,:)=0;
                                            detind=[detind,d];intt2(intt2==d)=[];
                                        else
                                            intt=intt2(1);indt=intt;
                                            link(intt).x=[link(intt).x, pointX(intt,:)];
                                            link(intt).y=[link(intt).y, pointY(intt,:)];
                                            line(intt,3:4)=[link(intt).x(end),link(intt).y(end)];
                                            templine(intt,:)=0;
                                            intt2(intt2==intt)=[];
                                        end
                                    end
                                    if any(templine,'all')
                                        continue
                                    else
                                        break
                                    end
                                end
                            elseif ss==0
                                detind=[detind,find(1:size(d,1)~=intt)];
                            else
                                ddist=[];
                                for ll=reshape(find((1:size(endlogi,1)~=intt)' & ss),1,[])
                                        ii = find(endlogi(ll,:));
                                        ddist=[ddist;ll,spdist([link(intt).x(end),link(intt).y(end)],templine(ll,2*ii-1:2*ii))];
                                end
                                [n,m]=min(ddist(:,2));%%%
                                if n<=a
                                    if find(endlogi(ddist(m,1),:))==1
                                        link(intt).x=[link(intt).x,pointX(ddist(m,1),:)];
                                        link(intt).y=[link(intt).y,pointY(ddist(m,1),:)];
                                        line(intt,3:4)=[link(intt).x(end),link(intt).y(end)];
                                        crackRaw{intt}=[crackRaw{intt};crackRaw{ddist(m,1)}];
                                        detind=[detind,find(1:size(d,1)~=intt)];
                                    else
                                        link(intt).x=[link(intt).x,fliplr(pointX(ddist(m,1),:))];
                                        link(intt).y=[link(intt).y,fliplr(pointY(ddist(m,1),:))];
                                        line(intt,3:4)=[link(intt).x(end),link(intt).y(end)];
                                        crackRaw{intt}=[crackRaw{intt};flipud(crackRaw{ddist(m,1)})];
                                        detind=[detind,find(1:size(d,1)~=intt)];
                                    end
                                else
                                    detind=[detind,find(1:size(d,1)~=intt)];
                                end
                            end





                        elseif endlogi(intt,2)==0 || ~isequal(cp,acp)   % Tail endpoint is a continuing point — keep only this segment
                            detind=[detind,find(1:size(d,1)~=intt)];
                        end

                    end

                elseif max(sum(endlogi,2))==2 && isequal(cp,acp)
                    ll= sum(endlogi,2)==max(sum(endlogi,2));
                    tt=reshape([templine(ll,[1,3]),templine(ll,[2,4])],[],2);

                    ee=endlogi(ll,:);
                    tt2=tt(ee(:),:);

                    [cd,ee1] = min(spdist(cp,tt2));[~,ee]=ismember(tt2(ee1,:),tt,'rows');

                    % Skip intersection points shared by multiple edges.
                    while(1)
                        if sum(ismember(tt,tt(ee,:),'rows'))>1
                            tt2(ismember(tt2,tt2(ee1,:),'rows'),:)=[];
                            [cd,ee1] = min(spdist(cp,tt2));[~,ee]=ismember(tt2(ee1,:),tt,'rows');
                        else
                            break
                        end
                    end

                    [~,ee] = dsearchn(tt(ee,:),realEndP);ee=ee<6;
                    [~,d]=dsearchn(realEndP(ee,:),[templine(:,1) templine(:,2);templine(:,3) templine(:,4)]);d=d<6;d=reshape(d,[],2);
                    d(~(sum(d,2)&ll),:)=0;ll=find(sum(d,2));pp1=[];

                    if length(ll)>1
                        for lll=ll'
                            [~,pp]=dsearchn(crackRaw{lll},acp);pp1=[pp1;pp];
                        end
                        d(ll(~(min(pp1)==pp1)),:)=0;
                    end

                    if any(d,'all')
                        intt= find(d(:,1));if isempty(intt); intt= find(d(:,2));end
                        if find(d(intt,:))==1
                            link(intt).x=pointX(intt,:);
                            link(intt).y=pointY(intt,:);
                        else
                            link(intt).x=fliplr(pointX(intt,:));
                            link(intt).y=fliplr(pointY(intt,:));
                            endlogi(intt,:)= fliplr(endlogi(intt,:));
                            crackRaw{intt}=flipud(crackRaw{intt});
                        end
                        detind=[detind,find(1:size(d,1)~=intt)];
                    else
                        detind=1:size(line,1);
                    end

                else
                    detind=1:size(line,1);
                end

                line(detind,:)=[];
                link(all(cell2mat(arrayfun(@(x) structfun(@isempty, x), link, 'UniformOutput', false)),1)) = []; %link(detind)=[];
                crackRaw(detind)=[];
                
                fP(detind)=[];
            end

            % Reintroduce the acute-crack edges that were held aside.
            if ~isempty(h_line)
                line = [line;h_line];
                crackRaw=[crackRaw,h_crackRaw];fP=[fP;h_fP];
                for kk=1:size(h_line,1)
                    link(end+1).x=fliplr(h_pointX(kk,:));
                    link(end).y=fliplr(h_pointY(kk,:));
                end
            end
            crackR=crackRaw;
            if ~isempty(link) && worki
                for p=1:length(link);plot(link(p).y,link(p).x,'LineWidth',2);end
            end
            
            if ~isempty(preCrack)
                if isequal(crackRaw,preCrack)
                    fline=[];  %to skip the whole function
                end
            end
            
          
           %% Minkowski-sum footprint offset — second pass after graph connection
            for l=1:length(crackRaw)
                if total_length(crackRaw{l})>2*a*sqrt(2)
                    aRan=polybuffer(crackRaw{l}([1,end],:),'points',a);
                    crackW=crackRaw{l}(~ismember(crackRaw{l},rmmissing(intersect(aRan,rmmissing(crackRaw{l}))),'rows'),:);
                    if total_length(crackW)>a
                        mov=crackW-crackW(1,:);
                    else
                        mov=crackRaw{l}-crackRaw{l}(1,:);
                    end
                else
                    mov=crackRaw{l}-crackRaw{l}(1,:);
                end
                v1=[0,yy];
                v2=mov(end,:);
                ang=ab2v(v1,v2);
                rotv=[cos(deg2rad(ang)) sin(deg2rad(ang)); -sin(deg2rad(ang)) cos(deg2rad(ang))];
                rv2=(rotv*[mov]')';
                [~,d]=min(rv2(:,1)); peakPts_idx=d;
                [~,d]=max(rv2(:,1)); peakPts_idx=[peakPts_idx;d];
                [~,d]=max(abs(rv2(peakPts_idx,1)));peakPts_idx=peakPts_idx(d);
                v1=-rv2(peakPts_idx,:);v2=rv2(end,:)-rv2(peakPts_idx,:);
                ang=ab2v(v1,v2)/2; % Half apex angle
                if (ang<45 || ang >360-45)&& ang~=0
                    fP(l)=tan(deg2rad(ang))*a;
                end
            end
           
           %% Build the crack graph
            if ~isempty(line)
                y=[line(:,1);line(:,3)];
                x=[line(:,2);line(:,4)];

                if worki;for b=1:length(link);plot(link(b).y,link(b).x,'LineWidth',2);end;end
                ss=struct2cell(link);
                llink = [horzcat(ss{1,1,:});horzcat(ss{2,1,:})]';clear ss

                % Re-validate: keep only edges that touch a real endpoint.
                endlogi = reshape(ismember([line(:,1) line(:,2);line(:,3) line(:,4)],realEndP,'rows'),[],2);
                [~,d]=dsearchn(realEndP,[line(:,1) line(:,2);line(:,3) line(:,4)]);d=d<6;endlogi=endlogi|reshape(d,[],2);
                line(~(endlogi(:,1)|endlogi(:,2)),:)=[];
                link(~(endlogi(:,1)|endlogi(:,2)))=[];

                %% Identify isolated graph nodes (no nearby neighbor)
                % Collect all crack endpoints; retain those whose nearest other
                % endpoint on a different edge lies farther than the footprint
                % radius (these become crack-graph nodes directly).
                endPoints=unique([line(:,1) line(:,2);line(:,3) line(:,4)],'rows');endNodes=[];endPoints_t=endPoints;
                e=1;
                while(~isempty(endPoints_t))
                    endlogi = reshape(ismember([line(:,1) line(:,2);line(:,3) line(:,4)],endPoints_t(e,:),'rows'),[],2);
                    endlogi=sum(endlogi,2);
                    TestPoints = endPoints_t(~ismember(endPoints_t,[line(logical(endlogi),1:2);line(logical(endlogi),3:4)],'rows'),:);
                    if all(spdist(endPoints_t(e,:),TestPoints)>a)
                        endNodes = [endNodes;endPoints_t(e,:)]; endPoints_t(e,:)=[];
                    else
                        endPoints_t(e,:)=[];
                    end
                end

                %% Minkowski-sum footprint buffers
                % Expand each crack polyline by the (possibly reduced) footprint
                % radius to form a filled polygon over the fillable area.
                n=size(line,1);

                for i=1:n
                    poly = polybuffer([link(i).y',link(i).x'],'line',fP(i));
                    [px,py]=poly.boundary;
                    if ~(px(1)==px(end)&&py(1)==py(end));px(end+1)=px(1);py(end+1)=py(1);end
                    poly.Vertices=DecimatePoly([px,py],[1 1],false);
                    S1(i).P.x = poly.Vertices(:,2);
                    S1(i).P.y = poly.Vertices(:,1);
                    S1(i).P.hole = poly.NumHoles;
                    S(i).P(1) = poly;
                    if i==1
                        u=poly;
                        v(i)=poly;
                    else
                        u = union(u,poly);  % Accumulate the union of all buffer polygons
                        v(i)=poly;
                    end

                end

                %% Compute overlapping areas and derive graph nodes
                % Where two buffer polygons overlap, the centroid of the
                % intersection is a candidate graph node.  Nodes closer than
                % the footprint radius to each other are merged.
                TF = overlaps(v);
                if any(TF.*~eye(size(TF)),'all')
                    Display_result = 0;
                    Accuracy       = 1e-3;
                    Geo=Polygons_intersection_modified(S,Display_result,Accuracy);
                    % Compute centroids of overlapping areas; discard nodes
                    % within footprint distance of an already-accepted node.
                    N=length(Geo);
                    Ind=[];
                    nodes=[];
                    for i=1:N
                        index=Geo(i).index;
                        inum=length(index);
                        if inum>=2                      % Keep only pairwise (or higher) overlaps
                            pnum = length(Geo(i).P);
                            for j = 1:pnum
                                Ind=[Ind;inum];
                                [cy,cx]=centroid(Geo(i).P(j));
                                nodes=[nodes;cx,cy];   % Centroid of this overlapping region
                            end
                        end
                    end

                    [~,Index]=sort(Ind,'descend');
                    SortNode=nodes(Index,:);

                    NodeCan = [SortNode;endNodes];node=[];

                    for n=1:size(NodeCan,1)
                        if ~isempty(node)
                            if all(spdist(NodeCan(n,:),node)>a)
                                node = [node;NodeCan(n,:)];
                            end
                        else
                            node = [node;NodeCan(n,:)];
                        end
                    end


                else
                    node=endNodes;
                end

                %% Trim crack endpoints to the inscribed-footprint boundary
                % For each graph node that is a real crack endpoint, clip
                % the corresponding crack polyline so it ends at the
                % footprint boundary (radius a/sqrt(2)).
                for e=1:size(node,1)
                    aRan=polybuffer(node(e,:),'points',a/sqrt(2));
                    logi_cell=cell2mat(cellfun(@(x) ~isempty(x),cellfun(@(x) intersect(aRan,rmmissing(x)),crackRaw,'un',0),'un',0));
                    if sum(logi_cell)==1&&ismember(node(e,:),[realEndP;fliplr(realEndP)],'rows')
                        if any(isinterior(aRan,crackRaw{logi_cell}([1,end],:)))
                            [~,out]=intersect(aRan,rmmissing(crackRaw{logi_cell}));
                            if ~isempty(out)
                                crackRaw{logi_cell}=out;ee=crackRaw{logi_cell}([1,end],:);
                                ddd=min(spdist(node(e,:),ee))==spdist(node(e,:),ee);
                                node(e,:)=ee(ddd,:);
                            end
                        end
                    end
                end
   
                %% Visibility graph
                % For each buffer polygon, sort the graph nodes that fall
                % inside it along the crack polyline; test each candidate
                % edge for line-of-sight.  Occluded edges are replaced by
                % shortest obstacle-avoiding paths via pathfinder.

                v={};vgNE={};vgEE={};iT=[];
                for i = 1: length(S)
                    in = inpolygon(node(:,1),node(:,2),S(i).P.Vertices(:,2),S(i).P.Vertices(:,1));
                    if any(in)
                        v(i).P = [find(in),node(in,:)];
                        mem=ismember(v(i).P(:,2:3),[link(i).x',link(i).y'],'rows');
                        [~,d]=dsearchn([link(i).x',link(i).y'],v(i).P(:,2:3));d=d<a;
                        mem = mem|d;
                        temp = sortrows([dsearchn([link(i).x',link(i).y'],v(i).P(mem,2:3)),v(i).P(mem,:)],1,'ascend');
                        if length(temp(:,1))-1>0
                            for j = 1:length(temp(:,1))-1
                                vgNE(i).S(j,:) = [temp(j,3:4) temp(j+1,3:4)];
                                vgEE(i).S(j,:) = [temp(j,2) temp(j+1,2)];
                            end
                        else
                            vgNE(i).S=[];vgEE(i).S=[];
                        end

                        if S(i).P.NumHoles>0
                            [~,d]=dsearchn([link(i).x(temp(end,1):end)',link(i).y(temp(end,1):end)'],v(i).P(mem,2:3));d=d<a;mem1=d;
                            temp1 = sortrows([dsearchn([link(i).x',link(i).y'],v(i).P(mem1,2:3)),v(i).P(mem1,:)],1,'ascend');
                            if length(temp1(:,1))-1>0
                                for j = 1:length(temp1(:,1))-1
                                    vgNE(i).S = [vgNE(i).S;temp1(j,3:4) temp1(j+1,3:4)];
                                    vgEE(i).S = [vgEE(i).S;temp1(j,2) temp1(j+1,2)];
                                end
                            end
                        end
                        temp2 = dsearchn(temp(:,3:4),v(i).P(~mem,2:3));
                        vgNE(i).S = [vgNE(i).S;v(i).P(~mem,2:3),temp(temp2,3:4)];
                        vgEE(i).S = [vgEE(i).S;v(i).P(~mem,1),temp(temp2,2)];
                    else
                        iT=[iT,i];vgNE(i).S=[];vgEE(i).S=[];
                    end
                end
                S(iT)=[];vgEE(iT)=[];vgNE(iT)=[];

                for i=1:length(S)
                    visibility(i).S = line_of_sight_poly(vgNE(i).S(:,1:2),vgNE(i).S(:,3:4),[S(i).P.Vertices(:,2)'; S(i).P.Vertices(:,1)']);
                    for j=1:length(visibility(i).S)
                        if ~visibility(i).S(j)
                            start=vgNE(i).S(j,1:2);goal=vgNE(i).S(j,3:4);boundary=[S(i).P.Vertices(:,2)'; S(i).P.Vertices(:,1)']';
                            startn=vgEE(i).S(j,1);goaln=vgEE(i).S(j,2);
                            slen = length(node(:,1));

                            [waypoints,~]=pathfinder(start,goal,boundary);nn=waypoints(2:end-1,:);
                            node = [node;nn];elen = length(node(:,1));
                            if ~isempty(nn)
                                vgNE(i).S = [vgNE(i).S;start,nn(1,:)];
                                vgEE(i).S = [vgEE(i).S;startn,slen+1];
                                if length(nn(:,1))>1
                                    for k = 2:size(nn,1)
                                        vgNE(i).S = [vgNE(i).S;vgNE(i).S(end,3:4),nn(k,:)];
                                        vgEE(i).S = [vgEE(i).S;vgEE(i).S(end,2),slen+k];
                                    end
                                end
                                vgNE(i).S = [vgNE(i).S;vgNE(i).S(end,3:4),goal];
                                vgEE(i).S = [vgEE(i).S;vgEE(i).S(end,2),goaln];
                            else
                                vgNE(i).S = [vgNE(i).S;start,goal];
                                vgEE(i).S = [vgEE(i).S;startn,goaln];
                            end
                        end
                    end
                    vgNE(i).S(~logical(visibility(i).S),:)=[];vgEE(i).S(~logical(visibility(i).S),:)=[]; % Removes the edges
                end

                edgeList = [];
                for i=1: size(vgEE,2)
                    edgeList = [edgeList ;vgEE(i).S];
                end

                if worki
                    text(node(:,2)+10,node(:,1)+10,int2str([1:length(node)]'));
                    plot([node(edgeList(:,1),2) node(edgeList(:,2),2)]',[node(edgeList(:,1),1) node(edgeList(:,2),1)]','--')
                end

                edgeList = unique(sort(edgeList, 2),'rows');

                G=graph(edgeList(:,1),edgeList(:,2));

                %% Connect isolated nodes
                % Any graph node not yet linked to an edge is joined to its
                % nearest neighbor so the graph remains connected.
                while any(~ismember(1:numel(node(:,1)),unique(G.Edges.EndNodes)))
                    nnode=find(~ismember(1:numel(node(:,1)),unique(G.Edges.EndNodes)));
                    [~,ddd]=min(spdist(cp,node(nnode,:)));
                    dd= spdist(node(nnode(ddd),:),node);dd(dd==0)=inf;[~,dd]=min(dd);
                    G=addedge(G,nnode(ddd),dd);
                    edgeList=[edgeList;nnode(ddd),dd];
                end

                %% Chinese-Postman tour
                % Compute an Eulerian (or near-Eulerian) traversal of the
                % crack graph.  Continuing endpoints (boundary-contact nodes)
                % are excluded from the odd-degree matching so the tour starts
                % and ends at a real crack endpoint nearest to the robot.

                if  size(edgeList,1)>1
                    adj=full(adjacency(G));b_n=mod(sum(adj,2),2)~=0;
                    remNode = find(ismember(node,contEndP,'rows'));eeP=[line(:,[1,2]);line(:,[3,4])];
                    [~,d]=min(spdist(cp,eeP));
                    % Select start node from nodes that have at least one neighbor.
                    nodewNeighbor=node(unique(G.Edges.EndNodes(:)),:);
                    [~,d]=min(spdist(eeP(d,:),nodewNeighbor));
                    [~,d]=ismember(nodewNeighbor(d),node);

                    Dist =squareform(pdist(node));
                    AdjMax=adj.*Dist;
                    [Path, ~, add,st]=ChinesePostman(adj,AdjMax,Dist,[],remNode,d);
                    if ~isempty(st); add(ismember(add,st,'rows'),:)=[]; end
                else
                    [~,d]=dsearchn(cp,node);[~,d]=min(d);
                    Path = [edgeList(d),edgeList(1:end~=d)];
                end

                %% Assemble the waypoint sequence

                waypoint=node(Path,:);

                waypoint_coords=waypoint;
                [~,d1]=dsearchn(waypoint_coords(1,:),realEndP);d1=d1<=a+5;[~,d2]=dsearchn(waypoint_coords(end,:),realEndP);d2=d2<=a+5;
                endPP = [waypoint_coords(1,:),any(d1);waypoint_coords(end,:),any(d2)];

                %% Densify waypoints and set fill flag
                % Orientation: ensure the tour starts at the end nearest
                % to the robot; then insert equidistant intermediate points
                % at footprint-radius spacing and mark each point with a
                % fill flag (1 = nozzle active, 0 = traversal only).
                if any(endPP(:,3)) && a==a
                    flag = 1;
                    if endPP(1,3)~=endPP(2,3)
                        if endPP(2,3)
                            waypoint_coords=flipud(waypoint_coords);
                            [marker_x1,marker_y1] = addPtsLin([cp(1),waypoint_coords(1,1)],[cp(2),waypoint_coords(1,2)],a+5);
                            marker_x2=[];marker_y2=[];
                            for h=1:size(waypoint_coords,1)-1
                                [m,n] = addPtsLin(waypoint_coords([h h+1],1)',waypoint_coords([h h+1],2)',a);
                                marker_x2=[marker_x2,[waypoint_coords(h,1) m]];marker_y2=[marker_y2,[waypoint_coords(h,2) n]];
                            end
                            if ~isempty(marker_x1);b=1;else;b=[];end
                            dd=false(1,length(marker_x2));
                            for ii=1:length(crackRaw)
                                [~,t]=dsearchn(crackRaw{ii},[marker_x2;marker_y2]');dd=dd|t'<=a*sqrt(2);
                            end
                            waypoint_coords=[[marker_x1;marker_y1;[zeros(1,size(marker_x1,2)-1) b]]';[marker_x2;marker_y2;dd]';[waypoint_coords(end,:),1]];
                        else
                            [marker_x1,marker_y1] = addPtsLin([cp(1),waypoint_coords(1,1)],[cp(2),waypoint_coords(1,2)],a+5);
                            marker_x2=[];marker_y2=[];
                            for h=1:size(waypoint_coords,1)-1
                                [m,n] = addPtsLin(waypoint_coords([h h+1],1)',waypoint_coords([h h+1],2)',a);
                                marker_x2=[marker_x2,[waypoint_coords(h,1) m]];marker_y2=[marker_y2,[waypoint_coords(h,2) n]];
                            end
                            if ~isempty(marker_x1);b=1;else;b=[];end
                            dd=false(1,length(marker_x2));
                            for ii=1:length(crackRaw)
                                [~,t]=dsearchn(crackRaw{ii},[marker_x2;marker_y2]');dd=dd|t'<=a*sqrt(2);
                            end
                            waypoint_coords=[[marker_x1;marker_y1;[zeros(1,size(marker_x1,2)-1) b]]';[marker_x2;marker_y2;dd]';[waypoint_coords(end,:),1]];
                        end
                    else
                        [~,d]=min(spdist(cp,endPP(:,[1,2])));
                        if d==2
                            waypoint_coords=flipud(waypoint_coords);
                            [marker_x1,marker_y1] = addPtsLin([cp(1),waypoint_coords(1,1)],[cp(2),waypoint_coords(1,2)],a+5);
                            marker_x2=[];marker_y2=[];
                            for h=1:size(waypoint_coords,1)-1
                                [m,n] = addPtsLin(waypoint_coords([h h+1],1)',waypoint_coords([h h+1],2)',a);
                                marker_x2=[marker_x2,[waypoint_coords(h,1) m]];marker_y2=[marker_y2,[waypoint_coords(h,2) n]];
                            end
                            if ~isempty(marker_x1);b=1;else;b=[];end
                            dd=false(1,length(marker_x2));
                            for ii=1:length(crackRaw)
                                [~,t]=dsearchn(crackRaw{ii},[marker_x2;marker_y2]');dd=dd|t'<=a*sqrt(2);
                            end
                            waypoint_coords=[[marker_x1;marker_y1;[zeros(1,size(marker_x1,2)-1) b]]';[marker_x2;marker_y2;dd]';[waypoint_coords(end,:),1]];
                        else
                            [marker_x1,marker_y1] = addPtsLin([cp(1),waypoint_coords(1,1)],[cp(2),waypoint_coords(1,2)],a+5);
                            marker_x2=[];marker_y2=[];
                            for h=1:size(waypoint_coords,1)-1
                                [m,n] = addPtsLin(waypoint_coords([h h+1],1)',waypoint_coords([h h+1],2)',a);
                                marker_x2=[marker_x2,[waypoint_coords(h,1) m]];marker_y2=[marker_y2,[waypoint_coords(h,2) n]];
                            end
                            if ~isempty(marker_x1);b=1;else;b=[];end
                            dd=false(1,length(marker_x2));
                            for ii=1:length(crackRaw)
                                [~,t]=dsearchn(crackRaw{ii},[marker_x2;marker_y2]');dd=dd|t'<=a*sqrt(2);
                            end
                            waypoint_coords=[[marker_x1;marker_y1;[zeros(1,size(marker_x1,2)-1) b]]';[marker_x2;marker_y2;dd]';[waypoint_coords(end,:),1]];
                        end

                    end
                else
                    flag = 0;
                end
                
            else
                waypoint_coords=[];
                flag = 0;
            end
        else
            waypoint_coords=[];crackR=[];
            flag = 0;
        end
    else
    	waypoint_coords=[];crackR=[];
        flag = 0;
    end
ttt=toc;
end


%% Local helper functions

function ang = ab2v(a,b)
% ab2v  Compute the unsigned angle (degrees, 0–360) from vector a to vector b.
%   Rotates both vectors into a canonical frame aligned with a, then measures
%   the CCW angle to b.  Returns a value in [0, 360).
    theta = rad2deg(atan2(norm(cross([1,0,0],[a,0])), dot([1,0,0],[a,0])));
    if a(2)<0 theta=360-theta; end

    aR = a*[1;1i]*exp(-1i*theta*pi/180); aR=[real(aR) imag(aR)];
    bR = b*[1;1i]*exp(-1i*theta*pi/180); bR=[real(bR) imag(bR)];

    ang=rad2deg(atan2(norm(cross([aR,0],[bR,0])), dot([aR,0],[bR,0])));
    if bR(2)<0 ang=360-ang; end
end

function y = bound(x,bl,bu)
% bound  Clip scalar x to the interval [bl, bu].
    y=min(max(x,bl),bu);
end

function [marker_x,marker_y] = addPtsLin(x,y,marker_dist)
% addPtsLin  Insert equidistant intermediate points along a polyline.
%   (x,y) are the polyline vertex coordinates; marker_dist is the desired
%   spacing.  Returns (marker_x, marker_y) with one point per spacing interval,
%   linearly interpolated along arc length.
    dist_from_start = cumsum( [0, sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2)] );marker_x=[];marker_y=[];
    marker_locs = marker_dist : marker_dist : dist_from_start(end);
    if length(marker_locs)>1
        marker_indices = interp1( dist_from_start, 1 : length(dist_from_start), marker_locs);
        marker_base_pos = floor(marker_indices);
        weight_second = marker_indices - marker_base_pos;
        marker_x = [marker_x, x(marker_base_pos) .* (1-weight_second) + x(marker_base_pos+1) .* weight_second];
        marker_y = [marker_y, y(marker_base_pos) .* (1-weight_second) + y(marker_base_pos+1) .* weight_second];
    end
end
