%% Reeb Path Planner (slave)
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [Path,wall_fol,adj]=ReebPath(adj,critP,reebEdge,cells,Start)
%
% Recursively traverses the Reeb graph from a given start node, ordering
% connected cells to minimise coverage overlap.
%
% INPUTS:
%   adj      - Adjacency matrix of the Reeb graph.
%   critP    - Critical point coordinates, one per row [x, y].
%   reebEdge - Reeb edge list (Nx2 index pairs).
%   cells    - Reeb cells (polyshape array indexed by reebCell).
%   Start    - Index of the starting critical point.
%
% OUTPUTS:
%   Path     - Ordered traversal: each row is [from, to, cellIndex].
%   wall_fol - Logical vector; 1 where the cell needs a wall-follow pass.
%   adj      - Updated adjacency matrix (visited edges zeroed).

    Path=[];wall_fol=[];
    A=adj;
    ind = Start;

    no=find(adj(ind,:));neighbor=[];
    repeat = adj(ind,logical(adj(ind,:)));
    for j= 1:length(no)
        neighbor = [neighbor;repmat(no(j),repeat(j),1)];
    end

    cellOrd=[];treebE=reebEdge;
    for j= 1:length(neighbor)
        [lia,lib]=ismember([ind,neighbor(j)],treebE,'rows');
        if lia
            treebE(lib,:)=[0,0];if lib>length(cells);lib=0;end
            cellOrd = [cellOrd;lib];
        end
        [lia,lib]=ismember([neighbor(j),ind],treebE,'rows');
        if lia
            treebE(lib,:)=[0,0];if lib>length(cells);lib=0;end
            cellOrd = [cellOrd;lib];
        end
    end

    neighborCon = sum(adj(neighbor,1:size(adj,2)~=ind),2);leftcell=zeros(length(neighbor(:,1)),1);
    leftcell(critP(neighbor,2) < critP(ind,2))=1;
    if cellOrd~=0; conn = [neighbor,neighborCon,leftcell,cellOrd,cells(cellOrd).area];
    else; conn = [neighbor,neighborCon,leftcell,cellOrd,0];end
    adj(ind,:)=0;adj(:,ind)=0;

    % Sort by connectivity count, then left cells first, then ascending area.
    for i=unique(conn(:,2))'
        temp2=conn(conn(:,2)==i,:);
        temp2=[sortrows(temp2(logical(temp2(:,3)),:),5,'ascend');sortrows(temp2(logical(~temp2(:,3)),:),5,'ascend')];
        if i>0
                Path=[Path;[repmat(ind,[length(temp2(:,1)),1]) ,temp2(:,[1,4])]];
                wall_fol = [wall_fol;temp2(:,2)==0];
                for j = temp2(1)'                                                       % Loops when two edges with connectivity > 0
                    Start = j;                                                          %|
                    [tPath,twall_fol,adj]=ReebPath(adj,critP,reebEdge,cells,Start);     %|Not
                    Path=[Path;tPath];                                                  %|Tested
                    wall_fol = [wall_fol;twall_fol];                                    %|
                end
        else
            Path=[Path;[repmat(ind,[length(temp2(:,1)),1]) ,temp2(:,[1,4])]];
            wall_fol = [wall_fol;temp2(:,2)==0];
        end
    end


    if sum(sum(adj(:,Path(end,2))))>0
        Start = Path(end,2);
        [tPath,twall_fol,adj]=ReebPath(adj,critP,reebEdge,cells,Start);            % Function Loops to find the path
        Path=[Path;tPath];
        wall_fol = [wall_fol;twall_fol];
    end

    if sum(sum(adj))==0;wall_fol(end)=0;end
end
