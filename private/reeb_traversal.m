%% Reeb Path Planner (master)
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [Path,wall_fol,adj]=reeb_traversal(adj,critP,reebEdge,cells,Start)
%
% Master Reeb-graph traversal. Calls ReebPath for each connected component
% and stitches disconnected components by jumping to the nearest unvisited
% critical point.
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

    Path=[];wall_fol=[];global spdist

    [tPath,twall_fol,adj]=ReebPath(adj,critP,reebEdge,cells,Start);
    Path=[Path;tPath];
    wall_fol = [wall_fol;twall_fol];wall_fol(end)=0;

    while sum(ismember(reebEdge,Path(:,1:2),'rows')|ismember(fliplr(reebEdge),Path(:,1:2),'rows'))~=size(reebEdge,1)

        Start = Path(end,2);
        t=find(~ismember(1:length(critP),reshape(Path(:,1:2),1,[])));

        [~,ind]=min(spdist(critP(Start,:),critP(t,:)));
        Path=[Path;[Start,t(ind),0]];
        wall_fol = [wall_fol;0];
        Start=t(ind);

        [tPath,twall_fol,adj]=ReebPath(adj,critP,reebEdge,cells,Start);            % Function Loops to find the path
        Path=[Path;tPath];
        wall_fol = [wall_fol;twall_fol];wall_fol(end)=0;

    end
end
