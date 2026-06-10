%% Export SCC results to CSV
% Run this script after SCC.m to save workspace variables as CSV files.
% Output files are written to Results/matlab_export/.

out_dir = 'Results/matlab_export';
if ~exist(out_dir,'dir'); mkdir(out_dir); end

% ── PathEdge (primary output) ────────────────────────────────────────────────
writematrix(PathEdge, [out_dir '/PathEdge.csv']);
fprintf('Saved PathEdge: %d x %d\n', size(PathEdge,1), size(PathEdge,2));

% ── Intermediate values ──────────────────────────────────────────────────────
writematrix(node,      [out_dir '/node.csv']);
writematrix(edgeList,  [out_dir '/edgeList.csv']);
writematrix(subcritP,  [out_dir '/subcritP.csv']);
writematrix(reebEdge,  [out_dir '/reebEdge.csv']);
writematrix(allNode,   [out_dir '/allNode.csv']);

% ── Metrics ──────────────────────────────────────────────────────────────────
writematrix(res,       [out_dir '/res.csv']);

fprintf('All exports complete → %s/\n', out_dir);
