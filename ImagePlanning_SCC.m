%% Image Planning Function for SCC
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [node,edgeList,ttt,crackRaw]=ImagePlanning_SCC(img_n)
%
% Performs Image processing, extacts cracks and metadata, and computes the waypoints for the crackGraph using Visibility Graph.
%
% INPUTS:
%   img_n = RGB or Binary image of the workspace. 
%
% OUTPUTS:
%   node = [x0; y0] = nodes list of all edges in crackGraph in undirected graph format. 
%   edgeList = [node1; node2] = edges list of all edges in crackGraph in undirected graph format.
%   ttt = Computation time.

spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);
midP= @(P1,P2) (P1(:)+ P2(:)).'/2;
total_length = @(Ps) sum(sqrt(sum(diff(Ps).*diff(Ps),2)));

inpxMap = @(x) fix(((x)*25.4)/2);
pxinMap = @(x) round((x)*2/25.4,1);

botD = 48;                                       % 48" Diameter
footD= 7;                                        % 7" Diameter
sensD= 4.5*12;                                   % 4.5*12" Diameter

r1=inpxMap(botD/2);                      % Robot Radius 
a=inpxMap(footD/2);                      % Footprint Radius
s=inpxMap(sensD/2);                      % Sensor Range Radius

worki=false;
plt = false;
% img_n='myCrack2_45';

if ischar(img_n)
    % A Gaussian map name (contains 'Gauss') loads its .mat directly (crack pixels are
    % already foreground; no imbinarize/bwareaopen/inversion needed). Accepts an explicit
    % 'Gaussian<b>/myCrackGauss_s<sig>_<den>' path or a bare name (searched in CrackMaps/Gaussian1..8).
    if contains(img_n, 'Gauss')
        if contains(img_n, '/') || contains(img_n, '\')               % explicit Gaussian<b>/name
            matpath = fullfile('CrackMaps', [img_n '.mat']);
        else                                                          % bare name -> search Gaussian1..8
            matpath = '';
            for gb = 1:8
                cand = fullfile('CrackMaps', sprintf('Gaussian%d', gb), [img_n '.mat']);
                if exist(cand, 'file'); matpath = cand; break; end
            end
            if isempty(matpath)
                error('ImagePlanning_SCC:mapNotFound', ...
                    'Gaussian map "%s" not found in CrackMaps/Gaussian1..8', img_n);
            end
        end
        L = load(matpath);
        BW = L.crackGen; BW(end,end) = 0;                             % crack pixels are foreground
        rawBW = BW;
    else                                                              % Uniform .png path
    BW = imread(['CrackMaps/Uniform/' img_n '.png']);
    BW = imbinarize(BW);
    BW = BW(:,:,1);

    BW = bwareaopen(BW, 50);        % Remove isolated specks smaller than 50 px
    BW = ~BW;                       % Invert so cracks are foreground

    IM2 = imcomplement(BW);
    rawBW = BW;
    end

else
	BW = img_n;
end

BW2=bwmorph(BW,'fill');             % Fill isolated single-pixel holes
BW2 = bwskel(BW2>0);

% Remove short side-branches (arc-length <= a); preserves true crack endpoints.
BW3 = pruneBranches(BW2, a);


if worki; figure, imshow(~BW3); hold on; end
if plt; figure, imshow(~BW3); hold on; end
[rowBW, colBW]= size(BW3);
% pos=get(gcf, 'Position');

skelBW = BW3;

endPoints = bwmorph(BW3, 'endpoints');  % Skeleton endpoints
[endP_row, endP_col]=find(endPoints);
if worki; plot(endP_col,endP_row,'b*'); end

intPoints = bwmorph(BW3, 'branchpoints');
[intP_row, intP_col]=find(intPoints);
if worki; plot(intP_col,intP_row,'r*'); end

I=BW3;

endP=[endP_row endP_col];

dir_map=[-1 -1;-1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0 ; 1 1];
dir_can=[1 2 4 3 6;2 1 3 4 5; 3 2 5 1 8; 4 1 6 2 7 ; 5 3 8 2 7; 6 4 7 1 8; 7 6 8 4 5; 8 5 7 3 6];
colors={'y','m','c','r','g','b'};

% Build the crack graph: trace all skeleton segments (including isolated cycles),
% fit a taut-string minimum-cover polyline to each segment, then merge endpoints
% and junctions into an undirected graph (node/edgeList).
tfun=tic;
crackRaw = trace_all_segments(BW3, dir_map);               % trace all skeleton segments + cycles
[node, edgeList] = taut_graph(crackRaw, a, 0.9);           % per-crack taut-string -> nodes+edges
ttt = toc(tfun);
end
%% Functions

function ang = ab2v(a,b) 
    theta = rad2deg(atan2(norm(cross([1,0,0],[a,0])), dot([1,0,0],[a,0])));
    if a(2)<0 theta=360-theta; end
    
    aR = a*[1;1i]*exp(-1i*theta*pi/180); aR=[real(aR) imag(aR)];
    bR = b*[1;1i]*exp(-1i*theta*pi/180); bR=[real(bR) imag(bR)];

    ang=rad2deg(atan2(norm(cross([aR,0],[bR,0])), dot([aR,0],[bR,0])));
    if bR(2)<0 ang=360-ang; end
end


function y = bound(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
end


function BW3 = pruneBranches(BW, a)
% Remove skeleton side-branches of arc-length <= a (short outlier stubs sticking off
% a larger crack; the robot footprint sweeps over them while covering the main crack,
% so they need no separate fill). A branch is an endpoint->branchpoint segment.
% A through-crack's true endpoint (segment ending at another endpoint, not a
% branchpoint) is never removed.
    BW3 = BW > 0;
    while true
        B = bwmorph(BW3,'branchpoints');
        E = bwmorph(BW3,'endpoints');
        if ~any(B(:)) || ~any(E(:)); break; end
        adjB = imdilate(B, true(3)) & BW3 & ~B;     % skeleton pixels adjacent to a branchpoint
        seg  = BW3 & ~B;                            % break the skeleton at branchpoints
        CC = bwconncomp(seg, 8);
        rm = false(size(BW3)); anyRm = false;
        for k = 1:CC.NumObjects
            px = CC.PixelIdxList{k};
            if any(E(px)) && any(adjB(px))          % side-branch: has an endpoint AND touched a branchpoint
                segM = false(size(BW3)); segM(px) = true;
                ep1  = find(bwmorph(segM,'endpoints'), 1);
                if isempty(ep1)
                    L = numel(px);
                else
                    D = bwdistgeodesic(segM, ep1, 'quasi-euclidean'); L = max(D(isfinite(D)));
                end
                if L <= a; rm(px) = true; anyRm = true; end
            end
        end
        if ~anyRm; break; end
        BW3(rm) = false;
    end
end


function crackRaw = trace_all_segments(BW3, dir_map)
% compCrack_branching (endpoint+branchpoint-seeded) + a completeness guard: any UNTRACED skeleton
% (cycles/loops, which have no endpoint/branchpoint so compCrack_branching skips them) is walked
% directly -> every crack pixel enters the graph (no dropped crack).
    BW3 = BW3 > 0;
    [sr, sc] = find(BW3);
    skpx = [sr, sc];                                   % all skeleton pixels (row,col)
    E = bwmorph(BW3, 'endpoints');
    B = bwmorph(BW3, 'branchpoints');
    [er, ec] = find(E | B);
    seed = [er, ec];                                   % endpoint+branchpoint seed
    crackRaw0 = compCrack_branching(BW3, seed, dir_map, {});
    crackRaw = crackRaw0(~cellfun(@isempty, crackRaw0));
    crackRaw = crackRaw(:)';
    if isempty(crackRaw); return; end
    traced = cell2mat(crackRaw(:));                    % all traced pixels (row,col)
    [~, dd] = dsearchn(traced, skpx);
    residual = skpx(dd > 1.5, :);                      % skeleton px not within 1.5 of any traced pt
    if ~isempty(residual)
        rmask = false(size(BW3));
        rmask(sub2ind(size(BW3), residual(:,1), residual(:,2))) = true;
        CC = bwconncomp(rmask, 8);
        for k = 1:CC.NumObjects
            [rr, cc] = ind2sub(size(BW3), CC.PixelIdxList{k});
            ls = walk_component([rr, cc]);
            crackRaw = [crackRaw, ls];                 %#ok<AGROW>
        end
    end
end


function [node, edgeList] = taut_graph(crackRaw, a, cover_frac)
% Each crack segment's centerline -> taut-string boundary min-cover polyline (covered
% within cover_frac*a). Its vertices become graph nodes, consecutive vertices become
% edges (downstream straight-chord emission reproduces the taut polyline). Endpoint/
% junction vertices shared within merge_tol are merged so crossings connect.
%   node     : Nx2 (row,col) ;  edgeList : Ex2 1-based undirected node-index pairs.
    merge_tol = 6.0;
    reps  = zeros(0,2);
    edges = zeros(0,2);
    for ci = 1:numel(crackRaw)
        c = crackRaw{ci};
        if size(c,1) < 2; continue; end
        v = min_cover_path_adaptive(c, a*cover_frac);
        ids = zeros(size(v,1),1);
        for p = 1:size(v,1)
            [reps, ids(p)] = nid(reps, v(p,:), merge_tol);
        end
        for k = 1:numel(ids)-1
            if ids(k) ~= ids(k+1)
                edges(end+1,:) = [ids(k), ids(k+1)];   %#ok<AGROW>
            end
        end
    end
    node = reps;
    if isempty(edges); edgeList = zeros(0,2); return; end
    key = sort(edges, 2);                              % undirected dedupe (keep first occurrence)
    [~, ia] = unique(key, 'rows', 'stable');
    edgeList = edges(ia, :);
end


function [reps, id] = nid(reps, p, merge_tol)
% Nearest-rep node id (1-based); merge a new vertex into an existing node within merge_tol.
    if ~isempty(reps)
        d2 = (reps(:,1)-p(1)).^2 + (reps(:,2)-p(2)).^2;
        [m, k] = min(d2);
        if m <= merge_tol^2
            id = k; return;
        end
    end
    reps(end+1,:) = p;
    id = size(reps,1);
end


function best = min_cover_path_adaptive(centerline, a)
% Per-crack OPTIMAL taut-string: binary-search the largest disk_frac whose path still
% has max-gap <= a (coverage guaranteed), recovering cost the fixed margin leaves behind.
% max-gap is monotone increasing in disk_frac, so bisection is exact to 2^-bisect.
    lo = 0.5; hi = 1.0; bisect = 6;
    cl = resample_poly(centerline, 3.0);
    best = min_cover_path(centerline, a, lo);          % lo assumed feasible
    if gapf(cl, best) > a                              % even lo violates -> fall back to centerline
        best = dp_simplify(centerline, 0.9*a); return;
    end
    L = lo; H = hi;
    for it = 1:bisect
        m = 0.5*(L + H);
        P = min_cover_path(centerline, a, m);
        if gapf(cl, P) <= a
            L = m; best = P;
        else
            H = m;
        end
    end
end


function g = gapf(cl, P)
% max distance of any centerline sample to the path polyline (coverage check).
    g = 0;
    for i = 1:size(cl,1)
        d = point_to_polyline_dist(cl(i,:), P);
        if d > g; g = d; end
    end
end


function out = min_cover_path(centerline, a, disk_frac)
% Deterministic taut-string (elastic band) through ordered radius-(disk_frac*a) disks on
% the crack centerline: resample -> Gauss-Seidel straighten + project into disk -> DP.
    iters = 400; simplify_tol = 1.0;
    r = disk_frac * a;                                 % disk radius < a: safety margin
    step = a / 3.0;
    c = resample_poly(centerline, step);               % ordered disk centers
    if size(c,1) < 3; out = c; return; end
    P = c;                                             % init on the crack (feasible)
    for it = 1:iters
        maxmove = 0; Pn = P;
        for i = 2:size(P,1)-1
            mid = 0.5*(Pn(i-1,:) + P(i+1,:));          % Gauss-Seidel straighten
            v = mid - c(i,:); d = hypot(v(1), v(2));
            if d <= r
                new = mid;                             % project into D(c_i, r)
            else
                new = c(i,:) + v/d*r;
            end
            mv = hypot(new(1)-Pn(i,1), new(2)-Pn(i,2));
            if mv > maxmove; maxmove = mv; end
            Pn(i,:) = new;
        end
        P = Pn;
        if maxmove < 0.05; break; end
    end
    out = dp_simplify(P, simplify_tol);
end


function out = resample_poly(poly, step)
% Resample a polyline at uniform arc-length spacing `step` (endpoint always included).
    if size(poly,1) < 2; out = poly; return; end
    seg = hypot(diff(poly(:,1)), diff(poly(:,2)));
    s = [0; cumsum(seg)];
    if s(end) < step; out = poly([1 end],:); return; end
    [s, ui] = unique(s, 'stable');        % drop zero-length-segment duplicate knots
    poly = poly(ui,:);                     %   (np.interp tolerates them; interp1 requires unique)
    if numel(s) < 2; out = poly([1 end],:); return; end
    u = (0:step:s(end))';
    if abs(u(end) - s(end)) > 1e-9; u = [u; s(end)]; end
    x = interp1(s, poly(:,1), u);
    y = interp1(s, poly(:,2), u);
    out = [x, y];
end


function out = dp_simplify(pts, tol)
% Iterative Douglas-Peucker polyline simplification (point-to-chord deviation > tol kept).
    n = size(pts,1);
    if n < 3; out = pts; return; end
    keep = false(n,1); keep(1) = true; keep(n) = true;
    stack = [1, n];
    while ~isempty(stack)
        s = stack(end,1); e = stack(end,2); stack(end,:) = [];
        if e <= s+1; continue; end
        ab = pts(e,:) - pts(s,:); L = hypot(ab(1), ab(2));
        idx = (s+1):(e-1);
        if L < 1e-9
            d = hypot(pts(idx,1)-pts(s,1), pts(idx,2)-pts(s,2));
        else
            d = abs(ab(1)*(pts(idx,2)-pts(s,2)) - ab(2)*(pts(idx,1)-pts(s,1))) / L;
        end
        [dm, ki] = max(d); k = idx(ki);
        if dm > tol
            keep(k) = true;
            stack(end+1,:) = [s, k];                   %#ok<AGROW>
            stack(end+1,:) = [k, e];                   %#ok<AGROW>
        end
    end
    out = pts(keep,:);
end


function lines = walk_component(P)
% Trace a connected skeleton component (Mx2 (row,col), 8-adjacency) into polylines.
% Handles CYCLES (no endpoint/branchpoint -> compCrack_branching skips them) so every pixel routes.
    lines = {};
    M = size(P,1);
    if M < 2; return; end
    nb = cell(M,1); deg = zeros(M,1);
    for i = 1:M
        d = max(abs(P(:,1)-P(i,1)), abs(P(:,2)-P(i,2)));   % chebyshev (8-conn)
        idx = find(d == 1);
        nb{i} = idx(:)'; deg(i) = numel(idx);
    end
    used = containers.Map('KeyType','char','ValueType','logical');
    order = [find(deg == 1); find(deg ~= 1)];              % deg-1 starts first, then the rest
    for oi = 1:numel(order)
        s0 = order(oi);
        allused = true;
        for q = nb{s0}
            if ~isKey(used, ekey(s0,q)); allused = false; break; end
        end
        if allused; continue; end
        cur = s0; prev = 0; lineIdx = s0;
        while true
            cand = [];
            for q = nb{cur}
                if q ~= prev && ~isKey(used, ekey(cur,q)); cand(end+1) = q; end %#ok<AGROW>
            end
            if isempty(cand); break; end
            [~, mi] = sortrows(P(cand,:));                 % min by (row,col) lexicographic
            pick = cand(mi(1));
            used(ekey(cur,pick)) = true;
            lineIdx(end+1) = pick;                         %#ok<AGROW>
            prev = cur; cur = pick;
            if cur == s0; break; end
        end
        if numel(lineIdx) >= 2
            lines{end+1} = P(lineIdx,:);                   %#ok<AGROW>
        end
    end
end


function k = ekey(u, v)
% Undirected edge key for walk_component's used-edge set.
    k = sprintf('%d_%d', min(u,v), max(u,v));
end


function dmin = point_to_polyline_dist(p, P)
% Minimum distance from point p to polyline P (min over its segments).
    if size(P,1) == 1
        dmin = hypot(p(1)-P(1,1), p(2)-P(1,2)); return;
    end
    A = P(1:end-1,:); B = P(2:end,:);
    AB = B - A; AP = [p(1)-A(:,1), p(2)-A(:,2)];
    denom = sum(AB.^2, 2);
    t = sum(AP.*AB, 2) ./ denom;
    t(denom == 0) = 0;
    t = min(max(t, 0), 1);
    proj = A + [t.*AB(:,1), t.*AB(:,2)];
    d = hypot(proj(:,1)-p(1), proj(:,2)-p(2));
    dmin = min(d);
end