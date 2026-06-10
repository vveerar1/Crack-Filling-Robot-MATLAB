% celled_rpp  Geometry-routed cell Rural-Postman coverage path.
%
% Builds a minimum-cost traversal over the MCD coverage cells and crack graph:
% min-weight T-join (intlinprog) to eulerianise odd-degree nodes, MST-based
% component reconnection (avoiding attachment at interior crack nodes), then an
% Euler circuit (Hierholzer) assembled with a DP sweep-direction chooser and
% wall clamping.
%
%   PE = celled_rpp(node, crackEdge, cells, s, a)
%     node     : (N,2) crack-graph nodes (row,col)
%     crackEdge: (E,2) crack edges, 1-based node indices
%     cells    : polyshape array (MCD coverage cells, splitReg_work)
%     s, a     : sensor / footprint radii (px)
%   PE         : (M,2) coverage path (x,y)
function PE = celled_rpp(node, crackEdge, cells, s, a)
    SP = s/sqrt(2);
    even_xy = even_crack_nodes(node, crackEdge);
    best_len = inf; PE = zeros(0,2);
    for corner = [false true]
        try
            [coords, reqE, reqW, reqCell, sweepPairs] = build_graph(node, crackEdge, cells, s, a, corner);
            euler = solve_rpp(coords, reqE, reqW, reqCell, sweepPairs);
        catch ME
            if ~isempty(getenv('CRPP_DBG')); fprintf(2,'[solve corner=%d] %s\n', corner, getReport(ME)); end
            continue
        end
        try
            [P, transit] = assemble(euler, cells, s, a);
            if has_midcrack(transit, even_xy)
                P = assemble_euler(euler, cells, s, a);
            end
            if size(P,1) >= 2
                P = wall_clamp(P, node, crackEdge, s);
                L = plen(P);
                if L < best_len; best_len = L; PE = P; end
            end
        catch ME
            if ~isempty(getenv('CRPP_DBG')); fprintf(2,'[assemble corner=%d] %s\n', corner, getReport(ME)); end
            continue
        end
    end
end

% --------------------------------------------------------------------------- %
function L = plen(P)
    if size(P,1) < 2; L = 0; else; L = sum(sqrt(sum(diff(P,1,1).^2,2))); end
end

function lp = loops_largest(V)
    % largest (outer) loop of a NaN-separated (x,y) vertex array
    V = reshape(V, [], 2);
    nanrow = isnan(V(:,1));
    starts = [1; find(nanrow)+1]; ends = [find(nanrow)-1; size(V,1)];
    bestA = -1; lp = zeros(0,2);
    for k = 1:numel(starts)
        if starts(k) > ends(k); continue; end
        seg = V(starts(k):ends(k), :); seg = seg(~isnan(seg(:,1)), :);
        if size(seg,1) < 3; continue; end
        A = abs(0.5*sum(seg(:,1).*seg([2:end 1],2) - seg([2:end 1],1).*seg(:,2)));
        if A > bestA; bestA = A; lp = seg; end
    end
end

function A = polyArea(p)
    A = abs(0.5*sum(p(:,1).*p([2:end 1],2) - p([2:end 1],1).*p(:,2)));
end

% --------------------------------------------------------------------------- %
function opts = cell_options(cellPoly, s_in, a_in)
    % 4 boustrophedon orientations -> struct array (entry, exit, poly, len).
    % BoustrophedonPath reads the globals s, allNode, a by those exact names.
    global s allNode a %#ok<GVMIS>
    s = s_in; a = a_in;
    lp = loops_largest(cellPoly.Vertices);
    [~,iL] = min(lp(:,1)); [~,iR] = max(lp(:,1));
    pL = lp(iL,:); pR = lp(iR,:);
    bp = s_in/sqrt(2);
    opts = struct('entry',{},'exit',{},'poly',{},'len',{});
    for startLeft = [true false]
        if startLeft; a0 = pL; a1 = pR; else; a0 = pR; a1 = pL; end
        for d = [0 1]
            allNode = [a0(2) a0(1); a1(2) a1(1)];      % (row,col); reebEdge=[1 2]
            subXY = BoustrophedonPath(cellPoly, cellPoly, [1 2], bp, d, 0, 0, true, false);
            subXY = reshape(subXY, [], 2);
            if ~isempty(subXY)
                opts(end+1) = struct('entry',subXY(1,:),'exit',subXY(end,:),...
                                     'poly',subXY,'len',plen(subXY)); %#ok<AGROW>
            end
        end
    end
    if isempty(opts)
        opts(1) = struct('entry',pL,'exit',pR,'poly',[pL;pR],'len',plen([pL;pR]));
    end
end

function [e0,x0] = corner_ports(cellPoly, s, a)
    opts = cell_options(cellPoly, s, a);
    keys = round([vertcat(opts.entry)], 1);
    [~,ix] = sortrows(keys);            % lexicographically smallest entry
    e0 = opts(ix(1)).entry; x0 = opts(ix(1)).exit;
end

% --------------------------------------------------------------------------- %
function [coords, reqE, reqW, reqCell, sweepPairs] = build_graph(node, crackEdge, cells, s, a, corner)
    nfx = fliplr(node);                 % (row,col)->(x,y)
    coords = nfx;
    nC = size(crackEdge,1);
    reqE = crackEdge;                   % crack edges first (1-based into coords)
    reqW = sqrt(sum((nfx(crackEdge(:,1),:)-nfx(crackEdge(:,2),:)).^2,2));
    reqCell = zeros(nC,1);              % 0 = crack edge
    sweepPairs = zeros(numel(cells),2);
    for c = 1:numel(cells)
        if corner
            [pl,pr] = corner_ports(cells(c), s, a);
        else
            lp = loops_largest(cells(c).Vertices);
            [~,iL]=min(lp(:,1)); [~,iR]=max(lp(:,1)); pl=lp(iL,:); pr=lp(iR,:);
        end
        w = polyArea(loops_largest(cells(c).Vertices)) / (s/sqrt(2));
        iL = size(coords,1)+1; coords(iL,:) = pl;
        iR = size(coords,1)+1; coords(iR,:) = pr;
        reqE(end+1,:) = [iL iR]; reqW(end+1) = w; reqCell(end+1) = c; %#ok<AGROW>
        sweepPairs(c,:) = [iL iR];
    end
end

% --------------------------------------------------------------------------- %
function euler = solve_rpp(coords, reqE, reqW, reqCell, sweepPairs)
    N = size(coords,1);
    [~,order] = sortrows(round(coords,3));      % canonical lexicographic relabel
    cid = zeros(N,1); cid(order) = (1:N)';      % cid(orig) = canonical id
    inv = order;                                % inv(canonical) = orig
    C = @(i) cid(i);

    % edge list in canonical ids: [u v], weight, cell, required
    E = [C(reqE(:,1)) C(reqE(:,2))];
    W = reqW(:); ECell = reqCell(:); EReq = true(size(E,1),1);

    % crack degree (canonical) for mid-crack avoidance
    crackDeg = zeros(N,1);
    for e = 1:size(E,1)
        if ECell(e)==0
            crackDeg(E(e,1)) = crackDeg(E(e,1))+1;
            crackDeg(E(e,2)) = crackDeg(E(e,2))+1;
        end
    end
    avoid = find(mod(crackDeg,2)==0 & crackDeg>0);

    dcoord = @(a,b) hypot(coords(inv(a),1)-coords(inv(b),1), coords(inv(a),2)-coords(inv(b),2));

    forbid = sort([C(sweepPairs(:,1)) C(sweepPairs(:,2))],2);   % don't match a cell's 2 ports

    % node degree
    deg = accumarray([E(:,1);E(:,2)], 1, [N 1]);
    odd = sort(find(mod(deg,2)==1));
    if ~isempty(odd)
        % shortest-path distances in the required-edge graph
        Greq = graph(E(:,1), E(:,2), W, N);
        D = distances(Greq);
        big = 0;
        for i=1:numel(odd); for j=i+1:numel(odd); big=max(big,dcoord(odd(i),odd(j))); end; end
        big = big + 1;
        % candidate pairs (not forbidden)
        pairs = []; cost = [];
        for i=1:numel(odd)
            for j=i+1:numel(odd)
                x=odd(i); y=odd(j); lo=min(x,y); hi=max(x,y);
                if any(forbid(:,1)==lo & forbid(:,2)==hi); continue; end
                e = dcoord(x,y); gd = D(x,y);
                pen = 0; if isfinite(gd) && e >= 0.85*gd; pen = big; end
                pairs(end+1,:) = [x y]; cost(end+1,1) = e + pen; %#ok<AGROW>
            end
        end
        % min-weight perfect matching via intlinprog
        np = size(pairs,1);
        Aeq = zeros(numel(odd), np);
        for p=1:np
            Aeq(odd==pairs(p,1), p) = 1; Aeq(odd==pairs(p,2), p) = 1;
        end
        beq = ones(numel(odd),1);
        opt = optimoptions('intlinprog','Display','off');
        x = intlinprog(cost, 1:np, [], [], Aeq, beq, zeros(np,1), ones(np,1), opt);
        sel = find(x > 0.5);
        msel = sortrows(sort(pairs(sel,:),2));    % deterministic order
        for k=1:size(msel,1)
            a=msel(k,1); b=msel(k,2);
            E(end+1,:)=[a b]; W(end+1)=dcoord(a,b); ECell(end+1)=0; EReq(end+1)=false; %#ok<AGROW>
        end
    end

    % connect remaining components (MST), avoiding mid-crack attach
    comp = conncomp(graph(E(:,1), E(:,2), [], N));
    ncomp = max(comp);
    if ncomp > 1
        cs = []; cd = []; cw = []; cpa = []; cpb = [];
        for i=1:ncomp
            for j=i+1:ncomp
                ni = find(comp==i); nj = find(comp==j);
                ci = setdiff(ni, avoid); if isempty(ci); ci = ni; end
                cj = setdiff(nj, avoid); if isempty(cj); cj = nj; end
                bestd = inf; ba=ci(1); bb=cj(1);
                for u=ci(:)'; for v=cj(:)'
                    dd=dcoord(u,v); if dd<bestd; bestd=dd; ba=u; bb=v; end
                end; end
                cs(end+1)=i; cd(end+1)=j; cw(end+1)=bestd; cpa(end+1)=ba; cpb(end+1)=bb; %#ok<AGROW>
            end
        end
        Gc = graph(cs, cd, cw, ncomp);
        T = minspantree(Gc);
        for te = 1:numedges(T)
            ij = sort(T.Edges.EndNodes(te,:));
            idx = find((cs==ij(1) & cd==ij(2)) | (cs==ij(2) & cd==ij(1)), 1);
            a=cpa(idx); b=cpb(idx);
            for rep=1:2     % add edge twice so all node degrees stay even
                E(end+1,:)=[a b]; W(end+1)=dcoord(a,b); ECell(end+1)=0; EReq(end+1)=false; %#ok<AGROW>
            end
        end
    end

    % Euler circuit (Hierholzer) from the canonical-smallest node
    eord = euler_edges(N, E, min(E(:)));
    euler = struct('typ',{},'cell',{},'p0',{},'p1',{});
    for t = 1:size(eord,1)
        e = eord(t,1); fromN = eord(t,2);
        u = fromN; v = E(e,1); if v==u; v=E(e,2); end
        p0 = coords(inv(u),:); p1 = coords(inv(v),:);
        if ECell(e) > 0
            typ = 'sweep'; cc = ECell(e);
        elseif EReq(e)
            typ = 'crack'; cc = 0;
        else
            typ = 'transit'; cc = 0;
        end
        euler(end+1) = struct('typ',typ,'cell',cc,'p0',p0,'p1',p1); %#ok<AGROW>
    end
end

function eord = euler_edges(N, E, startNode)
    % Hierholzer -> (edgeIdx, fromNode) in traversal order.
    M = size(E,1);
    adj = cell(N,1);
    for e=1:M; adj{E(e,1)}(end+1)=e; adj{E(e,2)}(end+1)=e; end
    used = false(M,1); ptr = ones(N,1);
    stackV = startNode; stackE = 0; stackF = 0;     % node, edge-to-here, from-node
    circ = zeros(0,2);
    while ~isempty(stackV)
        v = stackV(end);
        while ptr(v) <= numel(adj{v}) && used(adj{v}(ptr(v))); ptr(v)=ptr(v)+1; end
        if ptr(v) > numel(adj{v})
            if stackE(end) > 0; circ(end+1,:) = [stackE(end) stackF(end)]; end %#ok<AGROW>
            stackV(end)=[]; stackE(end)=[]; stackF(end)=[];
        else
            e = adj{v}(ptr(v)); used(e)=true;
            w = E(e,1); if w==v; w=E(e,2); end
            stackV(end+1)=w; stackE(end+1)=e; stackF(end+1)=v; %#ok<AGROW>
        end
    end
    eord = flipud(circ);
end

% --------------------------------------------------------------------------- %
function [PE, transit] = assemble(euler, cells, s, a)
    % DP over ordered required edges: choose cell sweep orientation and crack direction.
    elems = {}; cache = containers.Map('KeyType','double','ValueType','any');
    for k = 1:numel(euler)
        ev = euler(k);
        if strcmp(ev.typ,'crack')
            L = hypot(ev.p1(1)-ev.p0(1), ev.p1(2)-ev.p0(2));
            o(1) = mkopt(ev.p0, ev.p1, [ev.p0;ev.p1], L);
            o(2) = mkopt(ev.p1, ev.p0, [ev.p1;ev.p0], L);
            elems{end+1} = o; clear o %#ok<AGROW>
        elseif strcmp(ev.typ,'sweep')
            if ~isKey(cache, ev.cell); cache(ev.cell) = cell_options(cells(ev.cell), s, a); end
            opts = cache(ev.cell); o = struct('entry',{},'exit',{},'poly',{},'len',{});
            for q=1:numel(opts); o(q)=mkopt(opts(q).entry,opts(q).exit,opts(q).poly,opts(q).len); end
            elems{end+1} = o; clear o %#ok<AGROW>
        end
    end
    transit = {};
    if isempty(elems); PE = zeros(0,2); return; end
    Nn = numel(elems);
    dp = arrayfun(@(o) o.len, elems{1}); back = {-ones(1,numel(elems{1}))};
    for i=2:Nn
        cur = elems{i}; ndp = inf(1,numel(cur)); bi = -ones(1,numel(cur));
        for oi=1:numel(cur)
            for j=1:numel(elems{i-1})
                pj = elems{i-1}(j);
                c = dp(j) + hypot(cur(oi).entry(1)-pj.exit(1), cur(oi).entry(2)-pj.exit(2)) + cur(oi).len;
                if c < ndp(oi); ndp(oi)=c; bi(oi)=j; end
            end
        end
        dp = ndp; back{end+1}=bi; %#ok<AGROW>
    end
    chosen = zeros(1,Nn); [~,chosen(Nn)] = min(dp);
    for i=Nn:-1:2; chosen(i-1)=back{i}(chosen(i)); end
    segs = {}; prev = [];
    for i=1:Nn
        o = elems{i}(chosen(i));
        if ~isempty(prev)
            segs{end+1} = [prev; o.entry]; transit{end+1} = [prev; o.entry]; %#ok<AGROW>
        end
        segs{end+1} = o.poly; prev = o.exit; %#ok<AGROW>
    end
    PE = vertcat(segs{:});
end

function o = mkopt(entry, ex, poly, len)
    o = struct('entry',entry,'exit',ex,'poly',poly,'len',len);
end

function PE = assemble_euler(euler, cells, s, a)
    cache = containers.Map('KeyType','double','ValueType','any'); segs = {};
    for k=1:numel(euler)
        ev = euler(k);
        if strcmp(ev.typ,'sweep')
            if ~isKey(cache, ev.cell); cache(ev.cell) = cell_options(cells(ev.cell), s, a); end
            opts = cache(ev.cell);
            df = inf; dr = inf; pf=[]; pr=[];
            for q=1:numel(opts)
                f = hypot(opts(q).entry(1)-ev.p0(1),opts(q).entry(2)-ev.p0(2)) + hypot(opts(q).exit(1)-ev.p1(1),opts(q).exit(2)-ev.p1(2));
                r = hypot(opts(q).entry(1)-ev.p1(1),opts(q).entry(2)-ev.p1(2)) + hypot(opts(q).exit(1)-ev.p0(1),opts(q).exit(2)-ev.p0(2));
                if f<df; df=f; pf=opts(q).poly; end
                if r<dr; dr=r; pr=flipud(opts(q).poly); end
            end
            if dr < df; segs{end+1}=pr; else; segs{end+1}=pf; end %#ok<AGROW>
        else
            segs{end+1} = [ev.p0; ev.p1]; %#ok<AGROW>
        end
    end
    PE = []; prev=[];
    for k=1:numel(segs)
        poly = segs{k};
        if ~isempty(prev) && hypot(poly(1,1)-prev(1),poly(1,2)-prev(2)) > 1e-6
            PE = [PE; prev; poly(1,:)]; %#ok<AGROW>
        end
        PE = [PE; poly]; prev = poly(end,:); %#ok<AGROW>
    end
    if isempty(PE); PE = zeros(0,2); end
end

% --------------------------------------------------------------------------- %
function ev = even_crack_nodes(node, crackEdge)
    deg = accumarray([crackEdge(:,1);crackEdge(:,2)], 1, [size(node,1) 1]);
    nfx = fliplr(node);
    ev = nfx(mod(deg,2)==0 & deg>0, :);
end

function tf = has_midcrack(transit, even_xy)
    tf = false;
    if isempty(even_xy) || isempty(transit); return; end
    for k=1:numel(transit)
        seg = transit{k};
        if hypot(seg(2,1)-seg(1,1), seg(2,2)-seg(1,2)) < 15; continue; end
        for r=1:2
            if any(hypot(even_xy(:,1)-seg(r,1), even_xy(:,2)-seg(r,2)) <= 6); tf=true; return; end
        end
    end
end

function out = wall_clamp(PE, node, crackEdge, s)
    SP = s/sqrt(2); W = 3048; H = 2896;
    nfx = fliplr(node); cp = nfx;
    for e=1:size(crackEdge,1)
        A = nfx(crackEdge(e,1),:); B = nfx(crackEdge(e,2),:);
        n = max(2, floor(hypot(B(1)-A(1),B(2)-A(2))/40));
        cp = [cp; [linspace(A(1),B(1),n)' linspace(A(2),B(2),n)']]; %#ok<AGROW>
    end
    out = PE; cx = cp(:,1); cy = cp(:,2);
    for i=1:size(out,1)
        x = out(i,1); y = out(i,2);
        if y < SP        && ~any(cy < SP    & abs(cx-x) < SP); out(i,2) = SP;     end
        if y > H-SP      && ~any(cy > H-SP  & abs(cx-x) < SP); out(i,2) = H-SP;   end
        if x < SP        && ~any(cx < SP    & abs(cy-y) < SP); out(i,1) = SP;     end
        if x > W-SP      && ~any(cx > W-SP  & abs(cy-y) < SP); out(i,1) = W-SP;   end
    end
end
