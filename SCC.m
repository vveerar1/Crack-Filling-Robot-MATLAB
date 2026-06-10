%% SCC — Sensor-based Complete Coverage with known crack locations
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

clc
close all
clearvars -except k

ress=[];
for sa=1%:10
clearvars -except sa ress

global reebEdge reebCell crackEdge allNode vertical spdist spdist2 s a total_length
inpxMap = @(x) fix(((x)*25.4)/2);
pxinMap = @(x) round((x)*2/25.4,1);
mmpxMap = @(x) fix(((x)/2));
pxmmMap = @(x) round((x)*2,1);

warning off

GCC=0;

sim =false;

%% Visualization
animate     = false;   % show the live animation figure (plot even without writing frames/GIF)
makeGif     = false;   % stitch frames into Results/GIF/<map>.gif  (also triggers the figure)
writeFrames = false;   % save per-step PNG frames to Results/GIF/  (also triggers the figure)
fps         = 14;      % GIF frame rate
animate = animate || writeFrames || makeGif;   % figure is created when animate is set OR frames/GIF are wanted

% route: 'rpp' = geometry-routed cell Rural-Postman; 'native' = reeb-graph Chinese-Postman
route = 'rpp';

% Robot dimensions loaded from robot_config.json (edit that file to reconfigure the robot)
cfg   = config_loader('robot1');
botD  = cfg.base_diameter_in; footD = cfg.footprint_diameter_in; sensD = cfg.sensor_diameter_in;
r1    = cfg.r1;                         % Robot Radius
a     = cfg.a;                          % Footprint Radius
s     = cfg.s;                          % Sensor Range Radius

vertical = @(P) all(P(1,1)==P(:,1));
spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);
spdist2 = @(Ps1,Ps2) sqrt((Ps1(:,1)-Ps2(:,1)).^2 + (Ps1(:,2)-Ps2(:,2)).^2);
total_length = @(Ps) sum(sqrt(sum(diff(Ps).*diff(Ps),2)));


%% Importing Map

den = [35,45,50,65,80,90,95,100];   % density values
sig = [5,10,20];                    % Gaussian sigmas
Gau = 0;                            % 0 = Uniform, 1 = Gaussian
k   = 8;                            % density index into den
mapN = 1;                           % map number
if ~Gau
    img_n = ['myCrack' num2str(k) '_' num2str(den(k)) '_' num2str(mapN)];
else
    Gaussb = 6;                     % Gaussian folder 1..8
    b      = 3;                     % sigma index into sig (1..3)
    img_n  = ['Gaussian' num2str(Gaussb) '/myCrackGauss_s' num2str(sig(b)) '_' num2str(den(k))];
end
vizName = strrep(img_n,'/','_');
gifPath = ['Results/GIF/' vizName '.gif'];
if ~exist('Results/GIF','dir'); mkdir('Results/GIF'); end
delete(['Results/GIF/' vizName '_*.png']);
if exist(gifPath,'file'); delete(gifPath); end
[node,edgeList,ttt,crackRaw_c] = ImagePlanning_SCC(img_n);   % handles Uniform .png + Gaussian .mat by name
rowBW = 2896; colBW = 3048;
crackGen = zeros(rowBW,colBW);

%%

smflag=0;
ppn=0;iMor=0;

BW3=crackGen;
[rowBW, colBW]= size(BW3);
% Tight axis limits (image extent, no inner padding); even space is added OUTSIDE the axis box
xPad  = [0, colBW];
yPad  = [0, rowBW];
tol =0.00001;
colors={'y','m','c','r','g','b'};
if animate
    fig1 = figure('Color','w'); imshow(~crackGen); hold on
    set(gca,'Color','w','XColor','k','YColor','k');
    set(fig1,'Position',[100 100 900 900]);
    set(gca,'FontSize',20,'FontWeight','bold');
    xlabel('x (m)');ylabel('y (m)')
    axis([0-150 colBW+150 0-150 rowBW+150])
    axis on
    pbaspect([1 1 1])
    xlab=get(gca,'xtickLabel');ylab=get(gca,'ytickLabel');
    xlab={};for ppx=round(pxmmMap(linspace(0,3050,7)/1000));xlab=[xlab,num2str(ppx)];end
    ylab={};for ppx=round(pxmmMap(linspace(0,2898,7)/1000));ylab=[ylab,num2str(ppx)];end
    xticks(mmpxMap(round(pxmmMap(linspace(0,3050,7)/1000))*1000));set(gca,'xtickLabel',xlab);yticks(mmpxMap(round(pxmmMap(linspace(0,3050,7)/1000))*1000));set(gca,'ytickLabel',ylab)
    xlim(xPad);ylim(yPad)
    set(gca,'Units','normalized','Position',[0.09 0.09 0.82 0.82]);   % even canvas margin on all 4 sides
    ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps);
end

crackEdge=edgeList;

if animate
    kk=cellfun(@(x) plot(x(:,2),x(:,1),':k','LineWidth',3),crackRaw_c);
    ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps);

    CG=plot([node(crackEdge(:,1),2)';node(crackEdge(:,2),2)'],[node(crackEdge(:,1),1)';node(crackEdge(:,2),1)'],'g--','LineWidth',3);
    CG=[CG;cellfun(@(x) plot(x([1,end],2),x([1,end],1),'pr','MarkerSize',20,'MarkerFaceColor','r'),crackRaw_c)'];
    ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps);

end

if ~exist('obj_den','var')
    objCrack=union(cellfun(@(x) polybuffer(fliplr([node(x(:,1),:);node(x(:,2),:)]),'line',s),mat2cell(crackEdge,ones(size(crackEdge,1),1))));
end

Y = [0 colBW colBW 0]; X = [0 0 rowBW rowBW];
extBound = polyshape(Y, X);final=extBound;      % External Boundary
final_ws = final;
final = subtract(final,objCrack);           	% Subracting object from frame
final_ws = subtract(final_ws,objCrack);      	% Subracting object from frame
final_ws=polyclean(final_ws);
final_work=final_ws;
critP=0;WSarea=572.635;
PathEdge=[];

if animate
    objplt=plot(objCrack,'FaceColor','y','FaceAlpha',0.1,'LineStyle','none');
    reg=plot(final_ws,'FaceAlpha',0.1,'LineStyle','none');
    ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps);
end


%% Main
tmain=tic;

if s~=1524
%%% MCD
working = final_ws;iMor=iMor+1;
[subcritP,splitReg_work,splitReg,splitEdge]=MCD(working,final_work,final_ws,critP,smflag);

subcritP=fliplr(subcritP);

%%% Reeb Graph
splitReg_work = polyclean(splitReg_work);

[reebEdge,~,~,~,~]=Reeb(splitReg_work,subcritP,splitEdge);

rMiss=find(~ismember((1:length(subcritP))',unique(reebEdge),'rows'));r_rmin=[];indj=[];%rAreaN=[];
if ~isempty(rMiss)
    for i=1:length(rMiss)
        rdis=spdist(subcritP(rMiss(i),:),subcritP);rdis(~rdis)=inf;[~,rmin(i)]=min(rdis);

        if any(rdis<s)
            if numel(find(ismember(reebEdge,reebEdge(sum(ismember(reebEdge,rmin(i)),2)>0,:),'rows')))<3
                rAreaN=find(ismember(reebEdge,reebEdge(sum(ismember(reebEdge,rmin(i)),2)>0,:),'rows'))';%rMiss(i);rAreaN;
                if ~isempty(rAreaN)
                    rNodeN=reebEdge(rAreaN,:);
                    if size(rNodeN,1)>1
                        reebEdge(rAreaN,:)=[];r_rmin=[r_rmin;rmin(i)];%(i,:)
                        rArea=union(splitReg_work(rAreaN));splitReg_work(rAreaN)=[];%(i,:)
                        splitReg_work=[splitReg_work;rArea];reebEdge=[reebEdge;rNodeN(~ismember(rNodeN,rmin(i)))'];
                    end
                end
            end
        else
            indj=[indj,i];
        end
    end

    rMiss(indj)=[];
    subcritP(unique([find(~ismember(1:length(subcritP),unique(reebEdge))');rMiss;r_rmin]),:)=[];
end

[reebEdge,reebCell,reeb,reebwall,remreg]=Reeb(splitReg_work,subcritP,splitEdge);
if ~isempty(remreg)
    splitReg_work(remreg)=[];
    reebCell = reebCell - arrayfun(@(c) sum(remreg < c), reebCell);   % remap reebCell indices past removed cells
end

% For route='rpp': route directly over MCD cells using the geometry-routed cell Rural-Postman
% and skip the reeb critical-ordering + combined ChinesePostman + Boustrophedon_CellCon below.
% For route='native': use the reeb-graph ChinesePostman route (original method).
if strcmp(route,'rpp')
    PathEdge = celled_rpp(node, edgeList, splitReg_work, s, a);
    plpath=[]; plpathA=[]; plmcdcrtt=[]; reebt=[]; plmcdcrt=[]; plreeb=[]; plnode=[]; plnodet=[];  % empty handles so cleanup delete() below is a no-op
end
if ~strcmp(route,'rpp')   % native-only block (reeb plotting + critical ordering + remEdge)

if animate
    if iMor~=1;if exist('pllplot','var'); delete(pllplot);end;pllplot=plot(intersect(extBound,polybuffer(PathEdge(2:end,:),'line',s)),'FaceColor','y','FaceAlpha',0.1,'LineStyle','none');end
    delete(reg);reg=plot(splitReg_work,'FaceAlpha',0.1,'LineStyle','-','LineWidth',2);for r=1:length(reg);reg(r).EdgeColor=reg(r).FaceColor;end
    plreeb=plot(reeb(1:2:length(reeb(:,1)),:)',reeb((1:2:length(reeb(:,1)))+1,:)','LineWidth',2);
    nuum1=int2str((1:length(reeb(:,1))/2)');for rr=1:length(nuum1);nuum2(rr,:)=regexprep(nuum1(rr,:),'.','_$0');end
    reebt=text((reeb(1:2:length(reeb(:,1)),ceil(end/2))+20)',(reeb((1:2:length(reeb(:,1)))+1,ceil(end/2))+70)', [repmat('E',length(reeb(:,1))/2,1), nuum2],'FontSize',20);clear nuum1 nuum2;

    nuum1=int2str([1:length(subcritP)]');for rr=1:length(nuum1);nuum2(rr,:)=regexprep(nuum1(rr,:),'.','_$0');end
    plmcdcrtt=text(subcritP(:,2)+20,subcritP(:,1)+70,[repmat('C',length(subcritP),1) nuum2],'FontSize',20);clear nuum1 nuum2;
    plcurPP=[];
    if exist('plmcdcrt','var');delete(plmcdcrt);end;plmcdcrt=plot(subcritP(:,2),subcritP(:,1),'r.','MarkerSize',50);
    xlim(xPad);ylim(yPad)
    snapnow
end


if animate; ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps); end

%%% Combining the all egdes and Generatign the Path, LP Problem

allNode = [node; subcritP];No_node = length(allNode);

reebEdge=max(crackEdge(:))+reebEdge;          % Keep in mind

reebNodes=unique(reebEdge);cell=splitReg_work;
remEdge=[];
if length(cell)>1
    for i = 1:length(cell)
        cc=refinePoly(cell(i),1);
        [~,dd]=dsearchn(cc.Vertices,fliplr(allNode(reebNodes,:)));
        if sum(dd<a)>0
            remEdge=[remEdge;combnk(reebNodes(dd<a),2)];
        end
    end
    if ~isempty(remEdge)
        remEdge=[remEdge;remEdge(:,2),remEdge(:,1)];
    end
end
end   % end native-only block (reeb plotting + critical ordering + remEdge)
else
    allNode = node;No_node = length(allNode);
    GCC=1;
end
if ~strcmp(route,'rpp')   % native-only routing (combine edges + ChinesePostman + Boustrophedon_CellCon + start/end cut)
if GCC
    allEdge = [crackEdge];
else
    allEdge = [crackEdge;reebEdge];
end
No_edge = length(allEdge);

allEdge=[allEdge;allEdge(:,2),allEdge(:,1)];
adj = full(sparse(allEdge(:, 1), allEdge(:, 2), 1, No_node, No_node));
Dist =squareform(pdist(allNode));
% Weight each Reeb edge by the equivalent boustrophedon sweep length for its cell
% (cell area / sweep-strip width) rather than raw critical-point distance.
if ~GCC
    for rb = 1:length(reebEdge)
        Dist(reebEdge(rb,1), reebEdge(rb,2)) = splitReg_work(rb).area/(2*(s/sqrt(2)));
        Dist(reebEdge(rb,2), reebEdge(rb,1)) = splitReg_work(rb).area/(2*(s/sqrt(2)));
    end
end
AdjMax=adj.*Dist;

remEdge=[];

[Path, weight, add,st]=ChinesePostman(adj,AdjMax,Dist,[],remEdge);
Pathe=[Path(1:end-1)' Path(2:end)'];adds=Pathe(ismember(Pathe, [add; fliplr(add)], 'rows'),:);

if animate
    [~,mm]=max(spdist2(allNode(adds(:,1),:),allNode(adds(:,2),:)));addd=adds(mm,:);
    adds(mm,:)=[];
    plpath=plot([allNode(adds(:,1),2) allNode(adds(:,2),2)]',[allNode(adds(:,1),1) allNode(adds(:,2),1)]','--b','LineWidth',3);plpathA=[];
    for dA=1:length(adds(:,1))
        plpathA=[plpathA;drawArrowHead(fliplr(allNode(adds(dA,1),:)),fliplr(allNode(adds(dA,2),:)),'b')];
    end
    plnode=plot(node(:,2),node(:,1),'pr','MarkerSize',12,'MarkerFaceColor','r');   % crack-node markers
    nuum1=int2str([1:length(node)]');for rr=1:length(nuum1);nuum2(rr,:)=regexprep(nuum1(rr,:),'.','_$0');end
    plnodet=text(node(:,2)+20,node(:,1)+70,[repmat('N',length(node),1) nuum2],'FontSize',20);clear nuum1 nuum2;   % crack-node labels
    xlim(xPad);ylim(yPad)
end

if animate; ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps); end

Path = [Path(1:end-1)',Path(2:end)'];

if GCC
    ppPath=Path;
    ppPath(~(ismember(Path,add,'rows')|ismember(fliplr(Path),add,'rows')),:)=1;
    [~,rrem]=max(spdist2(allNode(ppPath(:,1),:),allNode(ppPath(:,2),:)));
    Path(rrem,:)=[];

    PathEdge=fliplr(allNode([Path(:,1);Path(end,2)],:));subXY=PathEdge;
else
    %%SCC
    %%%%%
    wall_fol= zeros(size(Path,1),1);
    %%%%%
    subXY = Boustrophedon_CellCon(Path,splitReg_work,[],[],critP,wall_fol,true,sim, r1-a, size(node,1));

    %%% New Start End Optimization
    subXY(end,:)=subXY(1,:); % Making/Verify Euler
    [~,m]=max(spdist2(subXY(1:end-1,:),subXY(2:end,:)));subXY(end,:)=[];
    subXY=circshift(subXY,size(subXY,1)-m);
    %%%

    PathEdge = [PathEdge;subXY];
end
end   % end native-only routing guard


ttt=[ttt,toc(tmain)];disp(ttt);ttt=sum(ttt);

if animate
    delete(plpath);delete(plpathA);delete(plmcdcrtt);delete(reebt);delete(plmcdcrt);delete(plreeb);delete(reg);delete(objplt);delete(kk);delete(CG)
    if exist('plnode','var');delete(plnode);end;if exist('plnodet','var');delete(plnodet);end   % clear crack-node markers/labels
    reg=plot(splitReg_work,'FaceColor','w','FaceAlpha',0.1,'LineStyle','-','LineWidth',2,'EdgeColor','#D95319');
    plot(PathEdge(:,1),PathEdge(:,2),'k--','LineWidth',2);%,'Color','#D95319')
    plot(PathEdge(end,1),PathEdge(end,2),'v','Color','#77AC30','MarkerSize',15,'MarkerFaceColor','#77AC30')
    plot(PathEdge(1,1),PathEdge(1,2),'^','Color','#0072BD','MarkerSize',15,'MarkerFaceColor','#0072BD')
    % black arrowheads at each straight-edge midpoint, pointing end-to-end (matches the final plot)
    verts = reducepoly(PathEdge, 0.012);
    for vi = 1:size(verts,1)-1
        v0 = verts(vi,:); v1 = verts(vi+1,:);
        if norm(v1-v0) < 160; continue; end
        plpathA = [plpathA; drawArrowHead(v0, 0.5*(v0+v1), 'k')];
    end
    xlim(xPad);ylim(yPad)
    ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps);
end

if animate
   for ren=1:length(PathEdge(:,1))
        curPt=PathEdge(ren,:);
        if exist('pctp','var');delete(pctp);delete(pctps);end
        pctps=plot(polybuffer(curPt,'points',s),'FaceColor','y','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);%pctps2=plot(polybuffer(curPt,'points',s*sqrt(2)),'FaceColor','r','FaceAlpha',0.1,'LineStyle','--','LineWidth',1);
        pctp=plot(polybuffer(curPt,'points',a*sqrt(2)),'FaceColor','r','FaceAlpha',0.3,'LineStyle','-','LineWidth',1);
        ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps);
   end
end


%% Calculations

pathLength=total_length(PathEdge);scclen=pxinMap(pathLength)/12;
areaCover = 2*pxinMap(s/sqrt(2))/12*scclen;
Ovp_Area =areaCover-WSarea;
coverPercent=areaCover/WSarea;
ovlapPercent=Ovp_Area/WSarea;

areaCover1=(area(polybuffer(PathEdge,'lines',s/sqrt(2)))/(2898*3050))*100;

res = [areaCover1,den(k),ttt,scclen,areaCover];disp(res)

ress=[ress;res];
end

finalPlot(PathEdge, crackRaw_c, crackGen, vizName, 'SCC');
close all


%% Functions

function polyout=polyclean(polyin)
    global s
    polyout=polyshape();
    if isscalar(polyin)
        poly = regions(polyin);
        poly=poly(fix(poly.area*1e-03)>10);%poly.area>100
        polyout = regJoin(poly);
    else
        poly = polyin;
        poly=poly(fix(poly.area*1e-03)>10);%poly.area>100
        polyout = poly;
    end
end

function polyout = regJoin(polyin)
    polyout=polyshape();
    for i=1:length(polyin)
        polyout=addboundary(polyout,polyin(i).Vertices);
    end
end

function ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps)
    ppn = ppn + 1;
    if ~writeFrames && ~makeGif, drawnow; return; end   % live animate update when nothing is written
    fr = getframe(gcf);
    if writeFrames
        imwrite(fr.cdata, ['Results/GIF/' vizName '_' num2str(ppn) '.png']);
    end
    if makeGif
        [A, cmap] = rgb2ind(fr.cdata, 256);
        if ppn == 1
            imwrite(A, cmap, gifPath, 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
        else
            imwrite(A, cmap, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
        end
    end
end

function finalPlot(PathEdge, crackRaw_c, BW, vizName, suffix)
% finalPlot  Styled "Final Path" figure: metre axes (Y inverted), dashed coverage
%   path + arrowheads, orange cracks, blue/green start/end. Saved to
%   Results/<suffix>/<vizName>_<suffix>.png (always, regardless of the animate flag).
    [rowBW, colBW] = size(BW);
    fig = figure('Color','w','Visible','off');
    ax = axes('Parent', fig, 'Color', 'w');
    hold(ax, 'on');
    for ci = 1:length(crackRaw_c)
        cr = crackRaw_c{ci};
        plot(ax, cr(:,2), cr(:,1), '-', 'Color', '#D95319', 'LineWidth', 1.5);
    end
    plot(ax, PathEdge(:,1), PathEdge(:,2), '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2);
    verts = reducepoly(PathEdge, 0.012);   % straight edges; one arrowhead per edge midpoint
    for vi = 1:size(verts,1)-1
        v0 = verts(vi,:); v1 = verts(vi+1,:);
        if norm(v1-v0) < 160; continue; end
        drawArrowHead(v0, 0.5*(v0+v1), [0.1 0.1 0.1]);
    end
    n = size(PathEdge,1);
    drawArrowHead(2*PathEdge(1,:) - PathEdge(min(5,n),:), PathEdge(1,:), [0 0.4470 0.7410]);
    drawArrowHead(PathEdge(max(1,n-4),:), PathEdge(end,:), [0.4660 0.6740 0.1880]);
    off = 0.025 * colBW;
    text(ax, PathEdge(1,1)+off, PathEdge(1,2)-off, 'Start Point', 'FontSize', 10, 'Color', 'k');
    text(ax, PathEdge(end,1)+off, PathEdge(end,2)+off, 'End Point', 'FontSize', 10, 'Color', 'k');
    set(ax, 'XTick', 0:500:colBW, 'XTickLabel', string(0:numel(0:500:colBW)-1));
    set(ax, 'YTick', 0:500:rowBW, 'YTickLabel', string(0:numel(0:500:rowBW)-1));
    xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)');
    set(ax, 'YDir', 'reverse', 'DataAspectRatio', [1 1 1]);
    % tight data limits (no inner padding) — the image fills the axis box
    xlim(ax, [0, colBW]); ylim(ax, [0, rowBW]);
    set(ax, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'LineWidth', 0.75);
    title(ax, 'Final Path', 'FontWeight', 'bold', 'Color', 'k');
    % even canvas margin around the axis (room for labels left/bottom, title top, matching space right)
    set(fig, 'Position', [100 100 900 900]);
    set(ax, 'Units', 'normalized', 'Position', [0.09 0.09 0.82 0.82]);
    savedir = ['Results/' suffix];
    if ~exist(savedir, 'dir'); mkdir(savedir); end
    saveas(fig, [savedir '/' vizName '_' suffix '.png']);
    close(fig);
end
