%% OnlineSCC — Sensor-based Complete Coverage with online crack detection
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

%% Initial Setup

clc
close all
clearvars -except obj_den k

ress=[];
for sa=1
clearvars -except sa ress

global reebEdge reebCell crackEdge allNode vertical spdist spdist2 s a sMask sRan total_length

inpxMap = @(x) fix(((x)*25.4)/2);
pxinMap = @(x) round((x)*2/25.4,1);
mmpxMap = @(x) fix(((x)/2));
pxmmMap = @(x) round((x)*2,1);

sim =false;

%% Visualization
animate     = false;   % show the live animation figure (plot even without writing frames/GIF)
makeGif     = false;   % stitch frames into Results/GIF/<map>.gif  (also triggers the figure)
writeFrames = false;   % save per-step PNG frames to Results/GIF/  (also triggers the figure)
fps         = 14;      % GIF frame rate
animate = animate || writeFrames || makeGif;   % figure is created when animate is set OR frames/GIF are wanted

warning off

% Robot dimensions from robot_config.json (edit there to reconfigure the robot)
cfg   = config_loader('robot1');
botD  = cfg.base_diameter_in; footD = cfg.footprint_diameter_in; sensD = cfg.sensor_diameter_in;
r1=cfg.r1;                               % Robot Radius
a=cfg.a;                                 % Footprint Radius
s=cfg.s;                                 % Sensor Range Radius
dir_map=[-1 -1;-1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0 ; 1 1];

sMask = @(BW,rr,cp) deal(BW(bound(round(cp(2)-(rr)),1,size(BW,1)):bound(round(cp(2)+(rr)),1,size(BW,1)),...
    bound(round(cp(1)-(rr)),1,size(BW,2)):bound(round(cp(1)+(rr)),1,size(BW,2))),[bound(round(cp(2)-(rr)-1),0,size(BW,1)),bound(round(cp(1)-(rr)-1),0,size(BW,2))]);

sMaskid = @(BW,rr,cp) deal([bound(round(cp(2)-(rr)),1,size(BW,1)):bound(round(cp(2)+(rr)),1,size(BW,1))],...
    [bound(round(cp(1)-(rr)),1,size(BW,2)):bound(round(cp(1)+(rr)),1,size(BW,2))]);

sRan = @(t,rr,cp) bound(round(cp-(rr)),1,t):bound(round(cp+(rr)),1,t);

crackGen = 0;numItr=0; ppn=0;

%% Importing Map

den  = [35,45,50,65,80,90,95,100];   % density values
sig  = [5,10,20];                    % Gaussian sigmas
Gau  = 0;                            % 0 = Uniform, 1 = Gaussian
dd   = 8;                            % density index into den
mapN = 1;                            % map number (Uniform variant)
if ~Gau
    img_n = ['myCrack' num2str(dd) '_' num2str(den(dd)) '_' num2str(mapN)];        % Uniform map
    BW = imread(['CrackMaps/Uniform/' img_n '.png']);                              % load the PNG
else
    Gaussb = 6;                      % Gaussian folder 1..8
    b      = 3;                      % sigma index into sig (1..3)
    img_n  = ['Gaussian' num2str(Gaussb) '/myCrackGauss_s' num2str(sig(b)) '_' num2str(den(dd))];  % Gaussian map (folder-qualified)
    load(['CrackMaps/' img_n '.mat']);                                             % loads crackGen (the crack mask)
    crackGen(end,end)=0;                                                           % clear sentinel
end

vizName = strrep(img_n,'/','_');
gifPath = ['Results/GIF/' vizName '.gif'];
if ~exist('Results/GIF','dir'); mkdir('Results/GIF'); end
delete(['Results/GIF/' vizName '_*.png']);
if exist(gifPath,'file'); delete(gifPath); end

%%

sm_flag=1; % Polygon Smoothing flag

if crackGen==0                     % ischar(img_n)
    BW = imbinarize(BW);           % Converting the image into a binary image
    BW = BW(:,:,1);

    BW = bwareaopen(BW, 50);       % Removes all objects that have fewer than 50 pixels from the binary image BW
    BW = padarray(~BW,[1 1]);      % Adds a boundary around the map

    IM2 = imcomplement(BW);        % Computes the complement of the image BW

else
    BW = crackGen;
end

%%%

BW2=bwmorph(BW,'fill');             % New Fills 1 pixel islands
BW2 = bwskel(BW2>0);                % Skeliton function removes the diagonal pixel
BW3 = bwmorph(BW2, 'spur', 10);     % Removes spur pixels
BW_working = BW3;

crackGen = zeros(size(BW_working));
crackGGen = zeros(size(BW_working));
[mBW,nBW,vBW]=find(sparse(BW_working));

if animate
    fig1 = figure('Color','w'); imshow(~crackGen); hold on
    set(gca,'Color','w','XColor','k','YColor','k');
    set(fig1,'Position',[100 100 900 900]);
    if ~isfolder(fullfile('Results', 'GIF'))
        mkdir(fullfile('Results', 'GIF'));
    end
end

%%%

[rowBW, colBW]= size(BW3);
% Tight axis limits (image extent, no inner padding); even space is added OUTSIDE the axis box
xPad  = [0, colBW];
yPad  = [0, rowBW];
tol =0.00001;
colors={'y','m','c','r','g','b'};

if animate
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
end

vertical = @(P) all(P(1,1)==P(:,1));
spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);
spdist2 = @(Ps1,Ps2) sqrt((Ps1(:,1)-Ps2(:,1)).^2 + (Ps1(:,2)-Ps2(:,2)).^2);
total_length = @(Ps) sum(sqrt(sum(diff(Ps).*diff(Ps),2)));

Y = [0 colBW colBW 0]; X = [0 0 rowBW rowBW];
extBound = polyshape(Y, X);

final=extBound;                                 % External Boundary
final_work=final;
critP=0;loop1=true;ss=polybuffer([0,0],'points',s);aa=polybuffer([0,0],'points',a);
PathEdge=[0,0];crackEdge=[];ttt=[];iMor=0;

if animate;reg=plot(final,'FaceAlpha',0.1,'LineStyle','-','LineWidth',2);reg.EdgeColor=reg.FaceColor;plot(s/sqrt(2),s/sqrt(2),'v','Color','#0072BD','MarkerSize',15,'MarkerFaceColor','#0072BD');end

%% Main
tStart = tic;

while any(final.regions.area>2*aa.area)
tic
if exist('plreeb','var')
    delete(plreeb);delete(reebt);delete(plmcdcrt);delete(plmcdcrtt);
end

% working holds the eroded/dilated free region; final_work holds the current free polygon
%%% MCD
working = polybuffer(final_work,-s/5);                                                % Impoting Final Work to working.
working = polybuffer(working,s/5);
bblob=regions(working);working=bblob(bblob.area>area(polybuffer([0,0],'points',a)));  % Removing regions < footpring area
if size(regions(final_work),1)>1                                                      % Split, remove regions < 1/5 of sensor range and rejoin them
    final_work=regions(final_work);
    final_work=final_work(~(area(polybuffer(final_work,-s/5))==0));
    final_work=regJoin(final_work);
end

if size(working,1)~=size(final_work,1); working=final_work;end                        % If Working ~ Final, resign the variable.

[subcritP,splitReg_work,splitReg,splitEdge]=MCD(working,final_work,final,critP,sm_flag);
splitReg_work = polyclean(splitReg_work);
%%%

subcritP=fliplr(subcritP);iMor=iMor+1;

if isempty(subcritP); break;end

%--------Reeb
[reebEdge,~,~,~,remreg]=Reeb(splitReg_work,subcritP,splitEdge);%subcritP=subcritP(unique(reebEdge),:);
if ~isempty(remreg);splitReg_work(remreg)=[];end

rMiss=find(~ismember((1:length(subcritP))',unique(reebEdge(:)),'rows'));r_rmin=[];indj=[];%rAreaN=[];

adj = full(sparse(reebEdge(:, 1), reebEdge(:, 2), 1, size(subcritP,1), size(subcritP,1)));
adj=adj+adj';
rCont= 1:length(subcritP);
rCont=rCont(sum(adj)==2);clear adj;

rCont(ismember(rCont,rMiss))=[];

if ~(isempty(rMiss)&&isempty(rCont))
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
                        splitReg_work=[splitReg_work;polybuffer(polybuffer(rArea,1),-1)];reebEdge=[reebEdge;rNodeN(~ismember(rNodeN,rmin(i)))'];
                    end
                end
            end
        else
            indj=[indj,i];
        end
    end

    for i=1:length(rCont)
        rAreaN=find(ismember(reebEdge,reebEdge(sum(ismember(reebEdge,rCont(i)),2)>0,:),'rows'))';
        if ~isempty(rAreaN)
            rNodeN=reebEdge(rAreaN,:);
            if size(rNodeN,1)>1
                reebEdge(rAreaN,:)=[];r_rmin=[r_rmin;rCont(i)];%(i,:)
                rArea=union(splitReg_work(rAreaN));splitReg_work(rAreaN)=[];%(i,:)
                splitReg_work=[splitReg_work;polybuffer(polybuffer(rArea,1),-1)];reebEdge=[reebEdge;rNodeN(~ismember(rNodeN,rCont(i)))'];
            end
        end
    end
    subcritP(unique([rMiss(:);rCont(:);r_rmin(:)]),:)=[];clear rmin
end

[reebEdge,reebCell,reeb,reebwall,remreg]=Reeb(splitReg_work,subcritP,splitEdge);
if ~isempty(remreg);splitReg_work(remreg)=[];end

%------------------------------------------

allNode = subcritP;
No_node = length(allNode);
allEdge = reebEdge;
G=graph(allEdge(:,1),allEdge(:,2));adj=full(adjacency(G));


if length(conncomp(graph(reebEdge(:,1),reebEdge(:,2)), 'OutputForm', 'cell'))>1
    t1=sum(reshape(subcritP(reebEdge(:),2)<PathEdge(end,1),[],2),2)>1;  % Cell with both the critical points to the left of the current point
    if 0%any(t1)%& t2
        rT = reebEdge(t1,:);[~,ind]=min(spdist(PathEdge(end,:),fliplr(subcritP(rT,:))));
        ind=rT(ind);
    else
        [~,ind]=min(cellfun(@(x) min(spdist(PathEdge(end,:),x)),{splitReg_work(:).Vertices}));rc=reebEdge(reebCell==ind,:)';
        [~,k]=min(abs(PathEdge(end,1)-subcritP(rc,2)));ind=rc(k);
    end
else
    [~,ind]=min(cellfun(@(x) min(spdist(PathEdge(end,:),x)),{splitReg_work(:).Vertices}));rc=reebEdge(reebCell==ind,:)';
    [~,k]=min(abs(PathEdge(end,1)-subcritP(rc,2)));ind=rc(k);
end

[Path,wall_fol,~]=reeb_traversal(adj,subcritP,reebEdge,splitReg_work,ind);se=[];seP=[];

if animate
    if iMor~=1;if exist('pllplot','var'); delete(pllplot);end;pllplot=plot(intersect(extBound,polybuffer(PathEdge(2:end,:),'line',s)),'FaceColor','y','FaceAlpha',0.1,'LineStyle','none');end
    delete(reg);reg=plot(splitReg_work,'FaceAlpha',0.1,'LineStyle','-','LineWidth',2);for r=1:length(reg);reg(r).EdgeColor=reg(r).FaceColor;end
    plreeb=plot(reeb(1:2:length(reeb(:,1)),:)',reeb((1:2:length(reeb(:,1)))+1,:)','LineWidth',2);
    if any(wall_fol)
        for w_reeb=find(wall_fol)'
            plreeb=[plreeb;plot(reebwall(w_reeb*2-1,:)',reebwall(w_reeb*2,:)','--','Color',plreeb(w_reeb).Color,'LineWidth',2)];
        end
    end
    nuum1=int2str((1:length(reeb(:,1))/2)');for rr=1:length(nuum1);nuum2(rr,:)=regexprep(nuum1(rr,:),'.','_$0');end
    reebt=text((reeb(1:2:length(reeb(:,1)),ceil(end/2))+20)',(reeb((1:2:length(reeb(:,1)))+1,ceil(end/2))+70)', [repmat('E',length(reeb(:,1))/2,1), nuum2],'FontSize',20);clear nuum1 nuum2;

    nuum1=int2str([1:length(subcritP)]');for rr=1:length(nuum1);nuum2(rr,:)=regexprep(nuum1(rr,:),'.','_$0');end
    plmcdcrtt=text(subcritP(:,2)+20,subcritP(:,1)+70,[repmat('C',length(subcritP),1) nuum2],'FontSize',20);clear nuum1 nuum2;
    plcurPP=[];
    if iMor~=1;plcurPP=[plcurPP;plot([curPt(:,1),subcritP(Path(1,1),2)],[curPt(:,2),subcritP(Path(1,1),1)],'b-.','LineWidth',2.5)];end
    if any(~Path(:,3))
        for hh=find((~Path(:,3)))'
            plcurPP=[plcurPP;plot(subcritP(Path(hh,1:2),2),subcritP(Path(hh,1:2),1),'b-.','LineWidth',2.5)];
        end
    end
    if exist('plmcdcrt','var');delete(plmcdcrt);end;plmcdcrt=plot(subcritP(:,2),subcritP(:,1),'r.','MarkerSize',50);
    xlim(xPad);ylim(yPad)
    snapnow
end


if animate; ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps); end

%%% Boustrophedon Path
if length(splitReg_work)==1
    %%% Boustrophedon Path Generation
    subXY = Boustrophedon(Path,splitReg_work,se,seP,critP,wall_fol,false,sim);
else
    %%% Boustrophedon Path Generation with Cell Connection Algorithm
    subXY = Boustrophedon_CellCon(Path,splitReg_work,se,seP,critP,wall_fol,false,sim);
end

if animate
    pp=plot(subXY(:,1),subXY(:,2));
    delete(pp)
end


clear wall_fol;
if exist('plreeb','var')
   delete(plmcdcrt);delete(plmcdcrtt);delete(plreeb);delete(reebt);if exist('plcurPP','var');delete(plcurPP);end
end

if loop1
    zigzag=total_length(subXY);ziglen=pxinMap(zigzag)/12;
    WSarea = 2*pxinMap(s/sqrt(2))/12*ziglen;
end

subXY=[subXY,zeros(size(subXY,1),1)];

fill=false;bww=sparse([],[],[],2*s,2*s);

i=1;crack=0;clear preCrack;
test=0;fl=false;fflg=false;intPoints=zeros(size(BW3));d=0;d_lock=0;ww=[];
ttt=[ttt,toc];

if animate; hold on; ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps); end

while i<=length(subXY)
    disp(length(subXY));
    tic
    curPt=subXY(i,[1,2]);disp(curPt)
    if i == 15
        disp(curPt)
    end
    if subXY(i,3)
        [bwr,bwc]=sMaskid(BW_working,a,curPt);
        bww(bwr,bwc)= 0;
    end

    %%%New Extraction method
    [in,~]=inpolygon(mBW,nBW,polybuffer(fliplr(curPt),'points',s).Vertices(:,1),polybuffer(fliplr(curPt),'points',s).Vertices(:,2));
    for bwh=find(in)';crackGen(mBW(bwh),nBW(bwh))=1;crackGGen(mBW(bwh),nBW(bwh))=1;end
    m=mBW(in);n=nBW(in);v=ones(length(find(in)),1);

    %%%%
    I=crackGen;I_te=I;
    while true
        intPoints = intPoints + neighbor_count_points(I,2); %2: Branchpoints
        I_te(intPoints>0)=0;BW3(intPoints>0)=0;
        pixelP=neighbor_count_points(I_te,0); %0: Single Pixels
        I_te(pixelP>0)=0;BW3(pixelP>0)=0;

        if ~isequal(I,I_te)
           I=I_te;
        else
           break;
        end
    end
    crackGen=I;
    %%%%

    %%%%
    [eP,rP,~] = endP_ident(crackGen,BW3);endlogi=0;
    ttt=[ttt,toc];
    if size(rP,1)>1
        [~,line,~,~] = compCrack(crackGen,eP,dir_map,[]);
        endlogi = reshape(sum(cell2mat(cellfun(@(s) spdist(s,[line(:,1) line(:,2);line(:,3) line(:,4)]),num2cell( rP, 2),'un',0)')<5,2)>0,[],2);
    end
    %%%%

    tic
    if ~isempty(m) && (~all(m==m(1)) && ~all(n==n(1))) || (max(sum(endlogi,2))==2 && any(spdist2(line(:,1:2),line(:,3:4))>a))

    if animate
        if exist('pcp','var');delete(pcp);end
        [p1,p2]=find(sparse(crackGGen));
        pcp=plot(p2,p1,'k.','LineWidth',2);
        if exist('pctp','var');delete(pctp);delete(pctps);end%delete(pctps2);
        pctps=plot(polybuffer(curPt,'points',s),'FaceColor','y','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);%pctps2=plot(polybuffer(curPt,'points',s*sqrt(2)),'FaceColor','r','FaceAlpha',0.1,'LineStyle','--','LineWidth',1);
        pctp=plot(polybuffer(curPt,'points',a*sqrt(2)),'FaceColor','r','FaceAlpha',0.3,'LineStyle','-','LineWidth',1);
    end

    [eP,rP,cP] = endP_ident(crackGen,BW3);

        %%% New Test
        if size(rP,1)>1
            ttt=[ttt,toc];
            [crackRaw,~,~,~] = compCrack(crackGen,eP,dir_map,[]);
            tic
            if ~isempty(ww)
                ww1=ww(1:find(ismember(ww,subXY(i,1:2),'rows')),:);int=[];[int(:,1),int(:,2)]=find(intPoints>0);
                if ~isempty(ww1)
                    pp=polybuffer(subXY(1:i,1:2),'line',a+25);
                    for c=1:size(crackRaw,2)
                        if ~isempty(rP)
                            if sum(ismember(crackRaw{c}([1,end],:),rP,'rows'))==2 || length(find(spdist(crackRaw{c}(1,:),rP)<a | spdist(crackRaw{c}(end,:),rP)<a))==2%(all(spdist(crackRaw{c}(1,:),rP)<a) && all(spdist(crackRaw{c}(end,:),rP)<a))
                                in=isinterior(pp,fliplr(crackRaw{c}([1,end],:)));endPoints=crackRaw{c}([1,end],:);
                                c_chq = isinterior(pp,fliplr(crackRaw{c}));c_chq=length(find(c_chq))/length(c_chq)*100;  % Checks for correct crack to remove

                                if sum(in)==2 && c_chq>90
                                    post=[];
                                    for ee = 1:2
                                        row=endPoints(ee,1);
                                        col=endPoints(ee,2);
                                        I = crackGen;
                                        b = [I(row-1,col-1:col+1) ...
                                                I(row, col-1) I(row, col+1)...
                                                I(row+1,col-1:col+1)];
                                        post = [post;sum(b)];
                                    end
                                    [m,n,v]=find(sparse(crackGen)); v(ismember([m,n],crackRaw{c},'rows'))=0;
                                    rind=ismember([mBW,nBW],crackRaw{c},'rows');mBW(rind)=[];nBW(rind)=[];vBW(rind)=[];
                                    crackGen=full(sparse([m;round(rowBW)],[n;round(colBW)],[v;1]));crackGen(round(rowBW),round(colBW))=0;fl=true;
                                    [m,n,v]=find(sparse(BW_working)); v(ismember([m,n],crackRaw{c},'rows'))=0;
                                    BW_working=full(sparse([m;round(rowBW)],[n;round(colBW)],[v;1]));BW_working(round(rowBW),round(colBW))=0;
                                     [m,n,v]=find(sparse(BW3)); v(ismember([m,n],crackRaw{c},'rows'))=0;
                                     BW3=full(sparse([m;round(rowBW)],[n;round(colBW)],[v;1]));BW3(round(rowBW),round(colBW))=0;
                                     if ~isempty(int);BW3(int(isinterior(pp,fliplr(int)),1),int(isinterior(pp,fliplr(int)),2))=0;end;BW3=bwmorph(BW3,'clean'); % Cleans int points
                                    for ii=1:2; if post(ii)>1;crackGen(endPoints(ii,1),endPoints(ii,2))=1;BW_working(endPoints(ii,1),endPoints(ii,2))=1;BW3(endPoints(ii,1),endPoints(ii,2))=1;end;end
                                    stt=fliplr(curPt);crackGen=bwmorph(crackGen,'thin', inf);fflg=true;d_lock=[];
                                    if ~isempty(bww);[m,n,v]=find(sparse(bww)); v(ismember([m,n],crackRaw{c},'rows'))=0;
                                    bww=full(sparse(m,n,v));end
                                    [eP,rP,cP] = endP_ident(crackGen,BW3);
                                end
                            end
                        end
                    end
                end
            end
        end

        %%%

        if ~subXY(i,3)
            stt=fliplr(curPt);
        else
            d_lock=d_lock+1; %Records when stt!=cp
        end

        if ~exist('preCrack','var')
            ttt=[ttt,toc];
            [WP,flag,preCrack,tt] = ImagePlanning_oSCC(crackGen,a,s,stt,fliplr(curPt),eP,rP,cP,subXY(1:i,1:2),[]);
            ttt=[ttt,tt];tic
        else
            ttt=[ttt,toc];
            [WP,flag,preCrack,tt] = ImagePlanning_oSCC(crackGen,a,s,stt,fliplr(curPt),eP,rP,cP,subXY(1:i,1:2),preCrack);
            ttt=[ttt,tt];tic
        end
        numItr=numItr+1;

        if flag
            if exist('plpath','var');delete(plpath);delete(plpathA);end
            if i==1
                d=1;
            else
                if ~subXY(i-1,3) || fl
                    d=1;fl=0;
                else
                    d=d+1;
                end
            end
            subXY(i+1:end,:)=[];crack = 1;
            if d==1
                subXY=[subXY;[fliplr(WP(d:end,[1,2])),WP(d:end,3)]];
            elseif ismember(subXY(end,1:2),fliplr(WP(:,1:2)),'rows')
                indWP=find(ismember(fliplr(WP(:,1:2)),subXY(end,1:2),'rows'))+1;
                subXY=[subXY;[fliplr(WP(indWP:end,[1,2])),WP(indWP:end,3)]];
            else
                % Find the nearest segment on the crack-fill waypoint polyline and
                % resume from its forward end vertex; if the robot has not yet passed
                % the segment's midpoint, resume from its start vertex instead.
                WPf=fliplr(WP(:,1:2));
                if size(WPf,1)<2
                    mx=1;
                else
                    A=WPf(1:end-1,:);B=WPf(2:end,:);AB=B-A;
                    denom=sum(AB.^2,2);denom(denom==0)=1;
                    last=subXY(end,1:2);
                    tt=sum((last-A).*AB,2)./denom;
                    ttc=min(max(tt,0),1);
                    proj=A+ttc.*AB;
                    dseg=sum((proj-last).^2,2);
                    [~,k]=min(dseg);
                    if tt(k)>0.5;mx=k+1;else;mx=k;end
                end
                subXY=[subXY;[fliplr(WP(mx:end,[1,2])),WP(mx:end,3)]];
            end

            ww=fliplr(WP(:,[1,2]));
            if animate;plpath=plot(subXY(:,1),subXY(:,2),'--','Color','#D95319','LineWidth',2.5);plpath=[plpath;plot(PathEdge(2:end,1),PathEdge(2:end,2),'--','Color','#D95319','LineWidth',2.5)];plpathA=drawArrowHead(subXY(end-2,1:2),subXY(end,1:2),'#D95319');xlim(xPad);ylim(yPad);end
            snapnow; if animate; ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps); end
        else
            if d_lock==1
                d=1;
            else
                d=d+1;
            end
            if animate;if exist('plpath','var');delete(plpath);delete(plpathA);end; plpath=plot(subXY(:,1),subXY(:,2),'--','Color','#D95319','LineWidth',2.5);plpath=[plpath;plot(PathEdge(2:end,1),PathEdge(2:end,2),'--','Color','#D95319','LineWidth',2.5)];plpathA=drawArrowHead(subXY(end-2,1:2),subXY(end,1:2),'#D95319');xlim(xPad);ylim(yPad);end
            if animate; ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps); end
        end
    else
        [eP,rP,cP] = endP_ident(crackGen,BW3);
        %%% New Test
        if size(rP,1)>1
            ttt=[ttt,toc];
            [crackRaw,~,~,~] = compCrack(crackGen,eP,dir_map,[]);
            tic
            if ~isempty(ww)
                ww1=ww(1:find(ismember(ww,subXY(i,1:2),'rows')),:);
                if ~isempty(ww1)
                    pp=polybuffer(subXY(1:i,1:2),'line',a+25);
                    for c=1:size(crackRaw,2)
                        if ~isempty(rP)
                            if sum(ismember(crackRaw{c}([1,end],:),rP,'rows'))==2
                                in=isinterior(pp,fliplr(crackRaw{c}([1,end],:)));endPoints=crackRaw{c}([1,end],:);
                                if sum(in)==2
                                    post=[];
                                    for ee = 1:2
                                        row=endPoints(ee,1);
                                        col=endPoints(ee,2);
                                        I = crackGen;
                                        b = [I(row-1,col-1:col+1) ...
                                                I(row, col-1) I(row, col+1)...
                                                I(row+1,col-1:col+1)];
                                        post = [post;sum(b)];
                                    end
                                    [m,n,v]=find(sparse(crackGen)); v(ismember([m,n],crackRaw{c},'rows'))=0;
                                    rind=ismember([mBW,nBW],crackRaw{c},'rows');mBW(rind)=[];nBW(rind)=[];vBW(rind)=[];
                                    crackGen=full(sparse([m;round(rowBW)],[n;round(colBW)],[v;1]));crackGen(round(rowBW),round(colBW))=0;fl=true;
                                    [m,n,v]=find(sparse(BW_working)); v(ismember([m,n],crackRaw{c},'rows'))=0;
                                    BW_working=full(sparse([m;round(rowBW)],[n;round(colBW)],[v;1]));BW_working(round(rowBW),round(colBW))=0;
                                    for ii=1:2; if post(ii)>1;crackGen(endPoints(ii,1),endPoints(ii,2))=1;BW_working(endPoints(ii,1),endPoints(ii,2))=1;end;end
                                    stt=fliplr(curPt);crackGen=bwmorph(crackGen,'thin', inf);fflg=true;
                                    if ~isempty(bww);[m,n,v]=find(sparse(bww)); v(ismember([m,n],crackRaw{c},'rows'))=0;
                                    bww=full(sparse(m,n,v));end
                                    [eP,rP,cP] = endP_ident(crackGen,BW3,org);
                                end
                            end
                        end
                    end
                end
            end
        end

        if animate
            if exist('pcp','var');delete(pcp);end
            [p1,p2]=find(sparse(crackGGen));
            pcp=plot(p2,p1,'k.','LineWidth',2);
            if exist('pctp','var');delete(pctps);delete(pctp);end
            pctps=plot(polybuffer(curPt,'points',s),'FaceColor','y','FaceAlpha',0.1,'LineStyle','-','LineWidth',1);
            pctp=plot(polybuffer(curPt,'points',a*sqrt(2)),'FaceColor','r','FaceAlpha',0.3,'LineStyle','-','LineWidth',1);
        end

        if animate;if exist('plpath','var');delete(plpath);delete(plpathA);end; plpath=plot(subXY(:,1),subXY(:,2),'--','Color','#D95319','LineWidth',2.5);plpath=[plpath;plot(PathEdge(2:end,1),PathEdge(2:end,2),'--','Color','#D95319','LineWidth',2.5)];plpathA=drawArrowHead(subXY(end-2,1:2),subXY(end,1:2),'#D95319');xlim(xPad);ylim(yPad);end
        if animate; ppn = snapFrame(ppn, gifPath, vizName, writeFrames, makeGif, fps); end

        %%%

    end
    i=i+1;

    if i>=length(subXY(:,1)) && fflg
        subXY=[subXY;subXY(end,:)];fflg=false;
    end
end


obj = polybuffer(fliplr(subXY(:,[2,1])),'line',s);
[px,py]=obj.boundary;
if ~(px(1)==px(end)&&py(1)==py(end));px(end+1)=px(1);py(end+1)=py(1);end
obj.Vertices=DecimatePoly([px,py],[1 1],false);

final = subtract(final,obj);final_work = subtract(final_work,obj);

splitReg_work=subtract(splitReg_work,obj);

PathEdge=[PathEdge;subXY(:,[1,2])];
final = polyclean(final);

if length(splitReg_work)==1
    final_work=splitReg_work;
    final_work = polyclean(final_work);
else
    splitReg_work = polyclean(splitReg_work);
    if ~isempty(splitReg_work)
        final_work = regCombine(splitReg_work);
    end
end

final_work = polyclean(final_work);

critP=PathEdge(end,:);

if animate
    delete(reg)
    reg=plot(final,'FaceAlpha',0.1,'LineStyle','none');loop1=false;
    snapnow
end

ttt=[ttt,toc];
end

tMul = sum(ttt);tEnd = toc(tStart);disp(tMul); disp(tEnd)
%% Calculations

pathLength=total_length(rmmissing(PathEdge(2:end,:)));oscclen=pxinMap(pathLength)/12;
areaCover = 2*pxinMap(s/sqrt(2))/12*oscclen;
Ovp_Area =areaCover-WSarea;
coverPercent=areaCover/WSarea;
ovlapPercent=Ovp_Area/WSarea;

res = [numItr,den(dd),tMul,oscclen,areaCover];res=round(res,3);disp(res)

finalPlot(PathEdge(2:end,:), crackGGen, BW3, vizName, 'OnlineSCC');
close all

ress=[ress;res];
end




function polyout=polyclean(polyin)
    polyout=polyshape();
    if isscalar(polyin)
        poly = regions(polyin);
        poly=poly(fix(poly.area*1e-03)>5);%poly.area>100 10
        polyout = regJoin(poly);
    else
        poly = polyin;
        poly=poly(fix(poly.area*1e-03)>5);%poly.area>100 10
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
    polyout = polybuffer(union(polybuffer(polyin,20)),-20);
end

function y = bound(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
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

function finalPlot(PathEdge, crackGen, BW, vizName, suffix)
% finalPlot  Styled "Final Path" figure: metre axes (Y inverted), dashed coverage
%   path + arrowheads, orange cracks, blue/green start/end. Saved to
%   Results/<suffix>/<vizName>_<suffix>.png (always, regardless of the animate flag).
    [rowBW, colBW] = size(BW);
    fig = figure('Color','w','Visible','off');
    ax = axes('Parent', fig, 'Color', 'w');
    hold(ax, 'on');
    [cr, cc] = find(crackGen);
    if ~isempty(cr)
        plot(ax, cc, cr, '.', 'Color', '#D95319', 'MarkerSize', 1);
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
