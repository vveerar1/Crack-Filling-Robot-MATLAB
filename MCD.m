%% Morse Cell Decomposition
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [critPT,polyout_work,polyout,splitEdge]=MCD(polyin_buffed,polyin_work,polyin,nodeend,flag)
%
% Computes a Morse Decomposition of the workspace. 
%
% INPUTS:
%   polyin_buffed = Polygon in the shape of the workspace with a minkowski sum of 'a'.  
%   polyin_work = Polygon in the shape of the workspace, working. 
%   polyin = Polygon in the shape of the workspace.
%   nodeend = Current Position of the robot. 
%   flag = Flag to smooth the edges of the imput polygon. 
%
% OUTPUTS:
%   critPT = Critcal points in the workspace. 
%   polyout_work = Decomposed cells in the workspace. 
%   polyout = Decomposed cells in the workspace.
%   splitEdge = Splitline at the critial points. 

    if length(polyin_buffed)==1; polyin_buffed=regions(polyin_buffed);end
    if length(polyin_buffed)>1;polyin_work=regions(polyin_work);end
    critPT =[];splitEdge=[];
    global MCD_CAP_ON MCD_CAP
    for p = 1:length(polyin_buffed)
        style='r*';nodeend=fliplr(nodeend);t=0.5;
        boundary=rmholes(polyin_buffed(p));
        spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);
        critP =[];cvex=[];
        % Smooth Polyshape Edge
        newBoundary = [];
       
        boundary = refinePoly(boundary,2,flag);
        [newBoundary(:,1),newBoundary(:,2)] = boundary.boundary;
        if ~isempty(MCD_CAP_ON); MCD_CAP{end+1}={'bnd_raw',p,newBoundary}; end

        % Forward boundary scan
        TF=islocalmin(round(newBoundary,5),'FlatSelection', 'center');sN=find(TF(:,1));
        while isempty(sN)
            clear newBoundary
            boundary = refinePoly(boundary,1,flag);[newBoundary(:,1),newBoundary(:,2)] = boundary.boundary;
            TF=islocalmin(round(newBoundary,5),'FlatSelection', 'center');sN=find(TF(:,1));
        end
        sN=sN(1);newBoundary=circshift(newBoundary,sN);
        TF=islocalmin(round(newBoundary,5),'FlatSelection', 'center');
        TF([find(isnan(newBoundary(:,1)))+1 ; find(isnan(newBoundary(:,1)))-1],:)=0;
        subcritP=newBoundary(TF(:,1),:);inPoly=polyshape(newBoundary);%inPoly.plot
        
        if isempty(subcritP)    % If first Critcal point is empty, try without smoothing. 
            critP =[];cvex=[];
            % Smooth Polyshape Edge
            newBoundary = [];
            boundary=rmholes(polyin_buffed(p));
            boundary = refinePoly(boundary,2,~flag);
            [newBoundary(:,1),newBoundary(:,2)] = boundary.boundary;

            % Forward boundary scan
            TF=islocalmin(round(newBoundary,5),'FlatSelection', 'center');sN=find(TF(:,1));
            while isempty(sN)
                clear newBoundary
                boundary = refinePoly(boundary,1,~flag);[newBoundary(:,1),newBoundary(:,2)] = boundary.boundary;
                TF=islocalmin(round(newBoundary,5),'FlatSelection', 'center');sN=find(TF(:,1));
            end
            sN=sN(1);newBoundary=circshift(newBoundary,sN);
            TF=islocalmin(round(newBoundary,5),'FlatSelection', 'center');
            TF([find(isnan(newBoundary(:,1)))+1 ; find(isnan(newBoundary(:,1)))-1],:)=0;
            subcritP=newBoundary(TF(:,1),:);inPoly=polyshape(newBoundary);
        end

        for j=find(TF(:,1))'
            in1 = isinterior(inPoly,[newBoundary(j,1)-5 newBoundary(j,2)]); 
            in2 = isinterior(inPoly,[newBoundary(j,1)+5 newBoundary(j,2)]);
            if in1&&~in2; cvex=[cvex;0];elseif ~in1&&in2 ; cvex=[cvex;1];else; cvex=[cvex;1]; end
        end
        critP = [critP;subcritP];
        if ~isempty(MCD_CAP_ON); MCD_CAP{end+1}={'fwd',p,subcritP,newBoundary(1,:)}; end

        % Holes boundary scan
        obj = holes(polyin_buffed(p));
        obj = polybuffer(obj, 10);%obj.plot
        obj = polybuffer(obj,-10,'JointType','miter','MiterLimit',4);%obj.plot

        poly=obj;

        for i= 1: length(poly)
            [values(:,1),values(:,2)] =poly(i).boundary;%spcrv(poly(i).Vertices',3);values=values';
            if ~isempty(MCD_CAP_ON); MCD_CAP{end+1}={'hbnd',p,values,i}; end
            TF=islocalmin(values(:,1));TF([find(isnan(values(:,1)))+1 ; find(isnan(values(:,1)))-1],:)=0;
            subcritP=values(TF(:,1),:);
            critP = [critP;subcritP];
            if ~isempty(MCD_CAP_ON); MCD_CAP{end+1}={'hmin',p,subcritP,i}; end
            for j=find(TF(:,1))'                            % Checking convex or concave
                in1 = ~isinterior(poly(i),[values(j,1)-5 values(j,2)]);
                in2 = ~isinterior(poly(i),[values(j,1)+5 values(j,2)]);
                if in1&&~in2; cvex=[cvex;0];elseif ~in1&&in2 ; cvex=[cvex;1];else; cvex=[cvex;1]; end
            end
            TF=islocalmax(values);TF([find(isnan(values(:,1)))+1 ; find(isnan(values(:,1)))-1],:)=0;
            subcritP=values(TF(:,1),:);
            critP = [critP;subcritP];
            if ~isempty(MCD_CAP_ON); MCD_CAP{end+1}={'hmax',p,subcritP,i}; end
            for j=find(TF(:,1))'
                in1 = ~isinterior(poly(i),[values(j,1)-5 values(j,2)]);
                in2 = ~isinterior(poly(i),[values(j,1)+5 values(j,2)]);
                if in1&&~in2; cvex=[cvex;1];elseif ~in1&&in2 ; cvex=[cvex;0];else; cvex=[cvex;1]; end
            end
            clear values
        end

        % Backward boundary scan
        TF=islocalmax(round(newBoundary,5),'FlatSelection', 'center');sN=find(TF(:,1));
        if ~(any(TF(:,1)))
           newBoundary=circshift(newBoundary,10);TF=islocalmax(round(newBoundary,5),'FlatSelection', 'center');
        end
        TF([find(isnan(newBoundary(:,1)))+1 ; find(isnan(newBoundary(:,1)))-1],:)=0;
        subcritP=newBoundary(TF(:,1),:);
        critP = [critP;subcritP];
        if ~isempty(MCD_CAP_ON); MCD_CAP{end+1}={'bwd',p,subcritP,newBoundary(1,:)}; end

        for j=find(TF(:,1))'
            in1 = isinterior(inPoly,[newBoundary(j,1)-5 newBoundary(j,2)]);
            in2 = isinterior(inPoly,[newBoundary(j,1)+5 newBoundary(j,2)]);
            if in1&&~in2; cvex=[cvex;1];elseif ~in1&&in2 ; cvex=[cvex;0];else; cvex=[cvex;1]; end
        end
        
        % Decomposing the cells
        ymin = round(min(polyin_work(p).Vertices(:,2)));ymax=round(max(polyin_work(p).Vertices(:,2)));  %%%polyin_buffed
        for c = 1:length(critP(:,1))
                bool=false;boo11=false;boo12=false;
                [in,~] = intersect(polyin_work(p),[critP(c,1),ymin;critP(c,1),ymax]);%in=fliplr(in);out=fliplr(out);    %%%polyin_buffed
                if ~isempty(in); if ymax-in(end,2)<50; in(end,2)=ymax;end; end                              % Keep in mind
                if ~isempty(in); if in(1,2)-ymin<50; in(1,2)=ymin;end; end                                  % Keep in mind
                if cvex(c)==0
                    if sum(isnan(in(:,1)))==0
                        splitEdge = [splitEdge ; polybuffer([in(1,:);in(end,:)],'lines',t)];
                    else
                        % Nudge NaN-adjacent rows to help identify the bounding segment
                        in(find(isnan(in(:,2)))-1,2)=in(find(isnan(in(:,2)))-1,2)-1;
                        in(find(isnan(in(:,2)))+1,2)=in(find(isnan(in(:,2)))+1,2)+1;

                        oo = rmmissing(in);
                        for k = 1:sum(isnan(in(:,1)))+1
                            if oo(2*k-1,2)<critP(c,2) && oo(2*k,2)>critP(c,2)
                                bool=true;%ii=[oo(2*k-1,:);oo(2*k,:)];
                                boo11=true;
                                break;
                            end
                        end
                        if ~bool
                           kk=find(diff((oo(:,2)-critP(c,2))>=0),1);kk=[kk;kk+1];
                           bool=true;
                           boo12=true;
                        end
                        if bool
                            if boo11; splitEdge = [splitEdge ; polybuffer([oo(2*k-1,:);oo(2*k,:)]+[0,-10;0,10],'lines',t)];end
                            if boo12; splitEdge = [splitEdge ; polybuffer([oo(kk(1)-1,:);oo(kk(2)+1,:)]+[0,-10;0,10],'lines',t)];end
                        end
                    end
                end
        end

        distt=[];
        for ii = 1:size(critP,1)
           distt = [distt,spdist(critP(ii,:),critP)];
        end
        distt(distt==0)=inf;
        if min(distt(:))<=25
            A=min(distt(:))==distt;A(find(cvex),:)=[];
            for ii=1:length(A(:,1))
               critP(A(ii,:),:)=[]; 
            end
        end
        critPT=[critPT;critP];
    end
    
    if ~isempty(splitEdge)
        for i=1:length(splitEdge)
            polyin_work=subtract(polyin_work,splitEdge(i));%splitEdge.plot          % Splits the Poly region into different cells, Decompsition 
            polyin=subtract(polyin,splitEdge(i));        
        end
    end
    polyin_work = sortregions(polyin_work,'centroid','ascend');
    polyin = sortregions(polyin,'centroid','ascend');

    if length(polyin_work)>1; polyin_work=regJoin(polyin_work);end
    polyout_work=regions(polyin_work);
    polyout=regions(polyin);

end


function polyout=polyclean(polyin)
    %global s a
    polyout=polyshape();
    if isscalar(polyin)
        poly = regions(polyin);% pp=polybuffer(poly,-s/4);
        poly=poly(fix(poly.area*1e-03)>5);%poly.area>100 10
        % poly=poly(fix(pp.area)>2*a*a*pi);
        polyout = regJoin(poly);
    else
        poly = polyin;% pp=polybuffer(poly,-s/4);
        poly=poly(fix(poly.area*1e-03)>5);%poly.area>100 10
        % poly=poly(fix(pp.area)>2*a*a*pi);
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