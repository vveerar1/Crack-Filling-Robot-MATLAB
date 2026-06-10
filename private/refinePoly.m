% refinePoly  Refine and upsample a polygon boundary by midpoint insertion.
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% October 2019, Last Revision: 20-October-2019

function polyout = refinePoly(polyin,times,smoothflag)
% refinePoly  Refine and upsample a polygon boundary.
%   polyout = refinePoly(polyin,times,smoothflag)
%   Inserts midpoints along every edge of polyin for 'times' passes.
%   If smoothflag is true, applies a light polybuffer smoothing pass as well.

    switch nargin
        case 2
            smoothflag=false;
    end
    polyout=polyshape();

    for n = 1:times
        for p = 1:length(polyin)
            boundary=rmholes(polyin(p));
            jj=regions(boundary);
            newBoundary = [];

            for i = 1:length(jj)
                if ~isempty(regions(jj(i)))
                    if smoothflag
                        jj(i)=polybuffer(polybuffer(jj(i),-2),2);
                        jj(i) = rmslivers(jj(i),0.00001);
                        x=jj(i).Vertices(:,1);
                        y=jj(i).Vertices(:,2);
                        smoothX = smooth(x);smoothX = smooth(smoothX);
                        smoothY = smooth(y);smoothY = smooth(smoothY);
                        jj(i) =polyshape([smoothX,smoothY]);
                        if length(regions(jj(i)))>1
                            jj(i)=polyshape([x,y]);
                        end
                    end
                    [Bound(:,1),Bound(:,2)] =jj(i).boundary;Bound(end,:)=[];
                    [values(:,1),values(:,2)] =jj(i).boundary;values(end,:)=[];
                    ii=0;
                    for j = 1:length(values(:,1))
                            P1=round(values(j,:),5);
                            P2=round(values(mod(j,length(values(:,1)))+1,:),5);
                            mid=(P1(:)+ P2(:)).'/2;
                            Bound=[Bound(1:j+ii,:) ;mid; Bound(1+j+ii:end,:)];
                            ii=ii+1;
                    end

                    if i==1
                        newBoundary = [newBoundary; Bound(end,:);Bound];clear values Bound;%
                        newBoundary=circshift(newBoundary,fix(size(newBoundary,1)/10));
                    else
                        newBoundary = [newBoundary; NaN,NaN;Bound(end,:);Bound];clear values Bound;%
                        newBoundary=circshift(newBoundary,fix(size(newBoundary,1)/10));
                    end
                end
            end

            ext = polyshape(newBoundary,'Simplify',false);

            jj=holes(polyin(p));
            newBoundary = [];

            for i = 1:length(jj)
                if smoothflag
                    jj(i)=polybuffer(polybuffer(jj(i),-1),1);
                    jj(i) = rmslivers(jj(i),0.00001);
                    x=jj(i).Vertices(:,1);
                    y=jj(i).Vertices(:,2);
                    smoothX = smooth(x);smoothX = smooth(smoothX);
                    smoothY = smooth(y);smoothY = smooth(smoothY);
                    jj(i) =polyshape([smoothX,smoothY]);
                end
                [Bound(:,1),Bound(:,2)] =jj(i).boundary;Bound(end,:)=[];
                [values(:,1),values(:,2)] =jj(i).boundary;values(end,:)=[];
                ii=0;
                for j = 1:length(values(:,1))
                        P1=round(values(j,:),5);
                        P2=round(values(mod(j,length(values(:,1)))+1,:),5);
                        mid=(P1(:)+ P2(:)).'/2;
                        Bound=[Bound(1:j+ii,:) ;mid; Bound(1+j+ii:end,:)];
                        ii=ii+1;
                end

                if i==1
                    newBoundary = [newBoundary; Bound(end,:);Bound];clear values Bound;%
                    newBoundary=circshift(newBoundary,fix(size(newBoundary,1)/10));
                else
                    newBoundary = [newBoundary; nan,nan;Bound(end,:);Bound];clear values Bound;%
                    newBoundary=circshift(newBoundary,fix(size(newBoundary,1)/10));
                end
            end

            if ~isempty(newBoundary); hol = polyshape(newBoundary,'Simplify',false);else;hol = polyshape(); end

            subpolyout = subtract(ext,hol,'KeepCollinearPoints',true);

            if length(polyout)>1
                polyout=[polyout;subpolyout];
            end
        end

        if length(polyin)==1
            polyout=subpolyout;
        end
        polyin=polyout;
    end
end
