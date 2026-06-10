% compCrack  Trace a crack skeleton into ordered polylines.
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [crackRaw,line,pointX,pointY] = compCrack(I,endP,dir_map,colors)
% compCrack  Trace a binary crack skeleton into ordered polylines.
%   [crackRaw,line,pointX,pointY] = compCrack(I,endP,dir_map,colors)
%   Walks the binary skeleton I starting from the endpoint list endP and
%   produces crackRaw (cell array of [row,col] polylines), line (Nx4
%   start/end pairs), and pointX/pointY (Nx24 sampled coordinate arrays).

    line=[];
    pointX=[];
    pointY=[];
    i=1;spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);
    endcheck=0;circc=false;
    tempXX1=[];
    tempYY1=[];

    while ~isempty(endP)
        start=endP(1,:);    % Selecting the first row of elements
        row=start(1);
        col=start(2);
        b = [I(row-1,col-1:col+1) ...
                I(row, col-1) I(row, col+1)...
                I(row+1,col-1:col+1)];
        number=sum(b);                         % adds all the elements in b
        if ~(number>=2)
            endP(1,:)=[];       % Removing the first row of elements
            I(row,col)=0;

        else
            endP=circshift(endP,-1);
        end

        tempX=[];
        tempY=[];
        crack=[];
        while (1)

            tempX=[tempX row];
            tempY=[tempY col];
            b = [I(row-1,col-1:col+1) ...
                I(row, col-1) I(row, col+1)...
                I(row+1,col-1:col+1)];             % b captures the 8 surounding nodes of the end point
            number=sum(b);                         % adds all the elements in b
            if number>=2 % intersection point

                %%% Intersection Point Check
                I_c=I;
                deleteInd=[];
                [nextall]=find(b);                 % Finds which elements in b has value 1
                for j=1:length(nextall)

                    rowb=row+dir_map(nextall(j),1);
                    colb=col+dir_map(nextall(j),2);
                    I_c(rowb,colb)=0;
                end
                count=length(nextall);
                for j=1:count %check false intersection point
                    dir=nextall(j);
                    rowb=row+dir_map(dir,1);
                    colb=col+dir_map(dir,2);
                    bb = [I_c(rowb-1,colb-1:colb+1) ...
                        I_c(rowb, colb-1) I_c(rowb, colb+1)...
                        I_c(rowb+1,colb-1:colb+1)];
                    [indb]=find(bb);
                    if isempty(indb)
                        deleteInd=[deleteInd j];
                        continue;
                    end

                end
                nextall(deleteInd)=[];
                if length(nextall)==1
                    I=I_c;
                    row=row+dir_map(nextall,1);
                    col=col+dir_map(nextall,2);
                    crack=[crack; row col];
                    continue;
                end
                %%%

                if length(tempX)==1
                    i=i+1;
                    break;
                end
                if circc
                    if any(spdist([tempX(end) tempY(end)],[tempXX1' tempYY1'])<2)
                       ttk = spdist([tempX(end) tempY(end)],[tempXX1' tempYY1']);
                       ttt=find(ttk<2 & min(ttk)==ttk);ttt=[ttt ttt+1];
                       tempX=[tempX tempXX1(ttt)];row = tempX(end);
                       tempY=[tempY tempYY1(ttt)];col = tempY(end);
                       circc=false;endcheck=0;
                    end
                end

                if ~isempty(endP)
                    [f,d]=dsearchn(endP,[row col]);
                    if d<10
                        EndPoint = endP(f,:);
                    else
                        EndPoint=[row col];
                    end
                else
                    EndPoint=[row col];
                end
                I(EndPoint(1),EndPoint(2))=1;
                if all(EndPoint==start)||spdist(EndPoint,start)<2
                    i=i+1;
                    break;
                end
                b = [I(row-1,col-1:col+1) ...
                        I(row, col-1) I(row, col+1)...
                        I(row+1,col-1:col+1)];
                number=sum(b);                         % adds all the elements in b
                if ~(number>=2)
                    endP(1,:)=[];       % Removing the first row of elements
                end

                i=i+1;
                break;

            end
            if number==1 % on going pt
                if ~isempty(endP)&&~isequal(start,[row col])
                    f=sum(endP==[row col],2)==2;
                    if any(f)
                        endP(f,:)=[];
                    end
                end

                    [ind]=find(b);
                    dir=ind;
                    row=row+dir_map(ind,1);
                    col=col+dir_map(ind,2);
                    crack=[crack; row col];
                    I(row,col)=0;

            end
            if number==0 % end point
                if circc
                    if any(spdist([tempX(end) tempY(end)],[tempXX1' tempYY1'])<2)
                       ttk=spdist([tempX(end) tempY(end)],[tempXX1' tempYY1']);
                       ttt=find(ttk<2 & min(ttk)==ttk);ttt=[ttt ttt+1];
                       tempX=[tempX tempXX1(ttt)];row = tempX(end);
                       tempY=[tempY tempYY1(ttt)];col = tempY(end);
                       circc=false;endcheck=0;
                    end
                end

                if ~isempty(endP)
                    [f,d]=dsearchn(endP,[row col]);
                    if d<10
                        EndPoint = endP(f,:);
                    else
                        EndPoint=[row col];
                    end
                else
                    EndPoint=[row col];
                end

                if all(EndPoint==start)||spdist(EndPoint,start)<2
                    i=i+1;
                    break;
                end
                rowb=EndPoint(1);colb=EndPoint(2);
                b = [I(rowb-1,colb-1:colb+1) ...
                        I(rowb, colb-1) I(rowb, colb+1)...
                        I(rowb+1,colb-1:colb+1)];
                ind=sum(endP==EndPoint,2)==2;
                if any(ind) && sum(b)==0
                    endP(ind,:)=[];
                end

                i=i+1;
                break;
            end

        end
        if ~isempty(crack)
            crackRaw{i-1}=[start;crack;EndPoint];
        else
            i=i-1;
        end

    end

    if exist('crackRaw','var')
        for o = 1:length(crackRaw)
            line=[line;[crackRaw{o}(1,:) crackRaw{o}(end,:)]];
            tempX=crackRaw{o}(:,1);tempY=crackRaw{o}(:,2);
            YI=tempY(round(linspace(1,length(tempY),24)));
            XI=tempX(round(linspace(1,length(tempY),24)));
            if ~isequal([XI(1),YI(1)],[tempX(1),tempY(1)]);XI(1)=tempX(1);YI(1)=tempY(1);end
            if ~isequal([XI(end),YI(end)],[tempX(end),tempY(end)]);XI(end)=tempX(end);YI(end)=tempY(end);end
            pointX(o,:)=XI;
            pointY(o,:)=YI;
        end
    else
        [crackRaw,line,pointX,pointY]=deal([]);
    end

end
