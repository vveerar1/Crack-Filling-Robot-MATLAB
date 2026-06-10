% compCrack_branching  Trace a crack skeleton into polylines, spawning at branch points.
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [crackRaw,line,pointX,pointY] = compCrack_branching(I,endP,dir_map,colors)
% compCrack_branching  Trace a binary crack skeleton into polylines.
%   [crackRaw,line,pointX,pointY] = compCrack_branching(I,endP,dir_map,colors)
%   Walks skeleton I from each endpoint; at each branch point the current
%   polyline is closed and the branch neighbors are queued as new start points.

    line=[];
    pointX=[];
    pointY=[];
    i=1;
    while ~isempty(endP)
        start=endP(1,:);    % Selecting the first row of elements
        row=start(1);
        col=start(2);
        endP(1,:)=[];       % Removing the first row of elements
        I(row,col)=0;
        tempX=[];
        tempY=[];
        crack=[];
        while (1)
            skp=true;
            tempX=[tempX row];
            tempY=[tempY col];
            b = [I(row-1,col-1:col+1) ...
                I(row, col-1) I(row, col+1)...
                I(row+1,col-1:col+1)];             % b captures the 8 surounding nodes of the end point
            number=sum(b);                         % adds all the elements in b
            if number>=2 % intersection point
                deleteInd=[];
                [nextall]=find(b);                 % Finds which elements in b has value 1
                for j=1:length(nextall)

                    rowb=row+dir_map(nextall(j),1);
                    colb=col+dir_map(nextall(j),2);
                    I(rowb,colb)=0;
                end
                count=length(nextall);
                for j=1:count %check false intersection point
                    dir=nextall(j);
                    rowb=row+dir_map(dir,1);
                    colb=col+dir_map(dir,2);
                    bb = [I(rowb-1,colb-1:colb+1) ...
                        I(rowb, colb-1) I(rowb, colb+1)...
                        I(rowb+1,colb-1:colb+1)];
                    [indb]=find(bb);
                    if isempty(indb)
                        deleteInd=[deleteInd j];
                        continue;
                    end

                end
                nextall(deleteInd)=[];
                if length(nextall)==1

                    row=row+dir_map(nextall,1);
                    col=col+dir_map(nextall,2);
                    crack=[crack; row col];
                    continue;
                end
                EndPoint=[row col];
                for j=1:length(nextall)
                    sP=size(endP,1);
                    nextrow=row+dir_map(nextall(j),1);
                    nextcol=col+dir_map(nextall(j),2);
                    [~,ind]=ismember([nextrow nextcol],endP,'rows');
                    if ind
                        continue;
                    end
                    endP(sP+1,:)=[nextrow nextcol];
                end


                if all(EndPoint==start)||abs(sum(EndPoint-start))<2
                    skp=false;
                    break;
                end

                YI=tempY(round(linspace(1,length(tempY),24)));
                XI=tempX(round(linspace(1,length(tempY),24)));
                if ~isequal([XI(1),YI(1)],start);XI(1)=start(1);YI(1)=start(2);end
                if ~isequal([XI(end),YI(end)],EndPoint);XI(end)=EndPoint(1);YI(end)=EndPoint(2);end

                pointX(i,:)=XI;
                pointY(i,:)=YI;
                i=i+1;
                break;

            end
            if number==1 % on going pt
                [ind]=find(b);
                dir=ind;
                row=row+dir_map(ind,1);
                col=col+dir_map(ind,2);
                crack=[crack; row col];
                I(row,col)=0;
            end
            if number==0 % end point
                EndPoint=[row col];
                if all(EndPoint==start)||abs(sum(EndPoint-start))<2
                    skp=false;
                    break;
                end
                [~,ind]=ismember(EndPoint,endP,'rows');
                if ind
                    endP(ind,:)=[];
                end
                YI=tempY(round(linspace(1,length(tempY),24)));
                XI=tempX(round(linspace(1,length(tempY),24)));
                if ~isequal([XI(1),YI(1)],start);XI(1)=start(1);YI(1)=start(2);end
                if ~isequal([XI(end),YI(end)],EndPoint);XI(end)=EndPoint(1);YI(end)=EndPoint(2);end

                pointX(i,:)=XI;
                pointY(i,:)=YI;
                i=i+1;
                break;
            end

        end
        if ~isempty(crack)&&skp
            crackRaw{i-1}=[start;crack;EndPoint];
        end

    end

    for o = 1:length(crackRaw)
        line=[line;[crackRaw{o}(1,:) crackRaw{o}(end,:)]];
    end

end
