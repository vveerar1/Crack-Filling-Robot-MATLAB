% endP_ident  Identify endpoint types in a binary crack skeleton.
% Author: Vishnu Veeraraghavan,
% Automated Control Systems and Robotics Lab.
% Email: vveerar1@binghamton.edu.
% July 2019, Last Revision: 25-Sep-2019

function [eP,rP,cP] = endP_ident(BW3,BW,org)
% endP_ident  Identify endpoint types in a binary skeleton.
%   [eP,rP,cP] = endP_ident(BW3,BW,org)
%   Returns all skeleton endpoints eP, the subset that are real (free-end)
%   endpoints rP, and the continuing (boundary-adjacent) endpoints cP.

    if sum(sum(BW3))>5

        spdist = @(P,Ps) sqrt((P(1,1)-Ps(:,1)).^2 + (P(1,2)-Ps(:,2)).^2);

        if ~isempty(BW)
            BW3=padarray(BW3,[0,1],'post');
            BW3=padarray(BW3,[1,0],'post');
        end

        [~,endP]=neighbor_count_points(BW3,1);  %1: End Points

        endPoints = endP;post=[]; eP=endP;

        %%%%
        %%%%
        I=BW;I_te=I;
        while true
            intPoints=neighbor_count_points(I,2); %2: Branchpoints
            I_te(intPoints>0)=0;
            pixelP=neighbor_count_points(I_te,0); %0: Single Pixels
            I_te(pixelP>0)=0;

            if ~isequal(I,I_te)
               I=I_te;
            else
               break;
            end
        end
        BW=I;
        %%%%
        %%%%

        if ~isempty(BW)
            I = BW;
            for i = 1:length(endPoints(:,1))
                ePp=endPoints(i,:);
                row=ePp(1);
                col=ePp(2);
                b = [I(row-1,col-1:col+1) ...
                        I(row, col-1) I(row, col+1)...
                        I(row+1,col-1:col+1)];
                post = [post;b];
            end

            post=sum(post,2);
            realEndP = endPoints(post~=2,:);rP=realEndP;
            contEndP = endPoints(post==2,:);cP=contEndP;
        else
            rP=[];
            cP=[];
        end

    else
        eP=[];rP=[];cP=[];
    end
