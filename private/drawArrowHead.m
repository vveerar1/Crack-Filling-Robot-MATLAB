function hArrow = drawArrowHead(p0,p1,color)
% drawArrowHead  Draw only the arrowhead (filled triangle) at the tip p1.
%   hArrow = drawArrowHead(p0,p1,color)
%   Renders a filled arrowhead patch pointing from p0 toward p1.
%   color is optional (default 'k').

if nargin == 2
   color = 'k';
end
% Parameters:
W1 = 0.08;   % half width of the arrow head, normalized by length of arrow
W2 = 0.014;  % half width of the arrow shaft
L1 = 0.18;   % Length of the arrow head, normalized by length of arrow
L2 = 0.13;  % Length of the arrow inset
% Unpack the tail and tip of the arrow
x0 = p0(1);
y0 = p0(2);
x1 = p1(1);
y1 = p1(2);
% Draw arrowhead only (no shaft), centred on the tip
P = [...
    (L1-L2), (L1-L1), L1, (L1-L1), (L1-L2);
    W2,     W1, 0,    -W1,    -W2];
P(1,:)=P(1,:)-L1;
% Scale,rotate, shift and plot:
dx = x1-x0;
dy = y1-y0;
Length = sqrt(dx*dx + dy*dy);
Angle = atan2(-dy,dx);
P = 342*P;   %Scale
P = [cos(Angle), sin(Angle); -sin(Angle), cos(Angle)]*P;  %Rotate
P = p1(:)*ones(1,5) + P;  %Shift
% Plot!
hArrow = patch('Faces',1:5,'Vertices',P','FaceColor',color,'EdgeColor',color);  axis equal;
end
