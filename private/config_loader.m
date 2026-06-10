function cfg = config_loader(robot)
% config_loader  Robot configuration loader for the SCC / OnlineSCC planners.
%   cfg = config_loader(robot) returns the dimensions and working radii
%   (pixels) for ROBOT ('robot1', ...), read from matlab/robot_config.json.
%   The JSON file is shared with ROS consumers and is the single source of
%   truth for robot geometry; edit it to change the robot without touching code.
%
%   cfg fields:
%     base_diameter_in       mobile base diameter (in)
%     footprint_diameter_in  nozzle / fill footprint diameter (in)
%     sensor_diameter_in     360-degree sensor diameter (in)
%     r1   base radius (px)        a   footprint radius (px)        s   sensor radius (px)
    if nargin < 1, robot = 'robot1'; end
    here = fileparts(mfilename('fullpath'));               % matlab/private
    raw  = jsondecode(fileread(fullfile(here, '..', 'robot_config.json')));  % matlab/robot_config.json
    d    = raw.(robot);
    inpxMap = @(x) fix(((x)*25.4)/2);                      % inches -> pixels (1 px ~= 2 mm)
    cfg.base_diameter_in      = d.base_diameter_in;
    cfg.footprint_diameter_in = d.footprint_diameter_in;
    cfg.sensor_diameter_in    = d.sensor_diameter_in;
    cfg.r1 = inpxMap(d.base_diameter_in/2);               % base radius
    cfg.a  = inpxMap(d.footprint_diameter_in/2);          % footprint radius
    cfg.s  = inpxMap(d.sensor_diameter_in/2);             % sensor radius
end
