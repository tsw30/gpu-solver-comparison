% Initialise case for TS3 solver
clearvars file settings
%% Inputs, settings
% input file location
file.grid_dir = '/home/tsw30/ifan/mesh';
file.grid_file_name = 'LSF5_template_2';
% output location
file.output_dir = '/rds/project/gp10006/rds-gp10006-hall-pullan-ifan/tsw30/hydra/LSF5/MESH/TS3';

% case settings
settings.NRT = 0.95;
settings.rpm_105 = -7224.315;
settings.T_in = 288.15;
settings.p_in = 101325;
settings.mcor_kg_s_des = 110.4711*1.00386;
settings.RPTR_des = 1.46523;
