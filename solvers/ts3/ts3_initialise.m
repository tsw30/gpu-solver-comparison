% Initialise case for TS3 solver
%   Reads the AutoGrid mesh, applies inter-stage smoothing, matches
%   periodic patches, adds an inlet block, and sets inlet profiles.
%
%   Condensed from the legacy pipeline (ts_autogrid_reader, etc.).
%   Endwall geometry functions (ifan_lsf5_endwalls / ifan_ultrafan_endwalls)
%   are user-provided separately.

clearvars file settings

%% Grid files
% file.grid_dir     = '/home/tsw30/ifan/mesh';
file.grid_dir = 'C:/Users/tsw30/Engineering/Rolls-Royce/iFan/Data/LSF5/GRIDS';
file.grid_file    = 'LSF5_template_4';
file.bcs          = fullfile(file.grid_dir, [file.grid_file '.bcs']);
file.p3d          = fullfile(file.grid_dir, [file.grid_file '.g']);

%% Output and case settings
file.output_dir   = '/rds/project/gp10006/rds-gp10006-hall-pullan-ifan/tsw30/hydra/LSF5/MESH/TS3';
file.inlet_x      = -0.23;      % upstream face of inlet block [m]

settings.NRT          = 0.95;
settings.rpm_105      = -7224.315;
settings.T_in         = 288.15;
settings.p_in         = 101325;
settings.mcor_kg_s_des = 110.4711 * 0.910296; % 110.4711 * 1.00386;
settings.RPTR_des     = 1.37057; % 1.46523;

%% ---- Pipeline ----

% 1. Load grid from AutoGrid files
grid = ts3_load_grid(file.bcs, file.p3d);

% 2. Smooth rotor outlet tangential distribution (needs bypass topology work)
% disp('Smoothing rotor outlet mesh');
% grid = ts3_smooth_interstage(grid);

% 3. Match periodic-patch coordinates
grid = ts3_match_patches(grid);

% 4. Add inlet block (endwall curves from user-provided geometry)
%    Uncomment your chosen endwall function and pass the curves:
%       xrhub = ...; xrcas = ...;  % from ifan_lsf5_endwalls or similar
xrhub = []; xrcas = [];
grid = ts3_interface_block(grid, file.inlet_x, xrhub, xrcas, ...
    'eratio', 1.03, 'straight');

% 5. Apply inlet stagnation profiles matching target mass flow
grid = ts3_inlet(grid, settings.mcor_kg_s_des, settings.p_in, settings.T_in);

fprintf('ts3_initialise: pipeline complete\n');