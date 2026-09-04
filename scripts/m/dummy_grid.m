function g = dummy_grid(plot_blocks, NI, NJ, NK)
% DUMMY_GRID Creates a dummy grid structure for testing purposes.
% Optional NI, NJ, NK (odd integers, default 11) scale the grid, e.g. to
% time ts3_extract on realistic block sizes.

if ~exist('plot_blocks', 'var')
    plot_blocks = false;
end
if ~exist('NI', 'var') || isempty(NI), NI = 11; end
if ~exist('NJ', 'var') || isempty(NJ), NJ = 11; end
if ~exist('NK', 'var') || isempty(NK), NK = 11; end
if mod(NI,2)~=1 || mod(NJ,2)~=1 || mod(NK,2)~=1
    error('dummy_grid: NI, NJ and NK must be odd');
end

% load measurement plane coordinates
% return real script path
script_path = fileparts(mfilename('fullpath'));
coords_file = fullfile(script_path,'../measurement_coords.json');
if ~exist(coords_file,'file')
    error('Measurement coordinates file not found: %s', coords_file);
end
% coords data storing (x,r) arrays for 'L00S', 'BPIN', 'ESIN' measurement planes
measurement_coords = jsondecode(fileread(coords_file));

NBLADE = 20;
PERIOD = 2*pi/NBLADE;
R1 = 1; R2 = 1.5;
ROW_R1 = [measurement_coords.L00S(1,2); measurement_coords.ESIN(1,2); measurement_coords.BPIN(1,2)];
ROW_R2 = [measurement_coords.L00S(end,2); measurement_coords.ESIN(end,2); measurement_coords.BPIN(end,2)];
ROW_SCALE = (ROW_R2 - ROW_R1)/(R2-R1);
DX = 1;
% offset each row so its rake leading point sits mid-way along the inlet block
ROW_DX = [measurement_coords.L00S(1,1); measurement_coords.ESIN(1,1); measurement_coords.BPIN(1,1)] - 0.5*DX*ROW_SCALE;
ROW_DR = ROW_R1 - R1*ROW_SCALE;
% 9 block grid, 3 rows containing inlet, blade o-mesh, and outlet blocks
cp = 1005;
ga = 1.4;
R = cp*(ga-1)/ga;
cv = cp/ga;
T_ref = 288.15;
p_ref = 101325;
ro_ref = p_ref/(R*T_ref);
v_ref = [10 0 0]; % reference velocity
M = v_rcpt0_return(sqrt(sum(v_ref.^2))/sqrt(cp*T_ref), ga);
ro = ro_ref / ro_eval(M, ga);
T = T_ref / (1 + (ga-1)/2*M^2);
roe = ro *(cv*T + 0.5*sum(v_ref.^2));

% inlet block (same coordinates for all rows with x-offset), simple annular grid
ni = NI; nj = NJ; nk = NK;
x = linspace(0, DX, ni);
r = linspace(R1, R2, nj);
theta = linspace(0, PERIOD, nk);
block_in = struct();
block_in.x = repmat(x', 1, nj, nk);
block_in.r = repmat(r, ni, 1, nk);
theta = repmat(reshape(theta, 1, 1, nk), ni, nj, 1);
block_in.rt = block_in.r .* theta;
block_in.attribute.bid = 0; block_in.attribute.np = 0;
block_in.attribute.ni = ni; block_in.attribute.nj = nj; block_in.attribute.nk = nk;
block_in = default_avs_bvs(block_in, cp, ga, NBLADE);
block_in = set_flow_vars(block_in, [ni,nj,nk], ro, v_ref, roe);
block_in.patch = {};
block_in = add_self_periodic(block_in, 0, 1);

g{1} = block_in;

% o-mesh block
ni_side = NI; % number of points on one side of the o-mesh
ni = 4*(ni_side-1) + 1; % o-mesh wraps around in i-direction
nk = (NK-1)/2 + 1; % leave internal face
block_o = struct();
dx_inner = 0.1;
dt_inner = PERIOD/5;
x_outer = DX+DX*[ones((ni_side-1)/2,1) ; linspace(1, 0, ni_side)'; zeros(ni_side-2,1); linspace(0, 1, ni_side)'; ones((ni_side-1)/2,1)];
t_outer = PERIOD/2 - PERIOD/2 * [linspace(0, 1, (ni_side-1)/2 + 1)'; ones(ni_side-2,1); linspace(1, -1, ni_side)'; -ones(ni_side-2,1); linspace(-1, 0, (ni_side-1)/2 + 1)'];
x_inner = DX+DX*[ones((ni_side-1)/2,1) - dx_inner ; linspace(1-dx_inner, dx_inner, ni_side)'; ones(ni_side-2,1)*dx_inner; linspace(dx_inner, 1-dx_inner, ni_side)'; ones((ni_side-1)/2,1) - dx_inner];
t_inner = PERIOD/2 - (PERIOD/2 - dt_inner) * [linspace(0, 1, (ni_side-1)/2 + 1)'; ones(ni_side-2,1); linspace(1, -1, ni_side)'; -ones(ni_side-2,1); linspace(-1, 0, (ni_side-1)/2 + 1)'];

% linearly interpolate between inner and outer surfaces
block_o.x = repmat(x_inner, 1, nj, nk) + (repmat(x_outer, 1, nj, nk) - repmat(x_inner, 1, nj, nk)) .* repmat(reshape(linspace(0, 1, nk), 1, 1, nk), ni, nj, 1);
block_o.r = repmat(linspace(R1, R2, nj), ni, 1, nk);
theta = repmat(t_inner, 1, nj, nk) + (repmat(t_outer, 1, nj, nk) - repmat(t_inner, 1, nj, nk)) .* repmat(reshape(linspace(0, 1, nk), 1, 1, nk), ni, nj, 1);
block_o.rt = block_o.r .* theta;
block_o.attribute.bid = 1; block_o.attribute.np = 0;
block_o.attribute.ni = ni; block_o.attribute.nj = nj; block_o.attribute.nk = nk;
block_o = default_avs_bvs(block_o, cp, ga, NBLADE);
block_o = set_flow_vars(block_o, [ni,nj,nk], ro, v_ref, roe);

block_o.patch = {};
% patch front face to inlet block rear face
p1 = struct(); p2 = struct();
p1.attribute.kind=5; p2.attribute.kind = 5;
p1.attribute.bid = 0; p1.attribute.pid = length(g{1}.patch);
p1.attribute.nxbid = 1; p1.attribute.nxpid = length(block_o.patch);
p1.attribute.ist = g{1}.attribute.ni-1; p1.attribute.ien = g{1}.attribute.ni;
p1.attribute.jst = 0; p1.attribute.jen = g{1}.attribute.nj;
p1.attribute.kst = 0; p1.attribute.ken = g{1}.attribute.nk;
p1.attribute.idir = 6; p1.attribute.jdir = 1; p1.attribute.kdir = 0;
p2.attribute.bid = 1; p2.attribute.pid = length(block_o.patch);
p2.attribute.nxbid = 0; p2.attribute.nxpid = length(g{1}.patch);
p2.attribute.ist = (ni_side-1)/2 + ni_side-1; p2.attribute.ien = (ni_side-1)/2 + 2*(ni_side-1)+1;
p2.attribute.jst = 0; p2.attribute.jen = block_o.attribute.nj;
p2.attribute.kst = block_o.attribute.nk-1; p2.attribute.ken = block_o.attribute.nk;
p2.attribute.idir = 2; p2.attribute.jdir = 1; p2.attribute.kdir = 6;
g{1}.patch{p1.attribute.pid+1} = p1; g{1}.attribute.np = length(g{1}.patch);
block_o.patch{p2.attribute.pid+1} = p2; block_o.np = length(block_o.patch);
% self periodic 
block_o = add_self_periodic(block_o,1,2);
% i=0,i=ni self periodic
p1 = struct(); p1.attribute.kind = 5; p1.attribute.bid = 1; p1.attribute.nxbid = 1;
p1.attribute.jst =0; p1.attribute.jen = nj; p1.attribute.kst = 0; p1.attribute.ken = nk;
p1.attribute.idir = 6; p1.attribute.jdir =1; p1.attribute.kdir = 2;
p2 = p1;
p1.attribute.pid = length(block_o.patch); p2.attribute.pid = p1.attribute.pid+1;
p1.attribute.nxpid = p2.attribute.pid; p2.attribute.nxpid = p1.attribute.pid;
p1.attribute.ist = 0; p1.attribute.ien = 1;
p2.attribute.ist = ni-1; p2.attribute.ien = ni;
block_o.patch{p1.attribute.pid+1} = p1; block_o.patch{p2.attribute.pid+1} = p2;
block_o.attribute.np = length(block_o.patch);

g{length(g)+1} = block_o;

% outlet block
ni = NI; nj = NJ; nk = NK;
block_out = struct();
block_out.x = block_in.x + 2*DX;
block_out.r = block_in.r;
block_out.rt = block_in.rt;
block_out.attribute.bid = 2; block_out.attribute.np=0;
block_out.attribute.ni = ni; block_out.attribute.nj = nj; block_out.attribute.nk = nk;
block_out = default_avs_bvs(block_out, cp, ga, NBLADE);
block_out = set_flow_vars(block_out, [ni,nj,nk], ro, v_ref, roe);

block_out.patch = {};
% patch o-mesh rear face(s) to outlet front face
p1 = struct;
p1.attribute.kind = 5;
p1.attribute.jst = 0; p1.attribute.jen = nj; p1.attribute.jdir = 1; 
p1a = p1; p1a.attribute.nxbid = 2;
p1a.attribute.bid = 1; p1a.attribute.idir = 5; p1a.attribute.kdir = 6;
p1a.attribute.kst = g{2}.attribute.nk-1; p1a.attribute.ken = g{2}.attribute.nk;
p1b = p1a;
p1a.attribute.pid = length(g{2}.patch); p1b.attribute.pid = p1a.attribute.pid+1;
p1a.attribute.nxpid = length(block_out.patch); p1b.attribute.nxpid = p1a.attribute.nxpid+1;
p1a.attribute.ist = 0; p1a.attribute.ien = (ni_side-1)/2 + 1;
p1b.attribute.ist = (ni_side-1)/2+3*(ni_side-1); p1b.attribute.ien = g{2}.attribute.ni;
p2a = p1; p2a.attribute.nxbid = 1;
p2a.attribute.bid = 2; p2a.attribute.idir = 6; p2a.attribute.kdir = 3;
p2a.attribute.ist = 0; p2a.attribute.ien = 1;
p2b = p2a;
p2a.attribute.pid = length(block_out.patch); p2b.attribute.pid = p1a.attribute.nxpid+1;
p2a.attribute.nxpid = length(g{2}.patch); p2b.attribute.nxpid = p1a.attribute.pid+1;
p2a.attribute.kst = (ni_side-1)/2; p2a.attribute.ken = block_out.attribute.nk;
p2b.attribute.kst = 0; p2b.attribute.ken = (ni_side-1)/2 + 1;
g{2}.patch{p1a.attribute.pid+1} = p1a;
g{2}.patch{p1b.attribute.pid+1} = p1b;
g{2}.attribute.np = length(g{2}.patch);
block_out.patch{p2a.attribute.pid+1} = p2a;
block_out.patch{p2b.attribute.pid+1} = p2b;
block_out.attribute.np = length(block_out.patch);
% self periodics
block_out = add_self_periodic(block_out, 2, 1);
g{length(g)+1} = block_out;

% Ignore outlet plane (default set as wall) - rows left uncoupled

% scale and add 3x rows, offsetting bid references recursively
geo_vars = {'x','r','rt'};
g_out = cell(9,1);
for row = 1:3
    bid_offset = (row-1)*3;
    for block = 1:3
        bid_lo = block-1 + bid_offset;
        g_out{bid_lo + 1} = g{block};
        for v = 1:length(geo_vars)
            var = geo_vars{v};
            g_out{bid_lo + 1}.(var) = g_out{bid_lo + 1}.(var)*ROW_SCALE(row);
        end
        g_out{bid_lo + 1}.x = g_out{bid_lo + 1}.x + ROW_DX(row);
        theta = g_out{bid_lo + 1}.rt./g_out{bid_lo+1}.r;
        g_out{bid_lo+1}.r = g_out{bid_lo+1}.r + ROW_DR(row);
        g_out{bid_lo+1}.rt = g_out{bid_lo+1}.r.*theta;
        g_out{bid_lo + 1}.attribute.bid = bid_lo;
        for ip = 1:g_out{bid_lo + 1}.attribute.np
            if isfield(g_out{bid_lo + 1}.patch{ip}.attribute,'nxbid')
                g_out{bid_lo + 1}.patch{ip}.attribute.nxbid = g_out{bid_lo + 1}.patch{ip}.attribute.nxbid + bid_offset;
            end
            g_out{bid_lo+1}.patch{ip}.attribute.bid = g_out{bid_lo+1}.patch{ip}.attribute.bid + bid_offset;
        end
    end
    if row == 1 % add inlet patch after scaling (patches carry no geometry)
        g_out{1} = add_inlet_patch(g_out{1}, 0, struct('pstag', p_ref, 'tstag', T_ref));
    end
end

g = g_out;

if plot_blocks
    figure;
    for b = 1:length(g)
        subplot(3,length(g)/3,b);
        block = g{b};
        % surf all 6 faces
        ni = block.attribute.ni; nj = block.attribute.nj; nk = block.attribute.nk;
        theta = block.rt./block.r;
        y = block.r.*cos(theta); z = block.r.*sin(theta);
        % face 1 (i=0)
        surf(squeeze(block.x(1,:,:)), squeeze(z(1,:,:)), squeeze(y(1,:,:)),'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
        hold on;
        % face 2 (i=ni-1)
        surf(squeeze(block.x(ni,:,:)), squeeze(z(ni,:,:)), squeeze(y(ni,:,:)),'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
        % face 3 (j=0)
        surf(squeeze(block.x(:,1,:)), squeeze(z(:,1,:)), squeeze(y(:,1,:)),'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
        % face 4 (j=nj-1)
        surf(squeeze(block.x(:,nj,:)), squeeze(z(:,nj,:)), squeeze(y(:,nj,:)),'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
        % face 5 (k=0)
        surf(squeeze(block.x(:,:,1)), squeeze(z(:,:,1)), squeeze(y(:,:,1)),'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
        % face 6 (k=nk-1)
        surf(squeeze(block.x(:,:,nk)), squeeze(z(:,:,nk)), squeeze(y(:,:,nk)),'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0,0,0]);
        axis equal;
        xlabel('X'); ylabel('Z'); zlabel('Y');
        title(['Block ' num2str(b)]);
    end
end
end

%% Helper functions

function M = v_rcpt0_return(v_rcpt0, ga)
% Computes the Mach number from the velocity and stagnation temperature.
V2 = v_rcpt0.*v_rcpt0;
M2 = V2/(ga-1)./(1-V2/2);
M = sqrt(M2);
end

function ro_rat = ro_eval(M, ga)
% Computes the density from the Mach number and stagnation temperature.
ro_rat = exp(1/(ga-1)*log(1+(ga-1)/2.*M.*M));
end

function self = default_avs_bvs(self,cp,ga,nblade)
% Set default avs and bvs for block
self.av.cp = cp; self.av.ga = ga; self.av.viscosity = 1.8e-5; self.av.viscosity_law = 0;
self.av.ilos = 1; % mixing length turbulence model
self.bv.rpm = 0; self.bv.nblade = nblade; self.bv.fblade = nblade;
self.bv.rpmi1 = 0; self.bv.rpmi2 = 0; 
self.bv.rpmj1 = 0; self.bv.rpmj2 = 0;
self.bv.rpmk1 = 0; self.bv.rpmk2 = 0;
end

function self = set_flow_vars(self, bshape, ro, v_ref, roe)
% Set uniform flow variables
ni = bshape(1); nj = bshape(2); nk = bshape(3);
self.ro = ones(ni, nj, nk) * ro; % uniform density
self.rovx = ones(ni, nj, nk) * ro * v_ref(1); % uniform velocity in x
self.rovr = ones(ni, nj, nk) * ro * v_ref(2); % uniform velocity in r
self.rorvt = ones(ni, nj, nk) * ro * v_ref(3); % uniform velocity in theta
self.roe = ones(ni, nj, nk) * roe; % uniform energy
end


function self = add_inlet_patch(self, bid, inlet_conds)
% Add inlet patch to block
nj = self.attribute.nj; nk = self.attribute.nk;
pid = length(self.patch);
p.attribute.bid = bid;
p.attribute.pid = pid;
p.attribute.ist = 0; p.attribute.ien = 1;
p.attribute.jst = 0; p.attribute.jen = nj;
p.attribute.kst = 0; p.attribute.ken = nk;
p.attribute.kind = 0; % inlet patch
p.pstag = ones(nj*nk,1) * inlet_conds.pstag; % uniform stagnation pressure
p.tstag = ones(nj*nk,1) * inlet_conds.tstag; % uniform stagnation temperature
p.pitch = zeros(nj*nk,1); % uniform pitch
p.yaw = zeros(nj*nk,1); % uniform yaw

self.patch{pid+1} = p;
self.attribute.np = length(self.patch);
end

function self = add_self_periodic(self, bid, block_kind)
% Add self-periodic patch to block

ni = self.attribute.ni; nj = self.attribute.nj; nk = self.attribute.nk;
pid = length(self.patch);
p.attribute.bid = bid;
p.attribute.pid = pid;
p.attribute.nxbid = bid;
p.attribute.nxpid = pid+1;
p.attribute.kind = 5;
p.attribute.jst = 0; p.attribute.jen = nj;
if block_kind == 2 % o-mesh block, wraps around in i-direction
    ni_side = (ni-1) / 4 + 1;
    p.attribute.ist = (ni_side-1)/2; p.attribute.ien = (ni_side-1)/2 + ni_side;
    p.attribute.kst = nk-1; p.attribute.ken = nk;
    p.attribute.idir = 3; p.attribute.jdir = 1; p.attribute.kdir = 6;
else
    p.attribute.ist = 0; p.attribute.ien = ni;
    p.attribute.kst = 0; p.attribute.ken = 1;
    p.attribute.idir = 0; p.attribute.jdir = 1; p.attribute.kdir = 6;
end
self.patch{pid+1} = p;
% 2nd patch
p.attribute.pid = pid+1;
p.attribute.nxpid = pid;
if block_kind == 2 % o-mesh block, wraps around in i-direction
    p.attribute.ist = (ni_side-1)/2 + 2*(ni_side-1); p.attribute.ien = (ni_side-1)/2 + 3*(ni_side-1)+1;
else
    p.attribute.kst = nk-1; p.attribute.ken = nk;
end
self.patch{pid+2} = p;
self.attribute.np = length(self.patch);
end