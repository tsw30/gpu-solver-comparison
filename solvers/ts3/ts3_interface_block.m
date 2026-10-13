function g = ts3_interface_block(g,xinlet,xrhub,xrcas,varargin)
%ts3_interface_block Add an inlet block forward of the first blade row.
%   g = ts3_interface_block(g,xinlet,xrhub,xrcas)
%   g = ts3_interface_block(g,xinlet,xrhub,xrcas,'eratio',1.03,'straight')
%
%   Condensed replacement for ifan_interface_block (original in ref/).
%   Endwall curves xrhub/xrcas are user-provided (e.g. ifan_lsf5_endwalls).
%
%   Parameters:
%     eratio   — hyperbolic bunching expansion ratio (default 1.03)
%     straight — use constant-radius casing (default: from xrcas)
%     frac     — fraction of domain covered by fine inlet cells (default 0.1)
%     hub_cell — target hub cell count fraction (default 0.4)

[B,bid_os] = ts3_get_rows(g);

% Parse optional args
opts.eratio = 1.03; opts.frac = 0.1; opts.hub_cell = 0.4; opts.straight = false;
for n = 1:2:numel(varargin)
    if strcmpi(varargin{n},'eratio'),   opts.eratio = varargin{n+1};
    elseif strcmpi(varargin{n},'frac'), opts.frac = varargin{n+1};
    elseif strcmpi(varargin{n},'hub_cell'), opts.hub_cell = varargin{n+1};
    elseif strcmpi(varargin{n},'straight'), opts.straight = true;
    end
end

% --- Find inlet patches on the first row ---
inbids = []; inpids = [];
for bid = B{1}
    if ~isfield(g{bid+1},'patch'), continue; end
    for pid = 0:numel(g{bid+1}.patch)-1
        p = g{bid+1}.patch{pid+1}.attribute;
        if p.kind == 0
            inbids = [inbids, bid];
            inpids = [inpids, pid];
        end
    end
end
if isempty(inbids), error('ts3_interface_block: no inlet patch found'); end

% --- Extract inlet face coordinates ---
xin = {}; rin = {}; tin = {}; xin1 = {}; rin1 = {}; tin1 = {};
for i = 1:numel(inbids)
    p = g{inbids(i)+1}.patch{inpids(i)+1}.attribute;
    I = p.ist+1:p.ien; J = p.jst+1:p.jen; K = p.kst+1:p.ken;
    xi = double(squeeze(g{inbids(i)+1}.x(I,J,K)));
    ri = double(squeeze(g{inbids(i)+1}.r(I,J,K)));
    rti = double(squeeze(g{inbids(i)+1}.rt(I,J,K)));
    if numel(K)==1 && K>1
        xi = xi'; ri = ri'; rti = rti';
        xi1 = double(squeeze(g{inbids(i)+1}.x(I,J,K-1)))';
        ri1 = double(squeeze(g{inbids(i)+1}.r(I,J,K-1)))';
        rti1 = double(squeeze(g{inbids(i)+1}.rt(I,J,K-1)))';
    elseif numel(I)==1 && I>1
        xi1 = double(squeeze(g{inbids(i)+1}.x(I-1,J,K)));
        ri1 = double(squeeze(g{inbids(i)+1}.r(I-1,J,K)));
        rti1 = double(squeeze(g{inbids(i)+1}.rt(I-1,J,K)));
    elseif numel(I)==1
        xi1 = double(squeeze(g{inbids(i)+1}.x(I+1,J,K)));
        ri1 = double(squeeze(g{inbids(i)+1}.r(I+1,J,K)));
        rti1 = double(squeeze(g{inbids(i)+1}.rt(I+1,J,K)));
    else
        xi1 = double(squeeze(g{inbids(i)+1}.x(I+1,J,K)));
        ri1 = double(squeeze(g{inbids(i)+1}.r(I+1,J,K)));
        rti1 = double(squeeze(g{inbids(i)+1}.rt(I+1,J,K)));
    end
    xin{i} = xi; rin{i} = ri; tin{i} = rti./ri;
    xin1{i} = xi1; rin1{i} = ri1; tin1{i} = rti1./ri1;
end

% Stitch inlet faces together by theta ordering
tmin = min(tin{1}(1));
order = 1:numel(inbids);
for i = 2:numel(inbids)
    if tin{i}(1) < tmin
        tmin = tin{i}(1);
        order = [i, order(order~=i)];
    end
end
xf = []; rf = []; tf = []; xf1 = []; rf1 = []; tf1 = [];
for i = order
    if isempty(xf)
        xf = xin{i}; rf = rin{i}; tf = tin{i};
        xf1 = xin1{i}; rf1 = rin1{i}; tf1 = tin1{i};
    else
        xf = [xf, xin{i}(:,2:end)];
        rf = [rf, rin{i}(:,2:end)];
        tf = [tf, tin{i}(:,2:end)];
        xf1 = [xf1, xin1{i}(:,2:end)];
        rf1 = [rf1, rin1{i}(:,2:end)];
        tf1 = [tf1, tin1{i}(:,2:end)];
    end
end
nj = size(rf,1); nk = size(rf,2);

% --- Grid dimensions ---
DX = mean(xf(:)) - xinlet;
dx0 = xf1 - xf;
dx = dx0; dx(:,2:end-1) = (dx0(:,1:end-2)+dx0(:,2:end-1)+dx0(:,3:end))/3;
dx(:,1) = dx(:,2); dx(:,end) = dx(:,end-1);
ds_mean = mean(dx(:));
% Two-stage hyperbolic bunching: fine cells near inlet, coarse upstream
if abs(opts.eratio-1) < 1e-3, er = 1+1e-3*sign(opts.eratio-1); else, er = opts.eratio; end
Sp1 = ds_mean/DX/opts.frac;
nia = max(round(1/Sp1^(2/(er+1))+1),3);
Sa = hyperbolic_bunch_1s(nia,Sp1);
ds_ena = (Sa(end)-Sa(end-1))*DX*opts.frac;
Sp2 = ds_ena/DX/(1-opts.frac);
nib = max(round(1/Sp2^(2/(er+1))+1),3);
Sb = hyperbolic_bunch_1s(nib,Sp2);
ni = nia + nib - 1;
S = [Sb(1:nib-1); Sa];

% --- Radial grid from endwall curves (transfinite interpolation) ---
if opts.straight
    xx = linspace(-0.5,xinlet,200)';
    xrcas = [xx, rf(end,1)*ones(size(xx))];
end
% Fallback endwalls when user provides empty arrays
if isempty(xrhub)
    rhub0 = rf(1,1);
    xrhub = [linspace(xinlet,xf(1,1),100)', linspace(rhub0,rf(1,1),100)'];
end
if isempty(xrcas)
    rcas0 = rf(end,1);
    xrcas = [linspace(xinlet,xf(end,1),100)', linspace(rcas0,rf(end,1),100)'];
end
xhub = linspace(xinlet,xf(1,1),ni)';
rhub_i = pchip(xrhub(:,1),xrhub(:,2),xhub);
rhub_i(rhub_i<0.0025) = 0.0025;
xcas = linspace(xinlet,xf(end,1),ni)';
rcas_i = pchip(xrcas(:,1),xrcas(:,2),xcas);
% Hub clustering
n_clust = ceil(nj/2);
rnorm = (rf(:,1)-rf(1,1))/(rf(end,1)-rf(1,1));
scale = (rcas_i(1)-rhub_i(1))/(rf(end,1)-rf(1,1));
ds_h = (rf(2,1)-rf(1,1))/(rf(end,1)-rf(1,1))/scale;
ds_c = (rf(end,1)-rf(end-1,1))/(rf(end,1)-rf(1,1))/scale;
rnorm2 = poly_dist_stretch(rnorm,ds_h,ds_c);
yplot = rnorm2*(rcas_i(1)-rhub_i(1)) + rhub_i(1);
rnorm2a = (yplot(1:n_clust)-yplot(1))/(yplot(n_clust)-yplot(1));
ds1 = opts.hub_cell*(rnorm2a(2)-rnorm2a(1));
ds2 = rnorm2a(end)-rnorm2a(end-1);
rnorm2b = hyperbolic_bunch_1s(n_clust,ds2/ds1)';
yplot(1:n_clust) = yplot(1) + rnorm2b*(yplot(n_clust)-yplot(1));
rnorm2a = (yplot(end-n_clust+1:end)-yplot(end-n_clust+1))/(yplot(end)-yplot(end-n_clust+1));
ds1 = rnorm2a(2)-rnorm2a(1);
ds2 = rnorm2a(end)-rnorm2a(end-1);
rnorm2b = hyperbolic_bunch_1s(n_clust,ds2/ds1)';
yplot(end-n_clust+1:end) = yplot(end-n_clust+1) + rnorm2b*(yplot(end)-yplot(end-n_clust+1));

% --- Build coordinates ---
xblock = zeros(ni,nj,nk,class(g{1}.x));
rblock = zeros(ni,nj,nk,class(g{1}.x));
tblock = zeros(ni,nj,nk,class(g{1}.x));
for j = 1:nj
    xot = reshape(xf(j,:),1,nk);
    xinlet_a = xinlet + opts.frac*(xot-xinlet);  % fine/coarse transition
    % i-index where fine section ends
    ia = find(S >= (xot-xinlet_a)/(xot-xinlet),1);
    if isempty(ia), ia = nia; end
    Sa_i = hyperbolic_bunch_1s(nia,dx(j,:)'/(xot(1)-xinlet_a(1)));
    Sb_i = hyperbolic_bunch_1s(nib,ds_ena/(DX*(1-opts.frac))*ones(1,nk));
    sst_a = ((xot-xinlet_a)./(xot-xinlet))';
    Sa_i = 1 - sst_a + sst_a .* reshape(Sa_i,nia,1,nk);
    Sb_i = reshape(Sb_i,nib,1,nk) .* (1 - sst_a);
    S_i = [Sb_i(1:nib-1,:,:); Sa_i];
    dxot = reshape(xot-xinlet,1,1,nk);
    xblock(:,j,:) = xinlet + S_i .* dxot;
    xblock(:,j,1) = 0.5*(xblock(:,j,1)+xblock(:,j,end));  % periodic
    xblock(:,j,end) = xblock(:,j,1);
    rblock(:,j,:) = interp1(yplot,xhub,reshape(xblock(:,j,:),ni*nk,1),'pchip');
    rblock = reshape(rblock,ni,nj,nk);
    t0 = repmat(reshape(linspace(tf(j,1),tf(j,end),nk),1,1,nk),ni,1,1);
    tblock(:,j,:) = t0 + S_i .* (repmat(reshape(tf(j,:),1,1,nk),ni,1,1) - t0);
end
rtblock = rblock .* tblock;

% --- Create the block ---
bid_new = numel(g);
blk = struct();
blk.av = g{inbids(1)+1}.av;
blk.bv = g{inbids(1)+1}.bv;
blk.patch = {};
blk.attribute = struct('ni',int32(ni),'nj',int32(nj),'nk',int32(nk),...
    'bid',int32(bid_new),'np',int32(0),'nc',int32(0),'ncl',int32(0));
blk.x = single(xblock); blk.r = single(rblock); blk.rt = single(rtblock);
g{bid_new+1} = blk;

% --- Add patches ---
bid = bid_new; g{bid+1}.patch = {};
% Inlet patch (kind 0)
p = struct('attribute',struct('kind',int32(0),'bid',int32(bid),'pid',int32(0),...
    'ist',int32(0),'ien',int32(1),'jst',int32(0),'jen',int32(nj),...
    'kst',int32(0),'ken',int32(nk)));
p.pstag = 101325*single(ones(nj*nk,1));
p.tstag = 288.15*single(ones(nj*nk,1));
p.yaw = single(zeros(nj*nk,1));
p.rfin = single(0.25); p.sfinlet = int32(0); p.pitch = single(zeros(nj*nk,1));
g{bid+1}.patch{1} = p;
g{bid+1}.attribute.np = int32(1);
% Connection to each old inlet
for i = 1:numel(inbids)
    % Old inlet becomes periodic connection to new block
    p0 = g{inbids(i)+1}.patch{inpids(i)+1};
    if isfield(p0,'pstag'), p0 = rmfield(p0,{'pstag','tstag','yaw','rfin','sfinlet','pitch'}); end
    p0.attribute.kind = int32(5);
    p0.attribute.nxbid = int32(bid); p0.attribute.nxpid = int32(i-1);
    p0.attribute.idir = int32(5); p0.attribute.jdir = int32(1); p0.attribute.kdir = int32(5);
    g{inbids(i)+1}.patch{inpids(i)+1} = p0;
    % New block's outlet connection
    p = struct('attribute',struct('kind',int32(5),'bid',int32(bid),'pid',int32(i),...
        'nxbid',int32(inbids(i)),'nxpid',int32(inpids(i)),...
        'ist',int32(ni-1),'ien',int32(ni),'jst',int32(0),'jen',int32(nj),...
        'kst',int32(0),'ken',int32(nk),'idir',int32(-1),'jdir',int32(1),'kdir',int32(2)));
    g{bid+1}.patch{end+1} = p;
    g{bid+1}.attribute.np = int32(g{bid+1}.attribute.np + 1);
end
% Periodic to self
p = struct('attribute',struct('kind',int32(5),'bid',int32(bid),'pid',int32(numel(g{bid+1}.patch)),...
    'nxbid',int32(bid),'nxpid',int32(numel(g{bid+1}.patch)+1),...
    'ist',int32(0),'ien',int32(ni),'jst',int32(0),'jen',int32(nj),...
    'kst',int32(0),'ken',int32(1),'idir',int32(0),'jdir',int32(1),'kdir',int32(-1)));
g{bid+1}.patch{end+1} = p;
p.attribute.pid = int32(numel(g{bid+1}.patch));
p.attribute.nxpid = int32(p.attribute.pid-1);
p.attribute.kst = int32(nk-1); p.attribute.ken = int32(nk);
g{bid+1}.patch{end+1} = p;
g{bid+1}.attribute.np = int32(g{bid+1}.attribute.np + 2);

% Match periodic patches of the new block
g = ts3_match_patches(g,bid);
end

%% ------------------------------------------------------------------------
function S = hyperbolic_bunch_1s(n,Sp1)
% One-sided hyperbolic bunching (n points, first spacing Sp1).
n = double(n(:)); Sp1 = Sp1(:);
if numel(n)==1 && numel(Sp1)>1, n = repmat(n,numel(Sp1),1); end
y = 1./Sp1./(n-1);
v = log(y);
w = 1./y - 0.028527431;
b = v + (1+1./v).*log(2*v) - 0.02041793 + 0.24902722*w + ...
    1.9496443*w.^2 - 2.6294547*w.^3 + 8.56795911*w.^4;
ybar = y(y<2.7829681)-1;
if any(y<2.7829681)
    b(y<2.7829681) = sqrt(6*ybar).*(1 - 0.15*ybar + 0.057321429*ybar.^2 - ...
        0.024907295*ybar.^3 + 0.0077424461*ybar.^4 - 0.0010794123*ybar.^5);
end
b(abs(b.^2) < 1e-7) = 1e-7;
R = repmat(linspace(0,1,max(n)),numel(b),1);
if numel(n)>1
    for i = 1:numel(n)
        if n(i) < size(R,2), R(i,n(i)+1:end) = NaN; end
    end
end
S = 1 + tanh(b.*(R-1))./tanh(b);
S = reshape(S,size(n,1),[]);
end

function Sout = poly_dist_stretch(Sin,s1,s2)
% Cubic polynomial offset to adjust first/last grid spacings.
s10 = Sin(2);
s20 = Sin(end-1);
if isempty(s1), s1 = s10; end
if ~exist('s2','var') || isempty(s2), s2 = 1-s20; end
d1 = s1-s10;
d2 = (1-s2)-s20;
d1a = s1 + s1*(Sin(3)-s10)/s10 - Sin(3);
d2a = 1-s2 - (s20-Sin(end-2))/(1-s20)*s2 - Sin(end-2);
A = [1,s10,s10^2,s10^3; 1,s20,s20^2,s20^3; ...
     1,Sin(3),Sin(3)^2,Sin(3)^3; 1,Sin(end-2),Sin(end-2)^2,Sin(end-2)^3];
p = flipud(A\[d1;d2;d1a;d2a]);
Sout = Sin;
Sout(2:end-1) = Sout(2:end-1) + polyval(p,Sin(2:end-1));
end