function g = ts3_smooth_interstage(g)
%ts3_smooth_interstage Redefine tangential mesh distribution between rotor
%TE and mixing plane to remove AutoGrid twisting.
%   Condensed replacement for ifan_smooth_interstage (original in ref/).

[B,bid_os] = ts3_get_rows(g);
row1 = B{1};
o1 = bid_os(1);  % first-row o-mesh block

% --- Find core, bypass, splitter, and inlet blocks ---
% Core: the o-mesh block among row1 whose periodic neighbours don't span
% full nj.  Bypass: the other o-mesh block.  Inlet: the block connected
% to both core and bypass via periodic patches but not an o-mesh block.
nj_o = g{o1+1}.attribute.nj;
core = o1; bypass = o1; splitter = []; inlet = [];
for bid = row1
    if bid == o1, continue; end
    if g{bid+1}.attribute.nj < nj_o
        % Check if it's a splitter (o-mesh-like) or inlet
        if g{bid+1}.attribute.nj == nj_o
            splitter = bid;
        end
        continue;
    end
end
% Find core and bypass by checking periodic connections
for bid = row1
    if ~isfield(g{bid+1},'patch'), continue; end
    for pid = 0:numel(g{bid+1}.patch)-1
        p = g{bid+1}.patch{pid+1}.attribute;
        if p.kind == 5 && p.nxbid ~= bid
            % This block is connected to another via periodic — likely core or bypass
            if g{bid+1}.attribute.nj >= nj_o
                if bid == o1, core = bid; else, bypass = bid; end
            end
            % The other end is the inlet block
            if ~ismember(p.nxbid, row1)
                inlet = [inlet, p.nxbid];
            end
        end
    end
end
inlet = unique(inlet);
if isempty(inlet), error('ts3_smooth_interstage: inlet block not found'); end
bid_in = inlet(1);

% Get hub/casing profiles
[xrh, xrc] = ts3_hub_casing(g, B{1});

% --- Radial target distribution from side surfaces ---
% Quadratic fit of spanwise theta on core and bypass faces
cface = g{core+1}.rt(:,1,1)./g{core+1}.r(:,1,1);
bface = g{bypass+1}.rt(:,1,1)./g{bypass+1}.r(:,1,1);
rs = [g{core+1}.r(:,1,1); g{bypass+1}.r(:,1,1)];
ts = [g{core+1}.rt(:,1,1)./g{core+1}.r(:,1,1); ...
      g{bypass+1}.rt(:,1,1)./g{bypass+1}.r(:,1,1)];
side_fit = polyfit(double(rs),double(ts),2);

% --- Tangential smoothing with sine transition (x-direction) ---
x_in = double(g{bid_in+1}.x(1,:,1));
x_out = double(g{bid_in+1}.x(end,:,1));
xblock = double(g{bid_in+1}.x);
rblock = double(g{bid_in+1}.r);
ni = size(xblock,1); nj = size(xblock,2); nk = size(xblock,3);
xnorm = (xblock - repmat(reshape(x_in,1,nj,1),ni,1,1)) ...
      ./ (repmat(reshape(x_out,1,nj,1),ni,1,nk) - repmat(reshape(x_in,1,nj,1),ni,1,1));
tblock = double(g{bid_in+1}.rt)./rblock;

% Target theta at outlet from side fit
t_out = repmat(reshape(polyval(side_fit,rblock(end,:,1)),1,nj,1),ni,1,nk);
Dt_out = t_out - repmat(reshape(tblock(end,:,1),1,nj,1),ni,1,nk);
dt = (sin(pi*(xnorm-0.5))+1)/2 .* Dt_out;
tblock = tblock + dt;
g{bid_in+1}.rt = single(tblock .* rblock);

% --- Pitchwise smoothing (k-direction uniform spacing) ---
t = double(g{bid_in+1}.rt)./double(g{bid_in+1}.r);
tst = repmat(t(:,:,1),1,1,nk);
ten = repmat(t(:,:,end),1,1,nk);
tnorm = repmat(reshape(linspace(0,1,nk),1,1,nk),ni,nj,1);
dtnorm = tnorm - (t - tst)./(ten - tst + eps);
dtnorm = (sin(pi*(xnorm-0.5))+1)/2 .* dtnorm;
tblock = tst + (tnorm + dtnorm).*(ten - tst);
g{bid_in+1}.rt = single(tblock .* double(g{bid_in+1}.r));
end