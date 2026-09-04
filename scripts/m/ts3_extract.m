function [results,dbg] = ts3_extract(g,return_dbg)
%ts3_extract Extract clean flow performance at rig measurement planes for
%a uniform inflow (single passage) CFD solution.
%
%   results = ts3_extract(g)
%
%   g is the solver block cell array (0-based block ids at g{bid+1}).
%   Optional [results,dbg] = ts3_extract(g,true) also returns dbg.(row)
%   structs holding the intermediate cut/plane/mix for inspection.
%
%   For each of the 3 rig rakes (L00S, ESIN, BPIN) defined in
%   scripts/measurement_coords.json:
%     1. cut the row blocks along the rake line (triangulated plane),
%     2. mass average To and integrate mass flow on the raw cut,
%     3. interpolate the cut onto a structured span x pitch grid,
%     4. mix out circumferentially for spanwise Po/To/Vx/Vr/Vt profiles,
%     5. sample the profiles at rake span stations and average.
%
%   Self-contained condensation of the legacy extract_measurements + ts_*
%   dependency chain (originals kept in ref/):
%     - rows from periodic patch connectivity, ordered by o-mesh leading
%       edge x (replaces ts_solution_info/ts_get_row_bids/ts_get_bids)
%     - hub/casing from kind-5 periodic patch surfaces (replaces ts_get_xr
%       and its ts_get_periodics/ts_get_patch_cut/ts_structured_cut chain)
%     - angled isosurface cut restricted to the row's own blocks, so rows
%       with different pitch/periodicity stay separate (replaces
%       ts_unstructured_cut; no plotting)
%     - span stations from the o-mesh radial grid (replaces ts_tinterp)
%     - theta grid spans one pitch from the low-theta periodic face as
%       traced by the cut, so skewed (non-constant-theta) periodic faces
%       are handled
%     - mix-out assumes forward axial flow: no reference frame velocity
%       shift (replaces ts_mix_out_average_la_new)
%     - mass flows are full annulus: per-passage triangle fluxes scaled
%       by bv.fblade (as the ref cut_averages)
%     - secondary variables limited to those needed (replaces ts_secondary)

rows = {'L00S','ESIN','BPIN'};
nk = 80; % pitchwise resolution of the interpolated cut plane

% measurement rake (x,r) coordinates for L00S, ESIN, BPIN
coords_file = fullfile(fileparts(mfilename('fullpath')),'../measurement_coords.json');
if ~exist(coords_file,'file')
    error('ts3_extract: measurement coordinates file not found: %s', coords_file);
end
coords = jsondecode(fileread(coords_file));

[B,bid_os] = get_rows(g);
if numel(B) < 3
    error('ts3_extract: expected at least 3 blade rows, found %d', numel(B));
end
[PAMB,TAMB] = get_ambient(g);
if nargin < 2, return_dbg = false; end
dbg = struct();

results = struct();
for n = 1:3
    nm  = rows{n};
    bl  = g{bid_os(n)+1};              % o-mesh block of row n
    crd = coords.(nm);
    [xrh, xrc] = hub_casing(g,B{n});   % hub and casing streamwise lines
    ang = 180/pi*atan2(crd(end,2)-crd(1,2), crd(end,1)-crd(1,1));

    cut = cut_plane(g,find_cut_blocks(g,crd,ang,B,n),bid_os(n),crd(1,:),ang);
    xref = [crd(1,1),crd(end,1)];
    plane = interp_plane(cut,span_stations(bl),xref,crd(1,:),ang,xrh,xrc,nk);
    mix = mix_out(plane,xrh,xrc);
    rn  = spanwise_values(mix,crd,n);
    [rn.To_mm,rn.mdot] = mass_avg_to(cut);
    rn.x = crd(:,1);
    rn.r = crd(:,2);
    rn.area = pi*(rn.r(end)^2 - rn.r(1)^2);
    rn.area_weights = area_weights(rn.r);
    rn.MO = mix;
    results.(nm) = rn;
    if return_dbg
        dbg.(nm) = struct('cut',cut,'plane',plane,'mix',mix);
    end
end
rpm = double(g{bid_os(1)+1}.bv.rpm);
rpm_corr = abs(rpm)/sqrt(TAMB/288.15);

NRT = rpm_corr/(7224.315/1.05);

results = evaluate_averages(results,PAMB,TAMB);
results.NRT = NRT;
end

%% Row identification

function [B,bid_os] = get_rows(g)
% Blade row block groups from periodic (kind 5 or 16) patch connectivity,
% following the same logic as ref ts_get_row_bids but condensed.
%
% An "o-mesh block" (bid_os) is a blade-passage block with a self-periodic
% wrap-around i-face (kind 5, nxbid==bid, ist==ni-1, ien==ni) AND a wall
% at k=0 (kwall==1).  Blocks without a full-j k=0 face (e.g. inlet plena)
% have kwall=0 and are rejected.
%
% Row groups are formed by flood-filling through periodic connections,
% then sorted by o-mesh leading-edge x.
nb = numel(g);
% --- kwall and periodic connectivity ---
kwall = ones(nb,1);
periodic = false(nb);
for bid = 0:nb-1
    bl = g{bid+1};
    if ~isfield(bl,'patch'), continue; end
    for pid = 0:numel(bl.patch)-1
        p = bl.patch{pid+1}.attribute;
        kst = p.kst; ken = p.ken;
        jst = p.jst; jen = p.jen;
        nj  = bl.attribute.nj;
        % kwall: a single-k-plane (kst==ken-1) full-j-span face at
        % kst==0 that is not a wall or sliding-plane makes kwall=0
        if kst==ken-1 && nj==jen-jst && p.kind~=8 && p.kind~=14
            if kst == 0
                kwall(bid+1) = 0;
            end
        end
        % periodic connectivity via kind 5 or 16
        if p.kind==5 || p.kind==16
            periodic(bid+1,p.nxbid+1) = true;
        end
    end
end
periodic = periodic | periodic.';
% --- o-mesh detection (kind 5, self, i-wrap, kwall==1) ---
os = false(nb,1);
for bid = 0:nb-1
    bl = g{bid+1};
    if ~isfield(bl,'patch') || kwall(bid+1)==0, continue; end
    for pid = 0:numel(bl.patch)-1
        p = bl.patch{pid+1}.attribute;
        if p.kind==5 && p.nxbid==bid && ...
                p.ist==bl.attribute.ni-1 && p.ien==bl.attribute.ni
            os(bid+1) = true;
            break
        end
    end
end
bid_os = find(os)-1;
if isempty(bid_os)
    B = {}; return
end
% --- flood-fill row groups ---
B = {};
used = false(nb,1);
for b0 = find(any(periodic,2))'                     % start from blocks with connections
    if used(b0), continue; end
    grp = b0; used(b0) = true;
    i = 1;
    while i <= numel(grp)
        nbrs = find(periodic(grp(i),:));
        nbrs = nbrs(~used(nbrs));
        used(nbrs) = true;
        grp = [grp; nbrs(:)]; %#ok<AGROW>
        i = i+1;
    end
    B{end+1} = grp-1; %#ok<AGROW>                  % 0-based
end
% --- assign each o-mesh block to its row, error if not one-per-row ---
nrows = numel(B);
bid_os_u = bid_os(:);
bid_os_in_row = cell(nrows,1);
for n = 1:nrows
    mask = ismember(bid_os_u,B{n});
    bid_os_in_row{n} = bid_os_u(mask);
end
os_per_row = cellfun(@numel,bid_os_in_row);
if any(os_per_row==0)
    error('ts3_extract: row has no o-mesh block');
end
bid_os = cellfun(@(x) x(1),bid_os_in_row);         % one per row
% --- order rows axially by o-mesh leading edge ---
xos = zeros(nrows,1);
for n = 1:nrows
    xos(n) = min(g{bid_os(n)+1}.x(:));
end
[~,ix] = sort(xos);
B = B(ix); bid_os = bid_os(ix);
end

function bids = find_cut_blocks(g,crd,ang,B,row)
% Blocks intersected by the cut plane, preferring blocks from the
% designated row B{row} but falling back to any intersecting block when
% the row's blocks do not contain the plane.
nrow = numel(B);
bids_all = find_intersecting(g,crd,ang);
% select blocks from the designated row if any intersect
bids = intersect(bids_all, B{row}(:));
if ~isempty(bids), return; end
% fallback: find the row whose mean radius best matches the rake
rm = @(bl) mean(bl.r(:));
rrake = mean(crd(:,2));
best = row;
dmin = inf;
for n = 1:nrow
    for k = 1:numel(B{n})
        d = abs(rm(g{B{n}(k)+1}) - rrake);
        if d < dmin, best = n; dmin = d; end
    end
end
bids = intersect(bids_all, B{best}(:));
if isempty(bids)
    bids = bids_all;                                  % desperate fallback
end
end

function bids = find_intersecting(g,crd,ang)
% All blocks whose volume the cut plane intersects (xrot spanning zero
% and radial overlap with the rake).
x0 = crd(1,1); r0 = crd(1,2);
cs = pi/2 - ang*pi/180;
rlo = min(crd(:,2)); rhi = max(crd(:,2));
bids = [];
for bid = 0:numel(g)-1
    bl = g{bid+1};
    if ~isfield(bl,'x') || ~isfield(bl,'r'), continue; end
    rr = double(bl.r);
    if max(rr(:)) < rlo || min(rr(:)) > rhi, continue; end
    xx = double(bl.x);
    xrot = (xx-x0)*cos(cs) - (rr-r0)*sin(cs);
    if min(xrot(:)) <= 0 && max(xrot(:)) >= 0
        bids = [bids, bid]; %#ok<AGROW>
    end
end
end

function [PAMB,TAMB] = get_ambient(g)
% Ambient conditions: max stagnation values over the inlet (kind 0) patch
for bid = 0:numel(g)-1
    bl = g{bid+1};
    if ~isfield(bl,'patch'), continue; end
    for pid = 0:numel(bl.patch)-1
        p = bl.patch{pid+1};
        if p.attribute.kind == 0
            PAMB = double(max(p.pstag(:)));
            TAMB = double(max(p.tstag(:)));
            return
        end
    end
end
error('ts3_extract: no inlet patch found');
end

%% Geometry helpers

function [xrh,xrc] = hub_casing(g,B_row)
% Hub and casing (x,r) polylines of a row from the kind-5 periodic patch
% surfaces of its blocks (ts_get_xr logic): j=1 of a patch surface is hub,
% j=nj is casing.  The o-mesh i-wrap self-periodic patch is skipped (it is
% a passage face, not an endwall).  Points are sorted by x and uniqued so
% pchip interpolation along each polyline is well posed.  Only the block
% coordinate arrays are read - no flow fields are copied.
xhub = []; rhub = []; xcas = []; rcas = [];
for bid = B_row(:)'
    bl = g{bid+1};
    if ~isfield(bl,'patch'), continue; end
    for pid = 0:numel(bl.patch)-1
        p = bl.patch{pid+1}.attribute;
        if p.kind ~= 5, continue; end
        if p.nxbid==bid && p.ist==bl.attribute.ni-1 && p.ien==bl.attribute.ni
            continue                     % o-mesh i-wrap self-periodic
        end
        % patch surface (x,r) window, singleton dimensions squeezed
        ist = p.ist+1; ien = p.ien;
        jst = p.jst+1; jen = p.jen;
        kst = p.kst+1; ken = p.ken;
        xw = double(bl.x(ist:ien,jst:jen,kst:ken));
        rw = double(bl.r(ist:ien,jst:jen,kst:ken));
        sw = size(xw); sw = sw(sw>1);
        xw = reshape(xw,sw); rw = reshape(rw,sw);
        if p.ist+1 == p.ien                % constant-i patch -> (j,k) surface
            xhub = [xhub; xw(1,:)'];  rhub = [rhub; rw(1,:)'];   % j=1 hub
            xcas = [xcas; xw(end,:)']; rcas = [rcas; rw(end,:)']; % j=nj casing
        else                               % constant-k patch -> (i,j) surface
            xhub = [xhub; xw(:,1)];   rhub = [rhub; rw(:,1)];    % j=1 hub
            xcas = [xcas; xw(:,end)];  rcas = [rcas; rw(:,end)];  % j=nj casing
        end
    end
end
[xhub,is] = sort(xhub); rhub = rhub(is);
[xhub,iu] = unique(xhub); rhub = rhub(iu);
[xcas,is] = sort(xcas); rcas = rcas(is);
[xcas,iu] = unique(xcas); rcas = rcas(iu);
xrh = [xhub(:), rhub(:)];
xrc = [xcas(:), rcas(:)];
end

function rnd = span_stations(bl)
% Normalised radial stations (0=hub, 1=casing) of the o-mesh j-grid.
% The o-mesh block's own j=1..nj grid defines hub and casing at its
% axial location, so span is exact without pchip extrapolation risk.
rj = double(squeeze(bl.r(1,:,1)));
rnd = (rj(:)-rj(1))/(rj(end)-rj(1));
end

%% Cut plane extraction

function [g,bids] = split_o(g,bids,bid_o)
% Split the o-mesh block at mid-i into two blocks (ref 'split_o') so its
% wrapped i-direction cannot produce spurious cut triangles.
bids = bids(:)';                                    % ensure row vector
bl = g{bid_o+1};
ni = bl.attribute.ni;
imid = (ni+1)/2;
bl1 = bl;
fns = fieldnames(bl);
for f = 1:numel(fns)
    a = bl.(fns{f});
    if isfloat(a) && isequal(size(a),[ni,bl.attribute.nj,bl.attribute.nk])
        bl.(fns{f})  = a(1:imid,:,:);
        bl1.(fns{f}) = a(imid:end,:,:);
    end
end
bl.attribute.ni  = imid;
bl1.attribute.ni = ni-imid+1;
bl1.attribute.bid = numel(g);
g{bid_o+1} = bl;
g{end+1} = bl1;
bids = [bids,bl1.attribute.bid];
end

function cut = cut_plane(g,bids,bid_o,xr0,ang)
% Triangulated cut through (x,r) point xr0 at angle ang [deg] from the
% x axis, across the row blocks bids (o-mesh split at mid-i first).
[g,bids] = split_o(g,bids,bid_o);
cs = pi/2 - ang*pi/180;              % rotation convention of the ref cut
vars = {'ro','rovx','rovr','rorvt','roe','x'};
cut = struct('faces',[],'y',[],'z',[]);
for v = 1:numel(vars), cut.(vars{v}) = []; end
for bid = bids
    bl = g{bid+1};
    % reshape all arrays to 3D if necessary (isosurface requires 3D)
    sz = [bl.attribute.ni, bl.attribute.nj, bl.attribute.nk];
    xx = double(bl.x); rr = double(bl.r); rrt = double(bl.rt);
    if ndims(xx) < 3, xx = reshape(xx,sz); rr = reshape(rr,sz); rrt = reshape(rrt,sz); end
    xrot = (xx-xr0(1))*cos(cs) - (rr-xr0(2))*sin(cs);
    if min(xrot(:))>=0 || max(xrot(:))<=0
        continue                     % cut plane misses this block
    end
    [yy,zz] = pol2cart(rrt./rr,rr);
    % ref convention: isosurface(y,z,var, xrot, 0, var) - X=y, Y=z, Z=var
    % (3D grid), V=xrot (iso-field), isovalue=0, data=var.  The faces and
    % the (y,z) vertex positions depend only on xrot (marching cubes), so
    % the same call with any other variable in place of var returns the
    % identical geometry with that variable edge-interpolated in the 3rd
    % column.  That is the ref's own way of sampling every field on the
    % cut, and far cheaper than a scatteredInterpolant over the whole
    % block per variable (which dominated run time on real grids).
    ro3 = double(bl.ro); if ndims(ro3) < 3, ro3 = reshape(ro3,sz); end
    [fc,vtx] = isosurface(yy,zz,ro3,xrot,0,ro3);
    if isempty(fc), continue; end
    n0 = numel(cut.ro);
    cut.faces = [cut.faces; fc+n0];
    cut.y = [cut.y; vtx(:,1)];
    cut.z = [cut.z; vtx(:,2)];
    cut.ro = [cut.ro; vtx(:,3)];
    for v = 2:numel(vars)
        fv = double(bl.(vars{v}));
        if ndims(fv) < 3, fv = reshape(fv,sz); end
        [~,vtx_v] = isosurface(yy,zz,fv,xrot,0,fv);
        if ~isequal(size(vtx_v),size(vtx))
            error('ts3_extract: isosurface geometry changed between variables');
        end
        cut.(vars{v}) = [cut.(vars{v}); vtx_v(:,3)];
    end
    if ~isfield(cut,'av')
        cut.av = bl.av; cut.bv = bl.bv;
    end
end
if isempty(cut.faces)
    error('ts3_extract: cut plane does not intersect the row blocks');
end
[tn,rn] = cart2pol(cut.y,cut.z);
cut.r = rn; cut.rt = tn.*rn;
end

%% Structured interpolation and mix-out

function plane = interp_plane(cut,rnd,xref,xr0,ang,xrh,xrc,nk)
% Re-interpolate the triangulated cut onto a structured (nj x nk) plane
% spanning one blade pitch: span stations from rnd; theta at each station
% starts on the low-theta periodic face (found by pitch-periodic matching
% of cut boundary vertices) and spans exactly one pitch, so k=1 and k=nk
% land on the periodic face pair even when that face is not a constant-
% theta surface (purpose of the ref ts_tinterp interpolant).
cs = pi/2 - ang*pi/180;
rhub = pchip(xrh(:,1),xrh(:,2),xref(1));
rcas = pchip(xrc(:,1),xrc(:,2),xref(2));
rj = rhub + rnd(:)*(rcas-rhub);      % span stations on the cut plane
xj = xr0(1) + (rj-xr0(2))*tan(cs);   % axial location on the cut line
pitch = 2*pi/double(cut.bv.fblade);
[rf,tf] = periodic_face_theta(cut,pitch);    % low-theta periodic face
if isempty(rf)
    thf_v = min(cut.rt./cut.r);      % fallback: lowest cut theta
    thf   = repmat(thf_v,numel(rj),1);
else
    thf_v = interp1(rf,tf,cut.r,'linear','extrap');
    thf   = interp1(rf,tf,rj,'linear','extrap');
end
% interpolate in (r, theta relative to the face) coordinates: the cut maps
% to a near-rectangular domain there, so hub/casing wall and periodic face
% target lines stay inside the interpolation hull (in (y,z), wall arcs sag
% outside the Delaunay hull and force nearest extrapolation)
tl = cut.rt./cut.r - thf_v;
[Rq,Tq] = ndgrid(rj,(0:nk-1)/(nk-1)*pitch);
vars = {'ro','rovx','rovr','rorvt','roe'};
for v = 1:numel(vars)
    % double casts: solver/bv fields may be single precision on real grids
    F = scatteredInterpolant(double(cut.r),double(tl),double(cut.(vars{v})),'linear','nearest');
    plane.(vars{v}) = F(double(Rq),double(Tq));
end
th = thf + Tq;                       % one blade pitch per row
plane.r  = Rq;
plane.x  = repmat(xj,1,nk);
plane.rt = plane.r.*th;
plane.av = cut.av;
plane.bv = cut.bv;
end

function [rf,tf] = periodic_face_theta(cut,pitch)
% (r,theta) polyline of the cut's low-theta periodic face, found by
% matching pitch-periodic boundary vertex pairs: a boundary vertex whose
% periodic image (same x,r at theta+pitch) is also a boundary vertex lies
% on the low-theta face - the passage spans exactly one pitch, so no
% interior or wall vertex has such a partner.  Returns empty when fewer
% than two matches are found.  Candidate partners lie in a +-0.1*pitch
% theta window, so the search walks a theta-sorted sliding window instead
% of scanning all pairs, which dominated run time on real-grid cuts.
f = cut.faces;
e = sort([f(:,[1 2]); f(:,[2 3]); f(:,[3 1])],2);
[eu,~,ic] = unique(e,'rows');
be = eu(accumarray(ic,1)==1,:);                 % boundary edges
vb = unique(be(:));
nb = numel(vb);
thlo = nan(nb,1);
if nb < 2, rf = []; tf = []; return; end
rb = cut.r(vb); tb = cut.rt(vb)./rb; xb = cut.x(vb);
rrng = max(cut.r)-min(cut.r)+1e-12;
xrng = max(cut.x)-min(cut.x)+1e-12;
% window [tb+0.9*pitch, tb+1.1*pitch] per vertex, found by binary search
% on the distinct sorted thetas (interp1 needs unique sample points)
[~,ord] = sort(tb);
tbs = tb(ord);
[tbsu,fst] = unique(tbs);                       % first occurrence indices
lst = [fst(2:end)-1; nb];                       % last occurrence indices
lo_u = interp1(tbsu,1:numel(tbsu),tb+0.9*pitch,'next');      % first >= lo
hi_u = interp1(tbsu,1:numel(tbsu),tb+1.1*pitch,'previous'); % last <= hi
win = ~isnan(lo_u) & ~isnan(hi_u);
lo = nan(nb,1); hi = nan(nb,1);
lo(win) = fst(lo_u(win));
hi(win) = lst(hi_u(win));
% walk the (short) candidate windows, vectorised over all vertices
dmin = inf(nb,1);
jbest = zeros(nb,1);
if any(win)
    maxw = max(hi(win)-lo(win)) + 1;
    for k = 0:maxw-1
        j = lo + k;
        act = win & j <= hi;
        jj = ord(j(act));
        ii = find(act);
        d = abs(tb(jj)-tb(ii)-pitch);
        good = d < 0.1*pitch & ...
            abs(rb(jj)-rb(ii)) < 1e-6*rrng & ...
            abs(xb(jj)-xb(ii)) < 1e-6*xrng & d < dmin(ii);
        ii = ii(good);
        dmin(ii) = d(good);
        jbest(ii) = jj(good);
    end
end
hit = jbest > 0;
if sum(hit) < 2, rf = []; tf = []; return; end
thlo(hit) = tb(jbest(hit)) - pitch;
sel = hit;
rf = rb(sel); tf = thlo(sel);
[rf,ix] = sort(rf);
tf = tf(ix);
[rf,iu] = unique(rf,'stable');
tf = tf(iu);
end

function mix = mix_out(plane,xrh,xrc)
% Circumferential mix-out average of a structured (span x pitch) cut plane
% per spanwise station, conserving mass, momentum and energy (condensed
% ref ts_mix_out_average_la_new, dimension 2 only; assumes forward axial
% flow, no reference frame velocity shift). Secondary variables returned.
cp = double(plane.av.cp); ga = double(plane.av.ga);
Rgas = cp*(ga-1)/ga; cv = Rgas/(ga-1);
% nodal velocities, static T/P and stagnation enthalpy
Vx = plane.rovx./plane.ro;
Vr = plane.rovr./plane.ro;
Vt = plane.rorvt./(plane.ro.*plane.r);
T  = (plane.roe./plane.ro - 0.5*(Vx.^2+Vr.^2+Vt.^2))/cv;
P  = plane.ro.*T.*Rgas;
ho = (plane.roe + P)./plane.ro;
t  = plane.rt./plane.r;
% cell-centred (4-node) averages
cc = @(f) 0.25*(f(1:end-1,1:end-1)+f(2:end,1:end-1)+f(1:end-1,2:end)+f(2:end,2:end));
rc  = cc(plane.r); tc = cc(t); roc = cc(plane.ro);
Vxc = cc(Vx); Vrc = cc(Vr); Vtc = cc(Vt); Pc = cc(P); hoc = cc(ho);
% cell area components from cross product of diagonals
x1  = plane.x(2:end,2:end)-plane.x(1:end-1,1:end-1);
x2  = plane.x(1:end-1,2:end)-plane.x(2:end,1:end-1);
r1  = plane.r(2:end,2:end)-plane.r(1:end-1,1:end-1);
r2  = plane.r(1:end-1,2:end)-plane.r(2:end,1:end-1);
rt1 = plane.r(2:end,2:end).*(t(2:end,2:end)-tc) - plane.r(1:end-1,1:end-1).*(t(1:end-1,1:end-1)-tc);
rt2 = plane.r(1:end-1,2:end).*(t(1:end-1,2:end)-tc) - plane.r(2:end,1:end-1).*(t(2:end,1:end-1)-tc);
Ax = 0.5*(r1.*rt2 - r2.*rt1);
Ar = 0.5*(x2.*rt1 - x1.*rt2);
neg = Ax<0; Ax(neg) = -Ax(neg); Ar(neg) = -Ar(neg);  % flux in +x sense
An = hypot(Ax,Ar);
% pitchwise (dim 2) flux sums -> conservation integrals per span station
dVA  = Vxc.*Ax + Vrc.*Ar;
Anj  = sum(An,2);
AxAn = sum(Ax,2)./Anj;
ArAn = sum(Ar,2)./Anj;
I1 = sum(roc.*dVA,2)./Anj;
I2 = sum(roc.*Vxc.*dVA + Pc.*Ax,2)./Anj;
I3 = sum(roc.*Vtc.*dVA.*rc,2)./Anj;
I4 = sum(roc.*Vrc.*dVA + Pc.*Ar,2)./Anj;
I5 = sum(roc.*hoc.*dVA,2)./Anj;
mix_r = sqrt(sum(roc.*dVA.*rc.^2,2)./I1./Anj);
% mixed-out static pressure from the conservation integrals
Am = AxAn.*I2 + ArAn.*I4;
Pm = real((Am + sqrt(Am.^2 + (ga^2-1)*(I2.^2 + (I3./mix_r).^2 + I4.^2 - 2*I5.*I1)))/(ga+1));
Vxm = (I2 - Pm.*AxAn)./I1;
Vrm = (I4 - Pm.*ArAn)./I1;
Vtm = I3./(mix_r.*I1);
rom = I1./(Vxm.*AxAn + Vrm.*ArAn);
% re-average cell values onto nodal span stations
nf = @(v) [v(1); 0.5*(v(1:end-1)+v(2:end)); v(end)];
mix.r = mean(plane.r,2);             % nodal radii -> span covers 0..1
mix.x = mean(plane.x,2);
rom = nf(rom); Vxm = nf(Vxm); Vrm = nf(Vrm); Vtm = nf(Vtm); Pm = nf(Pm);
% mixed-out primary and secondary variables
Tm = Pm./(rom*Rgas);
V2 = Vxm.^2 + Vrm.^2 + Vtm.^2;
FM = 1 + (ga-1)/2*V2./(ga*Rgas*Tm);
mix.ro = rom;
mix.rovx = rom.*Vxm; mix.rovr = rom.*Vrm; mix.rorvt = rom.*Vtm.*mix.r;
mix.roe = rom.*(cv*Tm + 0.5*V2);
mix.P = Pm; mix.T = Tm;
mix.Po = Pm.*FM.^(ga/(ga-1));
mix.To = Tm.*FM;
mix.Vx = Vxm; mix.Vr = Vrm; mix.Vt = Vtm;
% span fraction from hub_casing curves at each spanwise station's axial x
rhub = pchip(xrh(:,1),xrh(:,2),mix.x);
rcas = pchip(xrc(:,1),xrc(:,2),mix.x);
mix.span = (mix.r - rhub)./(rcas - rhub);
end

%% Averages and performance parameters

function [To_mm,mdot] = mass_avg_to(cut)
% Mass-averaged stagnation temperature and mass flow of a triangulated
% cut, scaled to full annulus (multiply per-passage flux by fblade, as
% ref cut_averages in ts_unstructured_cut).
cp = double(cut.av.cp); ga = double(cut.av.ga);
Rgas = cp*(ga-1)/ga; cv = Rgas/(ga-1);
fblade = double(cut.bv.fblade);
f1 = cut.faces(:,1); f2 = cut.faces(:,2); f3 = cut.faces(:,3);
% triangle areas from edge lengths (Heron, longest-side stable form)
a = hypot(cut.y(f2)-cut.y(f1), cut.z(f2)-cut.z(f1));
b = hypot(cut.y(f3)-cut.y(f2), cut.z(f3)-cut.z(f2));
c = hypot(cut.y(f1)-cut.y(f3), cut.z(f1)-cut.z(f3));
d = sort([a,b,c],2,'descend');
dA = real(0.25*sqrt(4*d(:,1).^2.*d(:,2).^2 - (d(:,1).^2+d(:,2).^2-d(:,3).^2).^2));
rovx_f = (cut.rovx(f1)+cut.rovx(f2)+cut.rovx(f3))/3;
dM = dA.*rovx_f * fblade;           % full-annulus mass flux
mdot = sum(dM);
% stagnation temperature at vertices, face averaged
Vx = cut.rovx./cut.ro; Vr = cut.rovr./cut.ro; Vt = cut.rorvt./(cut.ro.*cut.r);
V2 = Vx.^2 + Vr.^2 + Vt.^2;
T  = (cut.roe./cut.ro - 0.5*V2)/cv;
To = T.*(1 + (ga-1)/2*V2./(ga*Rgas*T));
To_f = (To(f1)+To(f2)+To(f3))/3;
To_mm = sum(To_f.*dM)/mdot;
end

function out = spanwise_values(mix,crd,row)
% Sample mixed-out profiles at rake radial stations.  pchip is evaluated
% only within the data range (no extrapolation); hub/casing end values
% are then linearly extrapolated from the adjacent interior points.
out = struct();
rmin = mix.r(1); rmax = mix.r(end);
for nm = {'Po','To','Vx','Vr','Vt'}
    % clamp only the query points that fall outside the data range
    rq = crd(:,2);
    rq(rq < rmin) = rmin;
    rq(rq > rmax) = rmax;
    w = pchip(mix.r,mix.(nm{1}),rq);
    % Hub: copy from second point (row 1) or linear extrapolation
    if row == 1
        w(1) = w(2);
    else
        w(1) = w(2) + (w(2)-w(3))/(rq(2)-rq(3)+eps)*(crd(1,2)-rq(2));
    end
    % Casing: linear extrapolation
    w(end) = w(end-1) + (w(end-1)-w(end-2))/(rq(end-1)-rq(end-2)+eps)*(crd(end,2)-rq(end-1));
    out.(nm{1}) = w;
end
end

function w = area_weights(r)
% Normalised annular-area weights of radial samples
w = zeros(numel(r),1);
w(1) = 0.25*(r(1)+r(2))^2 - r(1)^2;
w(end) = r(end)^2 - 0.25*(r(end)+r(end-1))^2;
w(2:end-1) = 0.25*(r(3:end)+r(2:end-1)).^2 - 0.25*(r(2:end-1)+r(1:end-2)).^2;
w = w/sum(w);
end

function results = evaluate_averages(results,PAMB,TAMB)
% Area/mass weighted plane averages and derived performance parameters
for nm = {'L00S','ESIN','BPIN'}
    p = results.(nm{1});
    results.(nm{1}).Po_avg = sum(p.Po.*p.area_weights);
    results.(nm{1}).To_avg = sum(p.To.*p.area_weights);
end
results.M1 = results.L00S.mdot;
results.M2 = results.ESIN.mdot;
results.MF = results.M1 - results.M2;
results.PIN = results.L00S.Po_avg;
results.TIN = results.L00S.To_avg;
results.PAMB = PAMB;
results.TAMB = TAMB;
results.PF = (results.ESIN.Po_avg*results.ESIN.area + results.BPIN.Po_avg*results.BPIN.area)/(results.ESIN.area + results.BPIN.area);
results.TF = (results.ESIN.To_avg*results.ESIN.area + results.BPIN.To_avg*results.BPIN.area)/(results.ESIN.area + results.BPIN.area);
results.TFSH = (results.ESIN.To_mm*results.ESIN.mdot + results.BPIN.To_mm*results.BPIN.mdot)/(results.ESIN.mdot + results.BPIN.mdot);
results.TRRTSH = results.TFSH/results.L00S.To_mm - 1;
results.TRRT = results.TF/results.TIN - 1;
results.RPTR = results.PF/results.PIN;
end
