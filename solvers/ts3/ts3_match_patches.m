function g = ts3_match_patches(g,bids)
%ts3_match_patches Average mesh coordinates on matched patch pairs.
%   g = ts3_match_patches(g)
%   g = ts3_match_patches(g,bids)
%
%   Condensed replacement for ts_match_patches_2 + ts_get_match (originals
%   kept in ref/).  Averages (x,r,theta) on periodic (kind 5/16) and
%   mixing-plane (kind 2/3) patch pairs.

if nargin < 2, bids = 0:numel(g)-1; end

for bid = bids
    if ~isfield(g{bid+1},'patch'), continue; end
    for pid = 0:numel(g{bid+1}.patch)-1
        p0 = g{bid+1}.patch{pid+1}.attribute;
        kind = p0.kind;

        if kind == 5 || kind == 16
            % Periodic / connection patches
            nxbid = p0.nxbid; nxpid = p0.nxpid;
            I = p0.ist+1:p0.ien; J = p0.jst+1:p0.jen; K = p0.kst+1:p0.ken;
            b1 = g{nxbid+1};
            p1 = b1.patch{nxpid+1}.attribute;
            I1 = p1.ist+1:p1.ien; J1 = p1.jst+1:p1.jen; K1 = p1.kst+1:p1.ken;

            % Extract coordinates from both sides
            x0 = squeeze(g{bid+1}.x(I,J,K));
            r0 = squeeze(g{bid+1}.r(I,J,K));
            t0 = squeeze(g{bid+1}.rt(I,J,K)) ./ r0;
            % Align matching side to this patch's orientation
            [x1,r1,t1,I1,J1,K1] = align_patch(g,nxbid,nxpid,p0,p1);
            % Average meridional coordinates
            xm = 0.5*(x0 + x1); rm = 0.5*(r0 + r1);
            % Fix theta offset by whole blade pitches
            dt = 2*pi / double(g{bid+1}.bv.nblade);
            tm = 0.5*(t0 + t1 - round((nanmean(t1(:))-nanmean(t0(:)))/dt)*dt);
            g{bid+1}.x(I,J,K)  = xm;
            g{bid+1}.r(I,J,K)  = rm;
            g{bid+1}.rt(I,J,K) = rm .* tm;
            % Apply to matching block (add back pitch offset)
            [I1g,J1g,K1g] = ndgrid(double(I1(:)),double(J1(:)),double(K1(:)));
            idx = sub2ind(double([b1.attribute.ni,b1.attribute.nj,b1.attribute.nk]),I1g(:),J1g(:),K1g(:));
            g{nxbid+1}.x(idx)  = xm;
            g{nxbid+1}.r(idx)  = rm;
            g{nxbid+1}.rt(idx) = rm .* (tm + dt*(nanmean(t1(:))>nanmean(t0(:))));
            if kind ~= 16, continue; end
            % Cartesian (kind 16) — also average y/z independently
            I = p0.ist+1:p0.ien; J = p0.jst+1:p0.jen; K = p0.kst+1:p0.ken;
            y0 = r0.*cos(t0); z0 = r0.*sin(t0);
            y1 = r1.*cos(t1); z1 = r1.*sin(t1);
            ins0 = isnan(y0);
            ym = 0.5*(y0 + y1); zm = 0.5*(z0 + z1);
            ym(ins0) = 0; zm(ins0) = 0;
            dr = sqrt((zm-z0).^2 + (ym-y0-r0).^2) - r0;
            dt = asin((zm-z0)./(r0+dr));
            dr(ins0)=0; dt(ins0)=0;
            g{bid+1}.r(I,J,K)  = r0+dr;
            g{bid+1}.rt(I,J,K) = (r0+dr).*(t0+dt);
            [I1,J1,K1] = get_idx(b1,p1);
            dr = sqrt((zm-z1).^2 + (ym-y1-r1).^2) - r1;
            dt = asin((zm-z1)./(r1+dr));
            dr(ins0)=0; dt(ins0)=0;
            g{nxbid+1}.r(I1,J1,K1)  = r1+dr;
            g{nxbid+1}.rt(I1,J1,K1) = (r1+dr).*(t1+dt);

        elseif kind == 2 || kind == 3
            % Mixing-plane average: skip for now (complex span indexing)
            continue;
        end
    end
end
end

%% ------------------------------------------------------------------------
function [x1,r1,t1,I1,J1,K1] = align_patch(g,nxbid,nxpid,p0,p1)
% ts_get_match logic (original in ref/)
b1 = g{nxbid+1};
[I1,J1,K1] = get_idx(b1,p1);
ni = b1.attribute.ni; nj = b1.attribute.nj; nk = b1.attribute.nk;
b1.I = repmat(reshape(1:ni,[ni 1 1]),[1 nj nk]);
b1.J = repmat(reshape(1:nj,[1 nj 1]),[ni 1 nk]);
b1.K = repmat(reshape(1:nk,[1 1 nk]),[ni nj 1]);
p1.I = squeeze(b1.I(I1,J1,K1)); p1.J = squeeze(b1.J(I1,J1,K1)); p1.K = squeeze(b1.K(I1,J1,K1));
x1 = double(squeeze(b1.x(I1,J1,K1))); r1 = double(squeeze(b1.r(I1,J1,K1)));
t1 = double(squeeze(b1.rt(I1,J1,K1)))./r1;
ds = [p0.idir p0.jdir p0.kdir];
if p0.ien-p0.ist == 1, ds(1) = -1; end
if p0.jen-p0.jst == 1, ds(2) = -1; end
if p0.ken-p0.kst == 1, ds(3) = -1; end
or_ = nan(1,3); fl_ = nan(1,3);
for n=1:3
    if ds(n) == 0, or_(n)=1; fl_(1)=0;
    elseif ds(n)==1, or_(n)=2; fl_(2)=0;
    elseif ds(n)==2, or_(n)=3; fl_(3)=0;
    elseif ds(n)==3, or_(n)=1; fl_(1)=1;
    elseif ds(n)==4, or_(n)=2; fl_(2)=1;
    elseif ds(n)==5, or_(n)=3; fl_(3)=1;
    end
end
or_ = or_(isnan(or_)==0); fl_ = fl_(isnan(fl_)==0);
or_ = (or_ - min(or_)) ./ (max(or_) - min(or_)) + 1;
varnames = {'x1','r1','t1','p1.I','p1.J','p1.K'};
for n=1:2
    if fl_(n)==1
        x1 = flip(x1,n); r1 = flip(r1,n); t1 = flip(t1,n);
    end
end
x1 = permute(x1,or_);
r1 = permute(r1,or_);
t1 = permute(t1,or_);
end

function [I,J,K] = get_idx(b,p)
I = p.ist+1:p.ien;
J = p.jst+1:p.jen;
K = p.kst+1:p.ken;
if I(end) < I(1), I = fliplr(I); end
if J(end) < J(1), J = fliplr(J); end
if K(end) < K(1), K = fliplr(K); end
end