function [xrh,xrc] = ts3_hub_casing(g,B_row)
%ts3_hub_casing Hub and casing (x,r) polylines of a blade row.
%   [xrh,xrc] = ts3_hub_casing(g,B_row)
%
%   Condensed replacement for the legacy ts_get_xr and its
%   ts_get_periodics/ts_get_patch_cut/ts_structured_cut chain (originals
%   kept in ref/).  Shared with ts3_extract.
%
%   Hub/casing are read from the kind-5 periodic patch surfaces of the
%   row's blocks: j=1 of a patch surface is hub, j=nj is casing.  The
%   o-mesh i-wrap self-periodic patch is skipped (it is a passage face,
%   not an endwall).  Points are sorted by x and uniqued so pchip
%   interpolation along each polyline is well posed.  Only the block
%   coordinate arrays are read - no flow fields are copied.
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