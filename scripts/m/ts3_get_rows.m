function [B,bid_os] = ts3_get_rows(g)
%ts3_get_rows Blade row block groups from periodic patch connectivity.
%   [B,bid_os] = ts3_get_rows(g)
%
%   Condensed replacement for the legacy ts_get_row_bids + ts_get_bids
%   chain (originals kept in ref/).  Shared with ts3_extract.
%
%   An "o-mesh block" (bid_os) is a blade-passage block with a self-
%   periodic wrap-around i-face (kind 5, nxbid==bid, ist==ni-1, ien==ni)
%   AND a wall at k=0 (kwall==1).  Blocks without a full-j k=0 face (e.g.
%   inlet plena) have kwall=0 and are rejected.
%
%   Row groups are formed by flood-filling through periodic connections,
%   then sorted by o-mesh leading-edge x.

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
    error('ts3_get_rows: row has no o-mesh block');
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