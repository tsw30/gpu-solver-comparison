function g = ts3_load_grid(bcs_file, p3d_file)
%ts3_load_grid Read AutoGrid plot3d + BCS files into Turbostream cell array.
%   Condensed replacement for ts_autogrid_reader (original in ref/).

%% Read plot3d coordinate file
fprintf('ts3_load_grid: Reading %s\n',p3d_file);
raw = cell2mat(textscan(fopen(p3d_file,'r'),'%f'));
nblk = raw(1);
g = cell(nblk,1);
st = 5 + (nblk-1)*3;
for n = 1:nblk
    bid = n-1;
    g{n}.attribute.bid = bid;
end
for n = 1:nblk
    fprintf('    Reading block %d\n',n-1);
    ni = raw(2+(n-1)*3);
    nj = raw(3+(n-1)*3);
    nk = raw(4+(n-1)*3);
    nt = nk*nj*ni;
    x = zeros(nk,nj,ni,'single');
    y = zeros(nk,nj,ni,'single');
    z = zeros(nk,nj,ni,'single');
    for k=1:nk
        for j=1:nj
            for i=1:ni
                x(k,j,i) = single(raw(st));
                y(k,j,i) = single(raw(st+nt));
                z(k,j,i) = single(raw(st+2*nt));
                st = st+1;
            end
        end
    end
    st = st + 2*nt;
    b = struct();
    b.av = struct();
    b.bv = struct();
    b.r = y; b.rt = z; b.x = x;
    b.attribute = struct('np',int32(0),'nc',int32(0),'ncl',int32(0),...
        'bid',int32(n-1),'procid',int32(0),'threadid',int32(0),...
        'ni',int32(ni),'nj',int32(nj),'nk',int32(nk));
    g{n} = b;
end

%% Read BCS patch file
fprintf('ts3_load_grid: Reading %s\n',bcs_file);
mixing_plane_map = cell(nblk,1);
patch_face_map = cell(nblk,1);
for n=1:nblk
    mixing_plane_map{n} = {};
    for iface=1:6
        patch_face_map{n}{iface} = {};
    end
end
fid = fopen(bcs_file,'r');
l = fgetl(fid);
for n=1:nblk
    fprintf('    Reading block %d\n',n-1);
    [g,fid,mixing_plane_map,patch_face_map] = ...
        read_patches(g,n-1,fid,mixing_plane_map,patch_face_map);
end
fclose(fid);
g = connect_patches(g,mixing_plane_map,patch_face_map);

%% Permute coordinates: plot3d (k,j,i) -> TS (i,j,k), flip i, (y,z) -> (r,theta)
for n=1:nblk
    bid = n-1;
    b = g{n};
    x = g{n}.x; y = g{n}.r; z = g{n}.rt;
    ni = b.attribute.ni; nj = b.attribute.nj; nk = b.attribute.nk;
    ni2 = b.attribute.nk; nj2 = b.attribute.nj; nk2 = b.attribute.ni;
    x2 = zeros(nk2,nj2,ni2,'single');
    y2 = zeros(nk2,nj2,ni2,'single');
    z2 = zeros(nk2,nj2,ni2,'single');
    for j=1:nj
        for k=1:nk
            for i=1:ni
                i2 = k; j2 = j; k2 = i;
                x2(k2,j2,i2) = x(k,j,i);
                y2(k2,j2,i2) = y(k,j,i);
                z2(k2,j2,i2) = z(k,j,i);
            end
        end
    end
    b.attribute.ni = int32(ni2);
    b.attribute.nj = int32(nj2);
    b.attribute.nk = int32(nk2);
    x3 = zeros(nk2,nj2,ni2,'single');
    y3 = zeros(nk2,nj2,ni2,'single');
    z3 = zeros(nk2,nj2,ni2,'single');
    for j=1:nj2
        for k=1:nk2
            for i=1:ni2
                i2 = ni2-i+1; j2 = j; k2 = k;
                x3(k2,j2,i2) = x2(k,j,i);
                y3(k2,j2,i2) = y2(k,j,i);
                z3(k2,j2,i2) = z2(k,j,i);
            end
        end
    end
    r = sqrt(y3.*y3 + x3.*x3);
    t = atan2(y3,x3);
    for j=1:nj2
        for k=1:nk2
            for i=1:ni2-1
                dt = t(k,j,i+1)-t(k,j,i);
                if abs(dt)>pi
                    isign = -1;
                    if dt < 1
                        isign = 1;
                    end
                    for ii = i+1:ni2
                        t(k,j,ii) = t(k,j,ii) + 2*pi*isign;
                    end
                end
            end
        end
    end
    for i=1:ni2
        for k=1:nk2
            for j=1:nj2-1
                dt = t(k,j+1,i)-t(k,j,i);
                if abs(dt)>pi
                    isign = -1;
                    if dt < 1
                        isign = 1;
                    end
                    for jj = j+1:nj2
                        t(k,jj,i) = t(k,jj,i) + 2*pi*isign;
                    end
                end
            end
        end
    end
    for j=1:nj2
        for i=1:ni2
            for k=1:nk2-1
                dt = t(k+1,j,i)-t(k,j,i);
                if abs(dt)>pi
                    isign = -1;
                    if dt < 1
                        isign = 1;
                    end
                    for kk = k+1:nk2
                        t(kk,j,i) = t(kk,j,i) + 2*pi*isign;
                    end
                end
            end
        end
    end
    rt = r.*t;
    z4 = permute(z3,[3,2,1]);
    r  = permute(r,[3,2,1]);
    rt = permute(rt,[3,2,1]);
    g{n} = b;
    g{n}.x = z4;
    g{n}.r = r;
    g{n}.rt = rt;
end

%% Swap i/k in patch attributes (AutoGrid reads k as axial, TS uses i as axial)
for i=1:nblk
    bid = i-1;
    for pid = ts_get_pids(bid,g)
        patch = ts_get_patch(bid,pid,g);
        ist = patch.attribute.ist;
        ien = patch.attribute.ien;
        kst = patch.attribute.kst;
        ken = patch.attribute.ken;
        patch.attribute.ist = kst;
        patch.attribute.ien = ken;
        patch.attribute.kst = ist;
        patch.attribute.ken = ien;
        if patch.attribute.kind==5 || patch.attribute.kind==15
            idir = patch.attribute.idir;
            kdir = patch.attribute.kdir;
            dirmap = int32([2,1,0,5,4,3,6]);
            patch.attribute.idir = kdir;
            patch.attribute.kdir = idir;
            patch.attribute.idir = dirmap(patch.attribute.idir+1);
            patch.attribute.jdir = dirmap(patch.attribute.jdir+1);
            patch.attribute.kdir = dirmap(patch.attribute.kdir+1);
        end
        g = add_patch(g,patch);
    end
end

%% Flip i-direction in patches (TS axial direction opposite to AutoGrid)
for i=1:nblk
    bid = i-1;
    b = ts_get_block(bid,g);
    for pid = ts_get_pids(bid,g)
        patch = ts_get_patch(bid,pid,g);
        ist = patch.attribute.ist;
        ien = patch.attribute.ien;
        patch.attribute.ist = b.attribute.ni - ien;
        patch.attribute.ien = b.attribute.ni - ist;
        if patch.attribute.kind==5 || patch.attribute.kind==15
            idirmap = int32([0,4,5,3,1,2,6]);
            jkdirmap = int32([3,1,2,0,4,5,6]);
            patch.attribute.idir = idirmap(patch.attribute.idir+1);
            patch.attribute.jdir = jkdirmap(patch.attribute.jdir+1);
            patch.attribute.kdir = jkdirmap(patch.attribute.kdir+1);
        end
        g = add_patch(g,patch);
    end
end

%% Clean up extra patch fields
for i=1:nblk
    bid = i-1;
    for pid = ts_get_pids(bid,g)
        g = clean_patch(g,bid,pid);
    end
end

%% Set defaults
g = ts_set_default_avs(g);
for i=1:nblk
    bid = i-1;
    nblade = g{i}.bv.nblade;
    g = ts_set_default_bvs(g,bid);
    g{i}.bv.nblade = int32(nblade);
    g{i}.bv.fblade = single(nblade);
end

end

%% ------------------------------------------------------------------------
%  Local helpers matching the original ref/ functions exactly

function [self,fid,mixing_plane_map,patch_face_map] = ...
    read_patches(g,bid,fid,mixing_plane_map,patch_face_map)
self = g;
patch_kind.PER = 5;
patch_kind.CON = 5;
patch_kind.INL = 0;
patch_kind.OUT = 1;
patch_kind.NMB = 15;
patch_kind.PERNB = 15;
patch_kind.ROT = 2;
patch_kind.EXT = 13;
patch_kind.MIR = 7;
iface_map = [0,4,3,1,5,2];
b = ts_get_block(bid,self);
l = fgetl(fid);
l = fgetl(fid);
l = strsplit(fgetl(fid));
if strcmpi(l{1},'NONE')
    nblade = 1;
else
    nblade = l{end-2};
    nblade = str2double(nblade);
end
self{bid+1}.bv.nblade = int32(nblade);
self{bid+1}.bv.fblade = single(nblade);
for iface = iface_map
    l = strsplit(fgetl(fid));
    np = str2double(l{1});
    for ipatch = 1:np
        patch = struct();
        l = strsplit(fgetl(fid));
        label = l{1};
        type = l{2};
        tcomp = sum(strcmp(type,{'INL','OUT','EXT','MIR','PER','CON','ROT','NMB','PERNB'}));
        if tcomp==1
            st0 = str2double(l{3});
            st1 = str2double(l{4});
            en0 = str2double(l{5});
            en1 = str2double(l{6});
        end
        tsum1 = sum(strcmp(type,{'PER','CON','NMB','PERNB'}));
        if tsum1==1
            nxbid = str2double(l{7})-1;
            nxiface = iface_map(str2double(l{8}));
            nxipatch = str2double(l{9})-1;
            ori = str2double(l{10});
            cor = str2double(l{11});
        end
        if tcomp==1
            patch.attribute.bid = bid;
            patch.attribute.kind = patch_kind.(type);
            if iface==0||iface==2
                if iface==0
                    patch.attribute.ist = 0;
                    patch.attribute.ien = 1;
                end
                if iface==2
                    patch.attribute.ist = b.attribute.ni-1;
                    patch.attribute.ien = b.attribute.ni;
                end
                patch.attribute.jst = st0-1;
                patch.attribute.jen = en0;
                patch.attribute.kst = st1-1;
                patch.attribute.ken = en1;
            end
            if iface==4||iface==5
                if iface==4
                    patch.attribute.jst = 0;
                    patch.attribute.jen = 1;
                end
                if iface==5
                    patch.attribute.jst = b.attribute.nj-1;
                    patch.attribute.jen = b.attribute.nj;
                end
                patch.attribute.ist = st0-1;
                patch.attribute.ien = en0;
                patch.attribute.kst = st1-1;
                patch.attribute.ken = en1;
            end
            if iface==1||iface==3
                if iface==3
                    patch.attribute.kst = 0;
                    patch.attribute.ken = 1;
                end
                if iface==1
                    patch.attribute.kst = b.attribute.nk-1;
                    patch.attribute.ken = b.attribute.nk;
                end
                patch.attribute.ist = st0-1;
                patch.attribute.ien = en0;
                patch.attribute.jst = st1-1;
                patch.attribute.jen = en1;
            end
            [self,patch.attribute.pid] = add_patch(self,patch);
        end
        if tsum1==1
            patch.attribute.nxbid = nxbid;
            patch.attribute.iface = iface;
            patch.attribute.ipatch = ipatch-1;
            patch.attribute.nxiface = nxiface;
            patch.attribute.nxipatch = nxipatch;
            patch.attribute.ori = ori;
            patch.attribute.cor = cor;
            self = add_patch(self,patch);
            patch_face_map{bid+1}{iface+1}{ipatch} = patch;
        end
        if strcmp(type,'ROT')
            patch.attribute.iface = iface;
            patch.attribute.ipatch = ipatch-1;
            labsplit = strsplit(label,'(');
            labsplit = strsplit(labsplit{2},')');
            patch.attribute.label = str2double(labsplit{1});
            [self,patch.attribute.pid] = add_patch(self,patch);
            mixing_plane_map{bid+1}{iface+1} = patch;
            patch_face_map{bid+1}{iface+1}{ipatch} = patch;
        end
    end
end
end

function self = connect_patches(g,mixing_plane_map,patch_face_map)
self = g;
for bl=1:length(self)
    bid = bl-1;
    for pid = ts_get_pids(bid,self)
        patch = ts_get_patch(bid,pid,self);
        if patch.attribute.kind==2
            for nb = 1:length(mixing_plane_map)
                for nf = 1:length(mixing_plane_map{nb})
                    nxp = mixing_plane_map{nb}{nf};
                    if isfield(nxp,'attribute')
                        if isfield(nxp.attribute,'label')
                            psum = sum([patch.attribute.bid,patch.attribute.iface]==[nxp.attribute.bid,nxp.attribute.iface]);
                            if patch.attribute.label == nxp.attribute.label && psum<2
                                nxp_keep = nxp;
                                break
                            end
                        end
                    end
                end
            end
            patch.attribute.nxbid = nxp_keep.attribute.bid;
            patch.attribute.nxpid = nxp_keep.attribute.pid;
            self = add_patch(self,patch);
        end
        if patch.attribute.kind == 5 || patch.attribute.kind == 15
            nxp = patch_face_map{patch.attribute.nxbid+1}{patch.attribute.nxiface+1}{patch.attribute.nxipatch+1};
            patch.attribute.nxpid = nxp.attribute.pid;
            if patch.attribute.iface==0||patch.attribute.iface==2
                dir0 = 1;
                dir1 = 2;
            end
            if patch.attribute.iface==4||patch.attribute.iface==5
                dir0 = 0;
                dir1 = 2;
            end
            if patch.attribute.iface==1||patch.attribute.iface==3
                dir0 = 0;
                dir1 = 1;
            end
            if patch.attribute.nxiface==0||patch.attribute.nxiface==2
                nxdir0 = 1;
                nxdir1 = 2;
            end
            if patch.attribute.nxiface==4||patch.attribute.nxiface==5
                nxdir0 = 0;
                nxdir1 = 2;
            end
            if patch.attribute.nxiface==1||patch.attribute.nxiface==3
                nxdir0 = 0;
                nxdir1 = 1;
            end
            if ismember([patch.attribute.ori,patch.attribute.cor],[0,1;0,3;1,0;1,2],'rows')
                temp = nxdir1;
                nxdir1 = nxdir0;
                nxdir0 = temp;
            end
            if ismember([patch.attribute.ori,patch.attribute.cor],[0,2;0,3;1,1;1,2],'rows')
                nxdir0 = nxdir0+3;
            end
            if ismember([patch.attribute.ori,patch.attribute.cor],[0,1;0,2;1,2;1,3],'rows')
                nxdir1 = nxdir1+3;
            end
            nxdir = [6,6,6];
            nxdir(dir0+1) = nxdir0;
            nxdir(dir1+1) = nxdir1;
            patch.attribute.idir = nxdir(1);
            patch.attribute.jdir = nxdir(2);
            patch.attribute.kdir = nxdir(3);
            self = add_patch(self,patch);
        end
    end
end
end

function self = clean_patch(self,bid,pid)
bln = bid+1;
patch = self{bln}.patch{pid+1}.attribute;
fnames = fieldnames(patch);
dnames = {'label','ipatch','iface','nxipatch','nxiface','ori','cor'};
comp = setdiff(dnames,fnames);
cnames = setdiff(dnames,comp);
for i=1:length(cnames)
    name = cnames{i};
    self{bln}.patch{pid+1}.attribute = rmfield(self{bln}.patch{pid+1}.attribute,name);
end
end

function [self,pid] = add_patch(self,patch)
bid = patch.attribute.bid;
bln = bid+1;
if ~isfield(self{bln},'patch')
    np = 0;
elseif isempty(self{bln}.patch)
    np = 0;
else
    np = length(self{bln}.patch);
end
if isfield(patch.attribute,'pid')
    pid = patch.attribute.pid;
else
    pid = np;
    patch.attribute.pid = np;
    self{bln}.attribute.np = int32(self{bln}.attribute.np + 1);
end
innames = {'bid','pid','kind','ist','ien','jst','jen','kst','ken','idir',...
    'jdir','kdir','nxbid','nxpid','nface','nt','label','ipatch','iface',...
    'nxipatch','nxiface','ori','cor'};
for i=1:length(innames)
    name = innames{i};
    if ~isfield(patch.attribute,name)
        patch.attribute.(name) = 0;
        if i==16, patch.attribute.(name) = 1; end
    end
    self{bln}.patch{pid+1,1}.attribute.(name) = int32(patch.attribute.(name));
end
end

%% ------------------------------------------------------------------------
%  Inline minified versions of ts_get_pids, ts_get_patch, ts_get_block,
%  ts_set_default_avs, ts_set_default_bvs

function pids = ts_get_pids(varargin)
if length(varargin) == 1
    Data = varargin{1};
    b = Data{1};
else
    bid = varargin{1};
    Data = varargin{2};
    b = Data{bid+1};
end
if isfield(b,'patch')
    pids = zeros(1,length(b.patch));
    for n=1:length(b.patch)
        pids(n) = b.patch{n}.attribute.pid;
    end
else
    pids = [];
end
end

function p = ts_get_patch(bid,pid,Data)
b = ts_get_block(bid,Data);
p = b.patch{pid+1};
if p.attribute.pid ~= pid
    for n=1:length(b.patch)
        if b.patch{n}.attribute.pid == pid
            p = b.patch{n}; break;
        end
    end
end
end

function b = ts_get_block(bid,Data)
b = Data{bid+1};
if b.attribute.bid ~= bid
    for n=1:length(Data)
        if Data{n}.attribute.bid == bid
            b = Data{n}; break;
        end
    end
end
end

function self = ts_set_default_avs(self)
av = {'adaptive_smoothing',int32(1);'cfl',single(0.4);'cp',single(1005);...
      'dampin',single(10);'dts',int32(0);'dts_conv',single(0);...
      'fac_sa_smth',single(4);'fac_sa_step',single(1);'fac_st0',single(1.5);...
      'fac_st1',single(1);'fac_st2',single(1);'fac_st3',single(1);...
      'fac_stmix',single(0);'fac_wall',single(1);'facsafe',single(0);...
      'facsecin',single(0.8);'frequency',single(0);'ga',single(1.4);...
      'if_no_mg',int32(0);'ifgas',int32(0);'ilos',int32(1);...
      'nchange',int32(1000);'ncycle',int32(0);'nlos',int32(5);...
      'nomatch_int',int32(1);'nstep',int32(10000);'nstep_cycle',int32(0);...
      'nstep_inner',int32(0);'nstep_save',int32(0);'nstep_save_probe',int32(0);...
      'nstep_save_start',int32(100000);'nstep_save_start_probe',int32(0);...
      'poisson_cfl',single(0.7);'poisson_nsmooth',int32(0);...
      'poisson_nstep',int32(1000);'poisson_restart',int32(0);...
      'poisson_sfin',single(0.05);'prandtl',single(1);'restart',int32(0);...
      'rfmix',single(0);'rfvis',single(0.2);'rg_cp0',single(1005);...
      'rg_cp1',single(0);'rg_cp2',single(0);'rg_cp3',single(0);...
      'rg_cp4',single(0);'rg_cp5',single(0);'rg_rgas',single(287.15);...
      'sfin',single(0.01);'turb_vis_damp',single(1);'turbvis_lim',single(3000);...
      'viscosity',single(1.8e-5);'viscosity_law',int32(0);'wall_law',int32(0);...
      'write_force',int32(0);'write_yplus',int32(0);'write_tdamp',int32(0);...
      'write_egen',int32(0);'ifsuperfac',int32(0);'poisson_limit',int32(0);...
      'tref',single(300);'pref',single(1e5);'nspecies',int32(1);...
      'sf_scalar',single(0.05);'schmidt_0',single(1);'schmidt_1',single(1);...
      'cp0_0',single(1005);'cp1_0',single(0);'cp2_0',single(0);...
      'cp3_0',single(0);'cp4_0',single(0);'cp5_0',single(0);...
      'rgas_0',single(287.15);'cp0_1',single(1005);'cp1_1',single(0);...
      'cp2_1',single(0);'cp3_1',single(0);'cp4_1',single(0);...
      'cp5_1',single(0);'rgas_1',single(287.15);'cfl_ko',single(0.4);...
      'cfl_st_ko',single(0.01);'cfl_en_ko',single(0.4);'sfin_sa',single(0.05);...
      'sfin_ko',single(0.05);'ko_dist',single(1e-4);'ko_restart',int32(0);...
      'fac_st0_option',int32(0);'sa_helicity_option',int32(0);...
      'sa_ch1',single(1);'sa_ch2',single(1);'precon',int32(0);...
      'smooth_scale_dts_option',int32(1);'smooth_scale_precon_option',int32(0);...
      'smooth_scale_directional_option',int32(0)};
for i=1:length(self)
    for r=1:size(av,1)
        self{i}.av.(av{r,1}) = av{r,2};
    end
end
end

function self = ts_set_default_bvs(self,bid)
bv = {'rpm',single(0);'rpmi1',single(0);'rpmi2',single(0);...
      'rpmj1',single(0);'rpmj2',single(0);'rpmk1',single(0);'rpmk2',single(0);...
      'fmgrid',single(0.4);'poisson_fmgrid',single(0.1);'xllim',single(0);...
      'vgridin',single(50);'pstatin',single(8e5);'tstagin',single(1200);...
      'vgridout',single(50);'pstatout',single(8e5);'tstagout',single(1200);...
      'ftype',int32(0);'nblade',int32(1);'fblade',single(1);...
      'fracann',single(1);'sfin_mul',single(1);'facsecin_mul',single(1);...
      'dampin_mul',single(1);'fsturb',single(1);'nimixl',int32(0);...
      'srough_i0',single(0);'srough_i1',single(0);'srough_j0',single(0);...
      'srough_j1',single(0);'srough_k0',single(0);'srough_k1',single(0);...
      'itrans',int32(0);'itrans_j1_st',int32(0);'itrans_j2_st',int32(0);...
      'itrans_k1_st',int32(0);'itrans_k2_st',int32(0);...
      'itrans_j1_en',int32(0);'itrans_j2_en',int32(0);...
      'itrans_k1_en',int32(0);'itrans_k2_en',int32(0);...
      'itrans_j1_frac',single(0);'itrans_j2_frac',single(0);...
      'itrans_k1_frac',single(0);'itrans_k2_frac',single(0);...
      'jtrans',int32(0);'jtrans_i1_st',int32(0);'jtrans_i2_st',int32(0);...
      'jtrans_k1_st',int32(0);'jtrans_k2_st',int32(0);...
      'jtrans_i1_en',int32(0);'jtrans_i2_en',int32(0);...
      'jtrans_k1_en',int32(0);'jtrans_k2_en',int32(0);...
      'jtrans_i1_frac',single(0);'jtrans_i2_frac',single(0);...
      'jtrans_k1_frac',single(0);'jtrans_k2_frac',single(0);...
      'ktrans',int32(0);'ktrans_i1_st',int32(0);'ktrans_i2_st',int32(0);...
      'ktrans_j1_st',int32(0);'ktrans_j2_st',int32(0);...
      'ktrans_i1_en',int32(0);'ktrans_i2_en',int32(0);...
      'ktrans_j1_en',int32(0);'ktrans_j2_en',int32(0);...
      'jtrans_i1_frac',single(0);'jtrans_i2_frac',single(0);...
      'ktrans_i1_frac',single(0);'ktrans_i2_frac',single(0);...
      'superfac',single(0);'fac_st0',single(1);'ndup_phaselag',int32(1);...
      'turb_intensity',single(5)};
blk = bid+1;
for r=1:size(bv,1)
    self{blk}.bv.(bv{r,1}) = bv{r,2};
end
end