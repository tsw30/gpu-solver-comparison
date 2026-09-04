%test_ts3_extract Verify ts3_extract against the dummy_grid test case.
%
%   dummy_grid builds a 9 block, 3 row single passage grid with identical
%   uniform axial flow in every row, so ts3_extract must return:
%     - spanwise profiles of Po/To/Vx/Vr/Vt equal to the uniform state
%       derived from the primary flow variables,
%     - mass flows equal to ro*Vx*(annulus area)/nblade per row,
%     - RPTR = 1, TRRT = 0, TRRTSH = 0 and PAMB/TAMB equal to the inlet
%       stagnation conditions.
%   A second part skews row 1's periodic faces by one full pitch and
%   verifies interp_plane still tracks them (non-constant theta faces).
%
%   Run from anywhere:  run('tests/test_ts3_extract.m')

clearvars; clc;
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root_dir,'scripts','m'));

%% build grid and extract
g = dummy_grid(false);
res = ts3_extract(g);

%% expected uniform state
% dummy_grid builds roe on cv*T, so the grid stagnation state is exactly
% the inlet stagnation conditions (p_ref, T_ref) with Vx = 10, Vr = Vt = 0.
cp = 1005; ga = 1.4;
Rgas = cp*(ga-1)/ga;
T_ref = 288.15; p_ref = 101325; Vx_exp = 10;
To_exp = T_ref; Po_exp = p_ref;
% density (mirrors dummy_grid) for the analytic mass flow check
vrc = Vx_exp/sqrt(cp*T_ref);
M2 = vrc^2/(ga-1)/(1-vrc^2/2);
ro = (p_ref/(Rgas*T_ref))/(1+(ga-1)/2*M2)^(1/(ga-1));

coords = jsondecode(fileread(fullfile(root_dir,'scripts','measurement_coords.json')));
rows = {'L00S','ESIN','BPIN'};
NB = 20;                                         % dummy grid blade count

%% checks
check(isfield(res,'L00S') && isfield(res,'ESIN') && isfield(res,'BPIN'), ...
    'row planes present');
for nm = {'M1','M2','MF','PIN','TIN','PAMB','TAMB','PF','TF','TFSH', ...
        'TRRTSH','TRRT','RPTR'}
    check(isfield(res,nm{1}), ['field ' nm{1} ' present']);
end

for n = 1:3
    nm = rows{n};
    p = res.(nm);
    crd = coords.(nm);
    check(numel(p.Po)==size(crd,1) && numel(p.To)==size(crd,1), ...
        [nm ': profile lengths match rake']);
    check(abs(p.MO.span(1))<1e-10 && abs(p.MO.span(end)-1)<1e-10, ...
        [nm ': mix-out span covers 0..1']);
    check(max(abs(p.Po-Po_exp))/Po_exp < 1e-6, [nm ': Po profile uniform']);
    check(max(abs(p.To-To_exp))/To_exp < 1e-6, [nm ': To profile uniform']);
    check(max(abs(p.Vx-Vx_exp))/Vx_exp < 1e-6, [nm ': Vx profile uniform']);
    check(max(abs(p.Vr))<1e-6 && max(abs(p.Vt))<1e-6, [nm ': Vr/Vt ~ 0']);
    mdot_exp = ro*Vx_exp*pi*(crd(end,2)^2-crd(1,2)^2); % full-annulus (fblade)
    check(abs(p.mdot-mdot_exp)/mdot_exp < 1e-3, ...
        [nm ': mdot matches annulus flux']);
    check(abs(p.To_mm-To_exp)/To_exp < 1e-6, [nm ': mass averaged To']);
    check(abs(sum(p.area_weights)-1) < 1e-12, [nm ': area weights normalised']);
end

check(res.PAMB==p_ref && res.TAMB==T_ref, 'PAMB/TAMB from inlet patch');
check(abs(res.RPTR-1) < 1e-6, 'RPTR = 1');
check(abs(res.TRRT) < 1e-6, 'TRRT = 0');
check(abs(res.TRRTSH) < 1e-6, 'TRRTSH = 0');
check(abs(res.MF-(res.M1-res.M2)) < 1e-12, 'MF = M1 - M2');

%% skewed (non-constant theta) periodic face
% Twist row 1's periodic faces by one full pitch hub to casing and build
% the flow as a ramp in face-relative theta: physics is identical to the
% untwisted ramp case, so L00S results must match. A fixed-theta target
% grid would sample off the cut support (nearest extrapolation) and bias
% the result by O(eps_r), so this verifies interp_plane's face tracking.
eps_r = 0.2;                                     % ramp amplitude
g_r  = dummy_grid(false);                        % ramp, constant theta face
g_sk = dummy_grid(false);                        % ramp, skewed face
for bid = 0:2
    c_tw = @(rr) 2*pi/NB*(rr-g_sk{1}.r(1,1,1))/(g_sk{1}.r(1,end,1)-g_sk{1}.r(1,1,1));
    b = g_r{bid+1};
    w = 1 + eps_r*(b.rt./b.r)/(2*pi/NB);
    b.ro = b.ro.*w; b.rovx = b.rovx.*w; b.roe = b.roe.*w;
    g_r{bid+1} = b;
    b = g_sk{bid+1};
    thl = b.rt./b.r;                     % face-relative theta once skewed
    w = 1 + eps_r*thl/(2*pi/NB);
    b.ro = b.ro.*w; b.rovx = b.rovx.*w; b.roe = b.roe.*w;
    b.rt = b.r.*(thl + c_tw(b.r));       % apply the face skew
    g_sk{bid+1} = b;
end
res_r  = ts3_extract(g_r);
[res_sk,dbg_sk] = ts3_extract(g_sk,true);
% the structured plane's k=1/nk lines must lie on the skewed face pair
th_plane = dbg_sk.L00S.plane.rt./dbg_sk.L00S.plane.r;
r_plane  = dbg_sk.L00S.plane.r(:,1);
check(max(abs(th_plane(:,1)-c_tw(r_plane))) < 1e-6, ...
    'skew face: plane k=1 lies on skewed face');
check(max(abs(th_plane(:,end)-(c_tw(r_plane)+2*pi/NB))) < 1e-6, ...
    'skew face: plane k=nk lies on image face');
check(max(abs(res_r.L00S.Po-Po_exp*(1+eps_r/2)))/Po_exp < 1e-3, ...
    'ramp: pitch averaged Po at mean state');
check(max(abs(res_sk.L00S.Po-res_r.L00S.Po))/Po_exp < 1e-3, ...
    'skew face: Po profile insensitive to face skew');
check(max(abs(res_sk.L00S.To-res_r.L00S.To))/To_exp < 1e-6, ...
    'skew face: To profile insensitive to face skew');
check(abs(res_sk.L00S.mdot-res_r.L00S.mdot)/res_r.L00S.mdot < 1e-3, ...
    'skew face: mdot insensitive to face skew');
check(abs(res_sk.L00S.To_mm-res_r.L00S.To_mm)/To_exp < 1e-6, ...
    'skew face: mass averaged To insensitive to face skew');

fprintf('\ntest_ts3_extract: all checks passed\n');

%% helper

function check(ok,msg)
assert(ok,'test_ts3_extract FAILED: %s',msg);
fprintf('  ok: %s\n',msg);
end

