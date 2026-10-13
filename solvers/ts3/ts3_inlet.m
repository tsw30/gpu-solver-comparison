function g = ts3_inlet(g,mdot_cor,P_in,T_in)
%ts3_inlet Apply inlet stagnation profiles matching a corrected mass flow.
%   g = ts3_inlet(g,mdot_cor)
%   g = ts3_inlet(g,mdot_cor,P_in,T_in)
%
%   Condensed replacement for ifan_inlet + ts_interpolate_inlet (originals
%   in ref/).  Sets pstag/tstag/yaw/pitch on kind-0 patches by
%   interpolating a radial profile, then scaling Mach number to match
%   mdot_cor [kg/s] at given P_in [Pa], T_in [K].
%
%   Defaults: P_in = 101325, T_in = 288.15

if nargin < 3, P_in = 101325; end
if nargin < 4, T_in = 288.15; end
cp = double(g{1}.av.cp);
ga = double(g{1}.av.ga);
Rgas = cp*(ga-1)/ga;

% Find first inlet patch for hub/casing radii
for bid = 0:numel(g)-1
    if ~isfield(g{bid+1},'patch'), continue; end
    for pid = 0:numel(g{bid+1}.patch)-1
        p = g{bid+1}.patch{pid+1}.attribute;
        if p.kind == 0
            ist = p.ist+1;
            rhub = double(g{bid+1}.r(ist,1,1));
            rcas = double(g{bid+1}.r(ist,end,1));
            break
        end
    end
    if exist('rhub','var'), break; end
end
if ~exist('rhub','var'), error('ts3_inlet: no inlet patch found'); end

% Precomputed spanwise Pstag distribution (LSF5 fan characteristic)
sdist = [0,0.0024211679,0.012394074,0.02236698,0.032339886,0.042312792,0.052285698,0.062258604,0.07223151,0.082204416,0.092177322,0.10215023,0.11218078,0.12215369,0.13212659,0.1420995,0.1520724,0.17201822,0.19196403,0.21190984,0.23185565,0.25180146,0.27174728,0.29169309,0.3116389,0.33164236,0.35158817,0.37153398,0.3914798,0.4114256,0.43137142,0.46129014,0.50118176,0.5410734,0.58102266,0.6308872,0.650833,0.7006975,0.74058915,0.76053496,0.78053844,0.80048424,0.82043004,0.84037584,0.8603217,0.8802675,0.9002133,0.9201591,0.9401049,0.9600507,0.97999656,1]';
Pvals = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.99931973,0.99863946,0.99727891,0.99047619,0.97891156,0.81836735];

% Solve for Mach scaling factor K to match target mass flow
rnorm = linspace(0,1,2000);
M0 = 0.458;  % constant Mach (average)
K0 = 1;
K = fzero(K_fun, K0, optimset('display','off', 'TolFun', 1e-8));

% Build inlet profile at 501 stations
r = linspace(0,1,501);
pstag = pchip(sdist,Pvals*P_in,r);
inlet.rnorm = r;
inlet.pstag = pstag;
inlet.tstag = repmat(T_in,size(r));
inlet.yaw = zeros(size(r));
inlet.pitch = zeros(size(r));

% Interpolate profiles onto each kind-0 patch (inline of ts_interpolate_inlet)
pmin = inf; pmax = -inf;
sic = [];
for bid = 0:numel(g)-1
    if ~isfield(g{bid+1},'patch'), continue; end
    for pid = 0:numel(g{bid+1}.patch)-1
        p = g{bid+1}.patch{pid+1}.attribute;
        if p.kind == 0
            ist = p.ist+1;
            rin = double(squeeze(g{bid+1}.r(ist,:,1)));
            if max(rin) < 0.5
                sic = [sic; bid, pid];
                pmin = min(pmin, min(rin));
                pmax = max(pmax, max(rin));
            end
        end
    end
end
for i = 1:size(sic,1)
    bid = sic(i,1); pid = sic(i,2);
    nk = double(g{bid+1}.attribute.nk);
    ist = g{bid+1}.patch{pid+1}.attribute.ist+1;
    rin = double(squeeze(g{bid+1}.r(ist,:,1)));
    rn = (rin - pmin) ./ (pmax - pmin);
    g{bid+1}.patch{pid+1}.pstag = single(repmat(interp1(r,inlet.pstag,rn,'pchip','extrap'),nk,1));
    g{bid+1}.patch{pid+1}.tstag = single(repmat(T_in*ones(size(rn)),nk,1));
    g{bid+1}.patch{pid+1}.yaw   = single(zeros(nk*numel(rn),1));
    g{bid+1}.patch{pid+1}.pitch = single(zeros(nk*numel(rn),1));
    g{bid+1}.patch{pid+1}.rfin  = single(0.25);
    g{bid+1}.patch{pid+1}.sfinlet = int32(0);
end
fprintf('ts3_inlet: applied (mdot=%.4f kg/s, P_in=%.1f, T_in=%.1f)\n',mdot_cor,P_in,T_in);

%% Nested function for mass flow matching
    function dm = K_fun(K)
        M0 = 0.458*K;
        M = M0;
        Ploc = pchip(sdist,Pvals*P_in,rnorm);
        Pstat = P_in ./ (1 + (ga-1)/2*M.^2).^(ga/(ga-1));
        Q = ga/sqrt(ga-1)*M.*(1+(ga-1)/2*M.^2).^(-0.5*(ga+1)/(ga-1));
        rs = rnorm*(rcas-rhub)+rhub;
        dA = pi*(rs(2:end).^2 - rs(1:end-1).^2);
        Pm = 0.5*(Ploc(2:end)+Ploc(1:end-1));
        Qm = 0.5*Q(2:end) + 0.5*Q(1:end-1);
        m = 1/sqrt(cp*T_in)*sum(Pm.*Qm.*dA);
        dm = mdot_cor - m;
    end
end