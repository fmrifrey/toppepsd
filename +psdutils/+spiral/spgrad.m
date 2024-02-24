function g = spgrad(sys,spdir,F,kxymax,kzmax,nshots,nnav)
% Input arguments:
%   sys: toppe systemspects structure
%   spdir: spiral direction ('in'=in, 'out'=out, 'io'=in-out, 'oi'=out-in)
%   F: VD spiral fov coefficients ([cm, cm^2, cm^3])
%   kxymax: maximum kxy radius (cm^-1)
%   kzmax: maximum kz radius (cm^-1)
%   nshots: number of shots for vds
%   nnav: number of navigator points (at beginning for in & out, in middle
%       for in-out & out-in)
%   
% Output arguments:
%   g = gradient waveform ([x;y;z] G/cm)
%

gam = 4258; % gyromagnetic ratio (Hz/G)

% Generate the xy spiral out trajectory
if length(spdir) > 1
    Nleaves = 2;
else
    Nleaves = 1;
end
[~,gxy_spout] = toppe.utils.spiral.mintgrad.vds(sys.maxSlew*1e3, ...
    sys.maxGrad, sys.raster*1e-6, nshots*Nleaves, F(1), F(2), F(3), kxymax); % 1/cm
gxy_spout = [real(gxy_spout);imag(gxy_spout)];

% Ramp down the xy gradients to 0
gxy_spout = [gxy_spout, ...
    linspace(1,0,ceil(norm(gxy_spout(:,end),2)/sys.maxSlew*1e-3/sys.raster*1e6)).*gxy_spout(:,end)];

% Rewind the spiral to the center of kspace
kxy0 = gam*sys.raster*1e-6*sum(gxy_spout,2);
gxy_spout = [gxy_spout, ...
    psdutils.generic.opttrap(norm(kxy0,2)/gam, sys.maxGrad, sys.maxSlew*1e3, ...
    sys.raster*1e-6) .* -kxy0 / norm(kxy0,2)]; % G/cm

% Create kz ramps
if kzmax > 0
    gz_ramp = psdutils.generic.opttrap(kzmax/gam, sys.maxGrad, sys.maxSlew*1e3, sys.raster*1e-6);
else
    gz_ramp = [];
end

% Begin gradient with kz ramp-out waveform
g = [zeros(2,size(gz_ramp,2)); gz_ramp];

% Append gradient with spiral kspace trajectory
switch spdir
    case 'o' % spiral out
        g = [g, ...
            zeros(3,nnav), ... % navigators
            [gxy_spout; zeros(1,size(gxy_spout,2))] ... % spiral out
            ];
    case 'i' % spiral in
        g = [g, ...
            zeros(3,nnav), ... % navigators
            [gxy_spout(:,end:-1:1); zeros(1,size(gxy_spout,2))] ... % spiral in
            ];
    case 'oi' % spiral out-in
        g = [g, ...
            [gxy_spout; zeros(1,size(gxy_spout,2))], ... % spiral out
            zeros(3,nnav), ... % navigators
            [gxy_spout(:,end:-1:1); zeros(1,size(gxy_spout,2))] ... % spiral in
            ];
    case 'io' % spiral in-out
        g = [g, ...
            [gxy_spout(:,end:-1:1); zeros(1,size(gxy_spout,2))], ... % spiral in
            zeros(3,nnav), ... % navigators
            [gxy_spout; zeros(1,size(gxy_spout,2))] ... % spiral out
            ];
end

% Append gradient with kz ramp-in waveform
g = [g, [zeros(2,size(gz_ramp,2)); -gz_ramp]];

end