function [im,smap] = sp3drecon(varargin)
% Input arguments:
%   ro_delay: kspace sampling delay (can be a fraction), default = 0
%   nclip: number of samples to clip from beginning and end of readout,
%       default = 50
%   smap: complex coil sensitivity maps (dimensions must be consistent),
%       default = []
%   niter: number of conjugate gradient iterations (0 = conjugate phase),
%       default = 0
%   compfrac: coil compression factor (fraction < 1), default = 1

% Set default arguments
arg.ro_delay = 0;
arg.nclip = 50;
arg.niter = 0;
arg.compfrac = 1;
arg.smap = [];

% Parse arguments
arg = toppe.utils.vararg_pair(arg, varargin);

% Load in sequence parameters & system specs
fprintf('loading in sequence parameters & systems...\n')
load seq.mat seq
sys = toppe.systemspecs('maxSlew', 20);

% Read in pfile
fprintf('reading in pfile...\n');
pfile = dir('./P*.7');
raw = toppe.utils.loadpfile(pfile(1).name);
raw = raw(end:-1:1,:,:,:); % first dimension is flipped
ndat = size(raw,1);
ncoils = size(raw,2);

% Read in the initial kspace trajectory
fprintf('reading in initial kspace trajectory...\n');
[~,gx_sp0,gy_sp0,gz_sp0] = toppe.readmod('readout.mod');
g_sp0 = [gx_sp0';gy_sp0';gz_sp0'];
k_sp0 = sys.gamma*1e-4*cumsum(g_sp0,2)*sys.raster*1e-6;

% Order views and data
fprintf('ordering views and data...\n');
kdata = zeros(ndat,seq.nechoes*seq.nshots,seq.nframes,ncoils);
klocs = zeros(size(k_sp0,2),seq.nechoes*seq.nshots,3);
scanloop = importdata('scanloop.txt','\t',3);
scanloop = scanloop.data;
rotmatrices = scanloop(scanloop(:,10) == 1, end-11:end-3) / 32767;
for framen = 1:seq.nframes
    for shotn = 1:seq.nshots
        for echon = 1:seq.nechoes
            R = ones(3);
            R(1,:) = rotmatrices((framen-1)*seq.nshots*seq.nechoes + (shotn-1)*seq.nechoes + echon,1:3);
            R(2,:) = rotmatrices((framen-1)*seq.nshots*seq.nechoes + (shotn-1)*seq.nechoes + echon,4:6);
            R(3,:) = rotmatrices((framen-1)*seq.nshots*seq.nechoes + (shotn-1)*seq.nechoes + echon,7:9);
            kdata(:,(shotn-1)*seq.nechoes + echon,framen,:) = permute(raw(:,:,1,(shotn-1)*seq.nechoes + echon),[1,3,4,2]);
            klocs(:,(shotn-1)*seq.nechoes + echon,:) = permute(R*k_sp0,[2,3,1]);
        end
    end
end

% Apply sampling delay and clip end points
klocs = interp1(1:size(klocs,1), klocs, (arg.nclip+1:size(kdata,1)-arg.nclip)+arg.ro_delay);
kdata = kdata(arg.nclip+1:end-arg.nclip,:,:,:);

% Reshape data and recon
fprintf('reshaping data and reconning...\n');
klocs = reshape(klocs,[],3);
kdata = reshape(kdata,[],seq.nframes,ncoils);
[im,smap] = reconutils.recon_cgsense(klocs(:,1:seq.ndims), kdata, ...
    seq.N*ones(1,seq.ndims), seq.fov*ones(1,seq.ndims), ...
    'compfrac', arg.compfrac,'niter',arg.niter,'smap',arg.smap);

end