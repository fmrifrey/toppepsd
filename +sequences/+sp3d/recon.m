function im = sp3dgre_recon(varargin)
% Input arguments:
%   compfrac: coil compression factor (fraction < 1)
%   overwritesmap: option to calculate new sensitivity map and write out 
%       as a .mat file
%   ro_delay: number of samples to delay the readout from gradients
%   niter: number of iterations for pcg

%% Set default arguments
arg.compfrac = 0.25;
arg.overwritesmap = 0;
arg.ro_delay = 10;
arg.niter = 0;
arg.maxframe = 'all';
arg = toppe.utils.vararg_pair(arg, varargin);

%% Read in pfile
fprintf('reading in pfile...\n');
pfile = dir('./P*.7');
raw = toppe.utils.loadpfile(pfile(1).name);
raw = raw(end:-1:1,:,:); % first dimension is flipped
raw = permute(raw,[1,3,2]);
ncoils = ceil(arg.compfrac*size(raw,3));

% Load in sequence parameters & system specs
fprintf('loading in sequence parameters & systems...\n')
load seq.mat seq
sys = toppe.systemspecs('maxSlew', 20);

% Compress coils
fprintf('compressing coils...\n')
if ncoils > 1 && arg.compfrac < 1
    raw = ir_mri_coil_compress(raw,'ncoil',ncoils);
end

% Read in the initial kspace trajectory
fprintf('reading in readout gradients...\n')
[~,gx_sp0,gy_sp0,gz_sp0] = toppe.readmod('readout.mod');
g_sp0 = [gx_sp0';gy_sp0';gz_sp0'];
k_sp0 = sys.gamma*1e-4*cumsum(g_sp0(:,arg.ro_delay+1:size(raw,1)+arg.ro_delay),2)*sys.raster*1e-6;

% Read in rotation matrices and get the entire trajectory
fprintf('transforming trajectories...\n');
scanloop = importdata('scanloop.txt','\t',3);
scanloop = scanloop.data;
rotmatrices = scanloop(scanloop(:,1) == 2, end-11:end-3) / 32767;
klocs = [];
for i = 1:seq.nechoes*seq.nshots
    R = ones(3);
    R(1,:) = rotmatrices(i,1:3);
    R(2,:) = rotmatrices(i,4:6);
    R(3,:) = rotmatrices(i,7:9);
    klocs = [klocs, R*k_sp0];
end

% Create Gmri object and density compensation
fprintf('creating Gmri object...\n')
nufft_args = {seq.N*ones(1,3),...
    6*ones(1,3),...
    2*seq.N*ones(1,3),...
    seq.N/2*ones(1,3),...
    'table',...
    2^10,...
    'minmax:kb'};
Gm = Gmri(klocs', true(seq.N*ones(1,3)), ...
    'fov', seq.D, 'basis', {'rect'}, 'nufft', nufft_args(:)');
dcf = pipe_menon_dcf(Gm.Gnufft);
W = Gdiag(dcf(:));


% Make regularizer
R = Reg1(ones(seq.N*ones(1,3)), 'beta', 2^-12 * numel(klocs)/3, ...
    'mask', true(seq.N*ones(1,3)));
C = R.C;

% Load sensitivity map from file
if isfile('smap.mat')
    load smap.mat smap
end
if ncoils == 1
    smap = ones([seq.N*ones(1,3),ncoils]);
elseif ~isfile('smap.mat') || size(smap,4)~=ncoils || arg.overwritesmap
    
    % Loop through coils & recon
    fprintf('reconning coil by coil images...\n');
    imc = zeros([seq.N*ones(1,3),ncoils]);
    for i = 1:ncoils
        fidx = 1; % use frame 1
        vidx = 1:seq.nshots*seq.nechoes;

        % Get data for current coil
        data = reshape( raw(:,(fidx-1)*length(vidx) + vidx,i) ,[],1);
        
        % Recon the data
        imc(:,:,:,i) = recondata(Gm,data,W,arg.niter,C,seq.N*ones(1,3));
    end
    
    % Estimate sensitivity map using bart
    fprintf('estimating sensitivity map using bart...\n');
    imc = reshape(imc,[seq.N*ones(1,3),ncoils]);
    smap = bart('ecalib -b0 -m1',fftc(imc,1:3));

    % Save the sensitivity map
    save smap.mat smap
end

%% Incorporate sensitivity map and reshape density compensation
if ncoils > 1
    fprintf('creating new Gmri object w/ SENSE encoding...\n');
    Gm = Asense(Gm,smap);
    W = Gdiag(repmat(dcf(:),1,ncoils));
end

%% Loop through frames and recon
if strcmpi(arg.maxframe, 'all')
    arg.maxframe = seq.nframes;
end
fprintf('reconning frame by frame images...\n');
im = zeros([seq.N*ones(1,3),arg.maxframe]);
for i = 1:arg.maxframe
    % Calculate view indices
    vidx = 1:seq.nshots*seq.nechoes;

    % Get data for current frame
    data = reshape(raw(:,(i-1)*length(vidx) + vidx,:),[],1);
    
    % Recon preconditioner using conjugate-phase
    im(:,:,:,i) = recondata(Gm,data,W,arg.niter,C,seq.N*ones(1,3));
end

fprintf('all done!\n');

end

function im = recondata(Gm,data,W,niter,C,dim)

    % Recon preconditioner using conjugate-phase
    im_cp = Gm' * reshape(W * data(:), [], 1);
    im_cp = ir_wls_init_scale(Gm, data(:), im_cp);
    im_cp = embed(im_cp,true(dim));

    % Recon using preconditioned conjugate gradient (iterative)
    if niter > 0
        im_pcg = qpwls_pcg1(im_cp(true(dim)), Gm, 1, data, C, ...
            'niter', niter,'isave', 1:niter); 
        im = embed(im_pcg(:,1),true(dim));
        
    else % ...or save image with CP recon
        im = im_cp;
    end

end

function Wi = pipe_menon_dcf(G,itrmax)

    % Set default for itrmax
    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end
    
    % If G is a Gmri object, use its Gnufft object
    if isfield(G.arg,'Gnufft')
        G = G.Gnufft;
    end
    
    % Initialize weights to 1 (psf)
    Wi = ones(size(G,1),1);
    
    % Loop through iterations
    for itr = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( G.arg.st.interp_table(G.arg.st, ...
            G.arg.st.interp_table_adj(G.arg.st, Wi) ) );
        Wi = Wi ./ d;
        
    end
    
    % Normalize weights
    Wi = Wi / sum(abs(Wi));
    
end