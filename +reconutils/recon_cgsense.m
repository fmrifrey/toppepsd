function [im,smap] = recon_cgsense(klocs,kdata,dim,fov,varargin)
% Input arguments:
%   klocs: kspace sampling locations (cm^-1, [ndat x ndim])
%   kdata: complex data at sampling locations ([ndat x nframes x ncoils])
%   dim: image dimensions ([ndim x 1] vector)
%   fov: image fov (cm, [ndim x 1] vector)
%   smap: complex coil sensitivity maps (dimensions must be consistent),
%       default = []
%   niter: number of conjugate gradient iterations (0 = conjugate phase),
%       default = 0
%   compfrac: coil compression factor (fraction < 1), default = 1

% Set default arguments
arg.smap = [];
arg.niter = 0;
arg.compfrac = 1;

% Parse arguments
arg = toppe.utils.vararg_pair(arg, varargin);

% Get dimensions
ncoils = ceil(arg.compfrac*size(kdata,3));
nframes = size(kdata,2);
ndim = length(dim);

% Compress coils
fprintf('compressing coils...\n')
if size(kdata,3) > 1 && arg.compfrac < 1
    kdata = ir_mri_coil_compress(kdata,'ncoil',ncoils);
end

if (ncoils > 1 && isempty(arg.smap)) || (size(arg.smap,ndim+1) ~= ncoils )
    warning('no valid sense map passed, creating new one...');
    % Recurse with coil-wise images as frames
    kdatac = squeeze(kdata(:,1,:));
    imc = reconutils.recon_cgsense(klocs,kdatac,dim,fov,varargin{:});

    if ndims(imc) < 4 % if 2D
        imc = permute(imc, [1,2,4,3]);
    end

    % Create sense map using bart
    smap = bart('ecalib -b0 -m1', reconutils.fftc(imc,1:3));
elseif ncoils == 1
    smap = ones(dim); % for single coil, uniform sensitivity
else
    smap = arg.smap; % use passed smap
end

% Create Gmri object and density compensation
fprintf('creating Gmri object...\n')
nufft_args = {dim,...
    6*ones(1,length(dim)),...
    2*dim,...
    1/2*dim,...
    'table',...
    2^10,...
    'minmax:kb'};
Gm = Gmri(klocs, true(dim), ...
    'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args(:)');
dcf = reconutils.pipe_menon_dcf(Gm.Gnufft);
W = Gdiag(repmat(dcf(:),1,ncoils));
if ncoils > 1
    Gm = Asense(Gm,smap);
end

% Make regularizer
R = Reg1(ones(dim), 'beta', 2^-12 * numel(klocs)/3, ...
    'mask', true(dim));
C = R.C;

% Loop through frames and recon
im = zeros(prod(dim),nframes);
for i = 1:nframes

    kdataf = squeeze(kdata(:,i,:));

    % Recon preconditioner using conjugate-phase
    im_cp = Gm' * reshape(W * kdataf, [], 1);
    im_cp = ir_wls_init_scale(Gm, kdataf(:), im_cp);
        
    % Recon using preconditioned conjugate gradient (iterative)
    if arg.niter > 0
        im_pcg = qpwls_pcg1(im_cp(true(dim)), Gm, W, kdataf(:), C, ...
            'niter', niter,'isave', 1:niter);
        im(:,i) = im_pcg(:,1);
    else % ...or save image with CP recon
        im(:,i) = im_cp;
    end

end

% Reshape image
im = reshape(im,[dim,nframes]);

end
