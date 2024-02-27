% Load in sequence parameters & system specs
fprintf('loading in sequence parameters & systems...\n')
load arg.mat arg
sys = toppe.systemspecs('maxSlew', 20);

% Read in pfile
fprintf('reading in pfile...\n');
pfile = dir('./P*.7');
raw = toppe.utils.loadpfile(pfile(1).name);
raw = raw(end:-1:1,:,:,:); % first dimension is flipped
ndat = size(raw,1);
ncoils = size(raw,2);
nviews = size(raw,4);
ndims = 2;

% Read in the initial kspace trajectory
[~,gx_sp0,gy_sp0,gz_sp0] = toppe.readmod('readout.mod');
g_sp0 = [gx_sp0';gy_sp0';gz_sp0'];
k_sp0 = sys.gamma*1e-4*cumsum(g_sp0,2)*sys.raster*1e-6;

kdata = zeros(ndat,arg.nechoes*arg.nshots,arg.nframes,ncoils);
klocs = zeros(size(k_sp0,2),arg.nechoes*arg.nshots,3);
for framen = 1:arg.nframes
    for shotn = 1:arg.nshots
        for echon = 1:arg.nechoes
            R = eul2rotm([((shotn-1)*arg.nechoes + echon)*pi*(3 - sqrt(5))/2,0,0],"ZYX");
            kdata(:,(shotn-1)*arg.nechoes + echon,framen,:) = permute(raw(:,:,1,(shotn-1)*arg.nechoes + echon),[1,3,4,2]);
            klocs(:,(shotn-1)*arg.nechoes + echon,:) = permute(R*k_sp0,[2,3,1]);
        end
    end
end

% Apply readout delay
ro_delay = 10;
klocs = klocs((1:size(kdata,1)) + ro_delay,:,:);

% Recon
[im,smap] = reconutils.recon_cgsense(klocs(:,1:ndims), kdata, arg.N*ones(1,ndims), arg.fov*ones(1,ndims), 'compfrac', 1);
lbview(im);