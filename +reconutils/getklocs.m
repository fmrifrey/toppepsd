function klocs = getklocs(sys)

if nargin < 1 || isempty(sys)
    sys = toppe.systemspecs;
end

% Read in the initial kspace trajectory
[~,gx_sp0,gy_sp0,gz_sp0] = toppe.readmod('readout.mod');
g_sp0 = [gx_sp0';gy_sp0';gz_sp0'];
k_sp0 = sys.gamma*1e-4*cumsum(g_sp0,2)*sys.raster*1e-6;

% Read in rotation matrices and get the entire trajectory
scanloop = importdata('scanloop.txt','\t',3);
scanloop = scanloop.data;
rotmatrices = scanloop(scanloop(:,10) == 1, end-11:end-3) / 32767;
klocs = zeros(size(k_sp0,2),size(rotmatrices,1),3);
for i = 1:size(rotmatrices,1)
    R = ones(3);
    R(1,:) = rotmatrices(i,1:3);
    R(2,:) = rotmatrices(i,4:6);
    R(3,:) = rotmatrices(i,7:9);
    klocs(:,i,:) = permute((R*k_sp0)',[1,3,2]);
end

end

