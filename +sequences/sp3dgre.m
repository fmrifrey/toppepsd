function sp3dgre(varargin)
% Input arguments:
%   tipdown_fa: tipdown flip angle (deg), default = 90
%   tipdown_slthick: tipdown slab thickness (cm), default = 20
%   tipdown_tbw: tipdown time bandwidth product, default = 16
%   tipdown_dur: tipdown duration (us), default = 3200
%   tipdown_nchop: number of dead samples before and after
%       tipdown, default = [25, 0]
%   tipdown_ftype: slr pulse type, default = 'min
%   tipdown_rf: tipdown rf waveform (mG), leave empty for SLR
%   tipdown_gz: tipdown gz waveform (G/cm)
%
%   spiral_F: FOV coefficients for VD spiral (see brian hargraeve's
%       vds paper for more info), default = 24 (cm)
%   spiral_kxymax: maximum in-plane kspace radius (cm^-1), default
%       = 1.33 (cm^-1)
%   spiral_kzmax: maximum through-plane kspace radius (cm^-1),
%       default = 0 (cm^-1)
%   spiral_dir: spiral direction ('i' = in, 'o' = out, 'io' =
%       in-out, 'oi' = out-in), default = 'o'
%   spiral_nshots: number of shots in spiral, default = 1
%   spiral_nnav: number of navigator points in spiral, default =
%       100
%
%   crush_ncycles: number of cycles of phase across tipdown slab
%

% Set up sequence
sys = toppe.systemspecs('maxSlew',20);
save sys
toppe.writeentryfile('toppeN.entry');

% Set default arguments
arg.fov = 16;
arg.N = 128;
arg.tipdown_fa = 70;
arg.tipdown_slthick = 1;
arg.tipdown_tbw = 16;
arg.tipdown_dur = 3200;
arg.tipdown_nchop = [25 0];
arg.tipdown_ftype = 'min';
arg.spiral_kxymax = arg.N/arg.fov/2;
arg.spiral_kzmax = 0;
arg.spiral_dir = 'o';
arg.spiral_nnav = 100;
arg.crush_ncycles = 6;
arg.padtime = 1000;
arg.nframes = 1;
arg.nshots = 1;
arg.nechoes = 128;
arg.ndd = 0;
arg.TE = 'min';
arg.TR = 50;
arg.Tdelay = 10;
arg.ndims = 2;

% Parse arguments
arg = toppe.utils.vararg_pair(arg, varargin);

save arg

% Set tipdown type
if arg.tipdown_fa < 30
    tipdown_type = 'st';
else
    tipdown_type = 'ex';
end

% Create SLR pulse
[tipdown_rf, tipdown_gz, tipdown_freq] = toppe.utils.rf.makeslr( ...
    arg.tipdown_fa, ... % deg
    arg.tipdown_slthick, ... % cm
    arg.tipdown_tbw, ...
    arg.tipdown_dur*1e-3, ... % ms
    0, ... % no spoiling
    sys, ...
    'type', tipdown_type, ...
    'discardPrephaser', true, ...
    'ftype', arg.tipdown_ftype, ...
    'writeModFile', false ...
    );

% Write tipdown out to mod file
toppe.writemod(sys, ...
    'rf', tipdown_rf, ...
    'gz', tipdown_gz, ...
    'nChop', arg.tipdown_nchop, ...
    'ofname', 'tipdown.mod');

% Create initial spiral trajectory
g_sp0 = psdutils.spiral.spgrad(sys, arg.spiral_dir, ...
    arg.fov*[1,0,0], arg.spiral_kxymax, arg.spiral_kzmax, ...
    arg.nechoes, arg.spiral_nnav);

% Write readout out to mod file
toppe.writemod(sys, ...
    'gx', g_sp0(1,:)', ...
    'gy', g_sp0(2,:)', ...
    'gz', g_sp0(3,:)', ...
    'nChop', [10,10], ... % needs enough time for adc to catch up
    'ofname', 'readout.mod');

% Create the crusher
gspoil = toppe.utils.makecrusher(arg.crush_ncycles, ...
    arg.tipdown_slthick, ...
    sys, ...
    0, ...
    sys.maxSlew, ...
    sys.maxGrad);

% Write out to mod file
toppe.writemod(sys, ...
    'gx', 0*gspoil, ...
    'gy', 0*gspoil, ...
    'gz', gspoil, ...
    'ofname', 'crusher.mod');

% Calculate all module durations
dur_tipdown = length(toppe.readmod('tipdown.mod'))*sys.raster + arg.padtime; % us
dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + arg.padtime; % us
dur_crusher = length(toppe.readmod('crusher.mod'))*sys.raster + arg.padtime; % us

% Set cores file text
coresfiletext = [ ...
    "2 1 0" % tipdown, delay
    "1 2" % readout
    "2 0 3" % delay, crush
    "1 0" % delay
    ];

% Set modules file text
modulesfiletext = [ ...
    sprintf("tipdown.mod %d 1 0 -1", dur_tipdown)
    sprintf("readout.mod %d 0 1 -1", dur_readout)
    sprintf("crusher.mod %d 0 0 -1", dur_crusher)
    ];

% Write out cores and modules files
psdutils.filegen.writemodulesfile(modulesfiletext);
psdutils.filegen.writecoresfile(coresfiletext);

% Set minimum TE
minTE = 0;
if strcmpi(arg.TE, 'min')
    arg.TE = minTE;
elseif arg.TE < minTE
    error('TE must be >= %fms\n', minTE);
end

% Set minimum TR
minTR = 1e-3*(dur_tipdown + dur_readout + ...
    dur_crusher) + arg.TE; % ms
if strcmpi(arg.TR, 'min')
    arg.TR = minTR;
elseif arg.TE < minTE
    error('TR must be >= %fms\n', minTR);
end

% Loop through frames, shots, echoes
toppe.write2loop('setup',sys,'version',6);
for framen = 1:arg.nframes
    for shotn = 1:arg.nshots
        for echon = 1:arg.nechoes

            viewn = (framen-1)*arg.nshots*arg.nechoes + (shotn-1)*arg.nechoes + echon;

            if echon > arg.ndd % non-disdaq echoes
                dabmode = 'on';
                gamp = ones(3,1);
                R = eul2rotm([((shotn-1)*arg.nechoes + echon)*pi*(3 - sqrt(5))/2,0,0],"ZYX");
            else
                dabmode = 'off';
                gamp = zeros(3,1);
                R = eye(3);
            end

            % Write tipdown to loop
            toppe.write2loop('tipdown.mod', sys, ...
                'RFoffset', tipdown_freq, ...
                'echo', 1, ...
                'slice', 1, ...
                'view', viewn, ...
                'RFspoil', 1, ...
                'version', 6, ...
                'trigout', 0, ...
                'core', 1);

            % Write echo time deadtime to loop
            toppe.write2loop('delay', sys, ...
                'textra', arg.TE, ... % ms
                'core', 1);

            % Write readout to loop
            toppe.write2loop('readout.mod', sys, ...
                'echo', 1, ...
                'slice', 1, ...
                'view', viewn, ...
                'RFspoil', 1, ...
                'version', 6, ...
                'trigout', 0, ...
                'dabmode', dabmode, ...
                'Gamplitude', gamp, ...
                'rotmat', R, ...
                'core', 2);

            % Write repetition time deadtime to loop
            toppe.write2loop('delay', sys, ...
                'textra', arg.TR - 1e-3*(dur_tipdown + dur_readout) - arg.TE, ...
                'core', 3);

            % Crusher
            toppe.write2loop('crusher.mod', sys, ...
                'echo', 1, ...
                'slice', 1, ...
                'view', viewn, ...
                'version', 6, ...
                'core', 3);
        end
            
        % Write repetition time deadtime to loop
        toppe.write2loop('delay', sys, ...
            'textra', arg.Tdelay, ...
            'core', 4);

    end
end

% Finish loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry','seqstamp.txt',sys);

end