function sp3dfse(varargin)
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
%   refocus_fa: refocuser flip angle (deg), default = 180
%   refocus_slthick: refocuser slab thickness (cm), default = 20
%   refocus_tbw: refocuser time bandwidth product, default = 2
%   refocus_dur: refocuser duration (us), default = 3200
%   refocus_nchop: number of dead samples before and after
%       refocuser, default = [25, 25]
%   refocus_ftype: slr pulse type, default = 'ls'
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
arg.tipdown_fa = 90;
arg.tipdown_slthick = 20;
arg.tipdown_tbw = 8;
arg.tipdown_dur = 3200;
arg.tipdown_nchop = [25 0];
arg.tipdown_ftype = 'min';
arg.spiral_F = [20, 0, 0];
arg.spiral_kxymax = 64/20/2;
arg.spiral_kzmax = 0;
arg.spiral_dir = 'io';
arg.spiral_nnav = 100;
arg.refocus_fa = 120;
arg.refocus_slthick = 1;
arg.refocus_tbw = 4;
arg.refocus_dur = 3200;
arg.refocus_nchop = [25, 25];
arg.crush_ncycles = 6;
arg.padtime = 1000;
arg.nframes = 1;
arg.nshots = 64;
arg.nechoes = 1;
arg.ndd = 0;
arg.TE = 'min';
arg.TR = 'min';
arg.Tdelay = 50;

% Parse arguments
arg = toppe.utils.vararg_pair(arg, varargin);
seq = arg;
save seq

% Set tipdown type
if arg.tipdown_fa < 30
    tipdown_type = 'st';
else
    tipdown_type = 'ex';
end

% Create tipdown SLR pulse
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

% Create refocuser SLR pulse
[refocus_rf, refocus_gz, refocus_freq] = toppe.utils.rf.makeslr( ...
    arg.refocus_fa, ... % deg
    arg.refocus_slthick, ... % cm
    arg.refocus_tbw, ...
    arg.refocus_dur*1e-3, ... % ms
    arg.crush_ncycles, ...
    sys, ...
    'type', 'se', ...
    'writeModFile', false ...
    );

% Write refocuser out to mod file
toppe.writemod(sys, ...
    'rf', refocus_rf, ...
    'gz', refocus_gz, ...
    'nChop', arg.refocus_nchop, ...
    'ofname', 'refocus.mod');

% Create initial spiral trajectory
g_sp0 = psdutils.spiral.spgrad(sys, arg.spiral_dir, ...
    arg.spiral_F, arg.spiral_kxymax, arg.spiral_kzmax, ...
    arg.nshots, arg.spiral_nnav);

% Write readout out to mod file
toppe.writemod(sys, ...
    'gx', g_sp0(1,:)', ...
    'gy', g_sp0(2,:)', ...
    'gz', g_sp0(3,:)', ...
    'nChop', [10,10], ... % needs enough time for adc to catch up
    'ofname', 'readout.mod');

% Calculate all module durations
dur_tipdown = length(toppe.readmod('tipdown.mod'))*sys.raster + arg.padtime; % us
dur_refocus = length(toppe.readmod('refocus.mod'))*sys.raster + arg.padtime; % us
dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + arg.padtime; % us

% Set cores file text
coresfiletext = [ ...
    "2 1 0" % tipdown, delay
    "2 2 0" % refocus, delay
    "2 3 0" % readout, delay
    "1 0" % delay
    ];

% Set modules file text
modulesfiletext = [ ...
    sprintf("tipdown.mod %d 1 0 -1", dur_tipdown)
    sprintf("refocus.mod %d 1 0 -1", dur_refocus)
    sprintf("readout.mod %d 0 1 -1", dur_readout)
    ];

% Write out cores and modules files
psdutils.filegen.writemodulesfile(modulesfiletext);
psdutils.filegen.writecoresfile(coresfiletext);

% Set minimum TE
minTE = 1e-3*(dur_refocus + dur_readout);
if strcmpi(arg.TE, 'min')
    arg.TE = minTE;
elseif arg.TE < minTE
    error('TE must be >= %fms\n', minTE);
end

% Golden ratios
phi1 = 0.4656; % 2D golden ratio 1
phi2 = 0.6823; % 2D golden ratio 2

% Loop through frames, shots, echoes
toppe.write2loop('setup',sys,'version',6);
for framen = 1:arg.nframes
    for shotn = 1:arg.nshots

        % Write tipdown to loop
        toppe.write2loop('tipdown.mod', sys, ...
            'RFoffset', tipdown_freq, ...
            'RFphase', 0, ...
            'echo', 1, ...
            'slice', 1, ...
            'view', 1, ...
            'RFspoil', 0, ...
            'version', 6, ...
            'trigout', 0, ...
            'core', 1);

        % Write time delay to loop
        toppe.write2loop('delay', sys, ...
            'textra', (arg.TE - dur_refocus*1e-3)/2, ...
            'core', 1);

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
            toppe.write2loop('refocus.mod', sys, ...
                'RFoffset', refocus_freq, ...
                'RFphase', pi/2, ...
                'echo', 1, ...
                'slice', 1, ...
                'view', viewn, ...
                'RFspoil', 0, ...
                'version', 6, ...
                'trigout', 0, ...
                'core', 2);

            % Write echo time deadtime to loop
            toppe.write2loop('delay', sys, ...
                'textra', (arg.TE - minTE)/2, ... % ms
                'core', 2);

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
                'core', 3);

            % Write repetition time deadtime to loop
            toppe.write2loop('delay', sys, ...
                'textra', (arg.TE - minTE)/2, ...
                'rotmat', R, ...
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