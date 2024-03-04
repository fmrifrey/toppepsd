function sp3dflex(varargin)
% Input arguments:
%   N: image matrix size, default = 128
%   ndims: number of dimensions in matrix size, default = 2
%   fov: image isotropic field of view
%   ro_mode: readout mode ('gre' or 'fse'), default = 'gre'
%   nframes:
%   nshots:
%   nechoes:
%
%   ndd:
%   TE:
%   TR:
%   Tdelay:
%   padtime:
%
%   tipdown_fa: tipdown flip angle (deg), default = 90
%   tipdown_slthick: tipdown slab thickness (cm), default = 20
%   tipdown_tbw: tipdown time bandwidth product, default = 16
%   tipdown_dur: tipdown duration (us), default = 3200
%   tipdown_ftype: slr pulse type, default = 'min
%   tipdown_rf: tipdown rf waveform (mG), leave empty for SLR
%   tipdown_gz: tipdown gz waveform (G/cm)
%
%   refocus_fa: refocuser flip angle (deg), default = 180
%   refocus_slthick: refocuser slab thickness (cm), default = 20
%   refocus_tbw: refocuser time bandwidth product, default = 2
%   refocus_dur: refocuser duration (us), default = 3200
%   refocus_ftype: slr pulse type, default = 'ls'
%
%   spiral_F: FOV coefficients for VD spiral as a fraction of nominal FOV 
%       (see brian hargraeve's vds paper for more info), default = [1,0,0]
%   spiral_dir: spiral direction ('i' = in, 'o' = out, 'io' =
%       in-out, 'oi' = out-in), default = 'o'
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
arg.fov = 20;
arg.N = 64;
arg.ndims = 2;
arg.ro_mode = 'gre';
arg.nframes = 1;
arg.nshots = 64;
arg.nechoes = 1;

arg.ndd = 0;
arg.TE = 'min';
arg.TR = 50;
arg.Tdelay = 500;
arg.padtime = 1000;

arg.tipdown_fa = 70;
arg.tipdown_slthick = 1;
arg.tipdown_tbw = 8;
arg.tipdown_dur = 3200;
arg.tipdown_ftype = 'min';

arg.refocus_fa = 120;
arg.refocus_slthick = 1;
arg.refocus_tbw = 4;
arg.refocus_dur = 3200;

arg.spiral_dir = 'o';
arg.spiral_nnav = 100;
arg.spiral_F = [1,0,0];

arg.crush_ncycles = 6;

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
    'ofname', 'refocus.mod');

% Create initial spiral trajectory
g_sp0 = psdutils.spiral.spgrad(sys, arg.spiral_dir, ...
    arg.fov*arg.spiral_F, arg.N/arg.fov/2, 0, ...
    arg.nshots, arg.spiral_nnav);

% Write readout out to mod file
toppe.writemod(sys, ...
    'gx', g_sp0(1,:)', ...
    'gy', g_sp0(2,:)', ...
    'gz', g_sp0(3,:)', ...
    'ofname', 'readout.mod');

% Create the crusher for GRE
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
dur_refocus = length(toppe.readmod('refocus.mod'))*sys.raster + arg.padtime; % us
dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + arg.padtime; % us
dur_crusher = length(toppe.readmod('crusher.mod'))*sys.raster + arg.padtime; % us

% Set cores file text
coresfiletext = [ ...
    "2 1 0" % tipdown, delay
    "2 2 0" % refocus, delay
    "2 3 0" % readout, delay
    "1 4" % crusher
    "1 0" % delay
    ];

% Set modules file text
modulesfiletext = [ ...
    sprintf("tipdown.mod %d 1 0 -1", dur_tipdown)
    sprintf("refocus.mod %d 1 0 -1", dur_refocus)
    sprintf("readout.mod %d 0 1 -1", dur_readout)
    sprintf("crusher.mod %d 0 0 -1", dur_crusher)
    ];

% Write out cores and modules files
psdutils.filegen.writemodulesfile(modulesfiletext);
psdutils.filegen.writecoresfile(coresfiletext);

% Set minimum TE
switch lower(arg.ro_mode)
    case 'fse'
        minTE = 1e-3*(dur_refocus + dur_readout);
    case 'gre'
        minTE = 0;
    otherwise
        error('invalid readout mode: %s\n', arg.ro_mode);
end
if strcmpi(arg.TE, 'min')
    arg.TE = minTE;
elseif arg.TE < minTE
    error('TE must be >= %fms\n', minTE);
end

% Set minimum TR
switch lower(arg.ro_mode)
    case 'fse'
        minTR = 0;
    case 'gre'
        minTR = 1e-3*(dur_tipdown + dur_readout + dur_crusher) + arg.TE; % ms;
    otherwise
        error('invalid readout mode: %s\n', arg.ro_mode);
end
if strcmpi(arg.TR, 'min')
    arg.TR = minTR;
elseif arg.TR < minTR
    error('TR must be >= %fms\n', minTR);
end

% Loop through frames, shots, echoes
toppe.write2loop('setup',sys,'version',6);
for framen = 1:arg.nframes
    for shotn = 1:arg.nshots

        if strcmpi(arg.ro_mode,'fse')
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
        end

        % Loop through echoes
        for echon = 1:arg.nechoes+arg.ndd

            if echon > arg.ndd % non-disdaq echoes
                dabmode = 'on';
                gamp = ones(3,1);
                viewn = (framen-1)*arg.nshots*arg.nechoes + (shotn-1)*arg.nechoes + echon - arg.ndd;
                rz = 2*pi * shotn / ((3-sqrt(5))/2);
                R = eul2rotm([rz,0,0],"ZYX");
            else
                dabmode = 'off';
                gamp = zeros(3,1);
                viewn = 1;
                R = eye(3);
            end

            if strcmpi(arg.ro_mode,'fse')

                % Write refocuser to loop
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

                % Write echo time deadtime to loop
                toppe.write2loop('delay', sys, ...
                    'textra', (arg.TE - minTE)/2, ...
                    'rotmat', R, ...
                    'core', 3);

            else

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
                    'core', 3);

                % Write repetition time deadtime to loop
                toppe.write2loop('delay', sys, ...
                    'textra', arg.TR - 1e-3*(dur_tipdown + dur_readout) - arg.TE, ...
                    'rotmat', R, ...
                    'core', 3);

                % Write crusher to loop
                toppe.write2loop('crusher.mod', sys, ...
                    'echo', 1, ...
                    'slice', 1, ...
                    'view', viewn, ...
                    'version', 6, ...
                    'core', 4);

            end


        end

        % Write repetition time deadtime to loop
        toppe.write2loop('delay', sys, ...
            'textra', arg.Tdelay, ...
            'core', 5);

    end
end

% Finish loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry','seqstamp.txt',sys);

end

