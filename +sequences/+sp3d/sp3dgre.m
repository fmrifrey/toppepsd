function sys = sp3dgre(varargin)
% Input arguments:
%   path2experiment: directory for experiment (where to write files)
%
%   D: field of view (cm)
%   N: image matrix size
%   nechoes: number of echoes in each echo train
%   nshots: number of shots per frame
%   nframes: number of frames
%   TE: echo time (ms)
%   TR: repitition time (tipdown to tipdown, ms)
%   frameTR: frame repitition time (echo train to echo train, ms)
%
%   spvd0: spiral center oversampling factor
%   spvd1: spiral edge oversampling factor
%   spdir: spiral direction ('i' = spiral in, 'o' = spiral out, 'io' =
%       spiral in-out, 'oi' = spiral out-in)
%   rtype: 3d transformation type ('sos' = stack of spirals, 'r1ax' =
%       single axis golden angle rotations, 'r2ax' = dual axis golden angle
%       rotations, 'fromlist' = read from rlist)
%   rlist: 3d transformation parameters (in format
%       [thetax(:), thetay(:), thetaz(:), dz(:)])
%   nnav: number of navigator samples
%
%   rf_excite_slthick: rf excitation slab thickness (cm)
%   rf_excite_dz: rf excitation slab z offset (cm)
%   rf_excite_fa: rf excitation flip angle (deg)
%   rf_excite_tbw: rf excitation time-bandwidth product
%   rf_excite_dur: rf excitation pulse duration
%
%   Gmax: maximum allowed gradient amplitude (G/cm)
%   Smax: maximum allowed slew rate (G/cm/s)
%

%% Set default parameters

% Image parameters
arg.D = 24;
arg.N = 64;

% Acquisition parameters
arg.nechoes = 16;
arg.nshots = 1;
arg.nframes = 4;
arg.TE = 5;
arg.TR = 100;
arg.frameTR = 5000;

% Spiral parameters
arg.spvd0 = 1;
arg.spvd1 = 1;
arg.spdir = 'o';
arg.rtype = 'r2ax';
arg.rlist = [];
arg.nnav = 20;

% RF parameters
arg.rf_excite_slthick = 20;
arg.rf_excite_dz = 0;
arg.rf_excite_fa = 13;
arg.rf_excite_tbw = 12;
arg.rf_excite_dur = 3200;

% System parameters
arg.Gmax = 4;
arg.Smax = 17000;

%% Save system
sys = arg.sys;
save sys

% Store arguments into seq
seq = toppe.utils.vararg_pair(arg, varargin);
save seq

%% Create the module and core files
toppe.writeentryfile;

% Calculate all module durations
dur_excite = length(toppe.readmod('excite.mod'))*sys.raster + 100; % us
dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + 100; % us
dur_crusher = length(toppe.readmod('crusher.mod'))*sys.raster + 100; % us

% Calculate deadtimes
deadtime_TR = seq.TR - seq.TE - 1e-3*(dur_excite + dur_readout + dur_crusher);
deadtime_frameTR = seq.frameTR - seq.nechoes*seq.TR;

% Write the modules.txt file
modFileText = ['' ...
    'Total number of unique cores\n' ...
    '3\n' ...
    'fname  duration(us)    hasRF?  hasDAQ? hastrigout?\n' ...
    sprintf('excite.mod %d 1 0 -1\n', dur_excite) ...
    sprintf('readout.mod %d 0 1 -1\n', dur_readout) ...
    sprintf('crusher.mod %d 0 0 -1\n', dur_crusher)];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);
psdutils.filegen.writemodulesfile(mods);

% Write the cores.txt file
coreFileText = {[1, 0] % tip, delay
    2 % readout
    [0, 3] % delay, crusher
    0}; % deadtime
toppe.writecoresfile(coreFileText);

%% Loop through the sequence
toppe.write2loop('setup',sys,'version',6);

viewn = 1;
for framen = 1:seq.nframes
    for shotn = 1:seq.nshots

        for echon = 1:seq.nechoes

            % Generate rotation matrix
            R = genview(shotn, seq.nshots, echon, seq.nechoes, seq.rtype, seq.rlist);

            % Excitation
            toppe.write2loop('excite.mod',sys, ...
                'RFoffset', freq_excite, ...
                'echo', 1, ...
                'slice', viewn, ...
                'view', 1, ...
                'RFspoil', 1, ...
                'version', 6, ...
                'trigout', 0, ...
                'core', 1);

            % Echo time delay
            toppe.write2loop('delay', sys, ...
                'textra', seq.TE, ... % ms
                'core', 1)    

            % Readout
            toppe.write2loop('readout.mod', sys, ...
                'echo', 1, ...
                'slice', viewn, ...
                'view', 1, ...
                'RFspoil', 1, ...
                'version', 6, ...
                'trigout', 0, ...
                'dabmode', 'on', ...
                'rotmat', R, ...
                'core', 2);

            % Repetition time delay
            toppe.write2loop('delay', sys, ...
                'textra', deadtime_TR, ... % ms
                'core', 3);

            % Crusher
            toppe.write2loop('crusher.mod', sys, ...
                'echo', 1, ...
                'slice', viewn, ...
                'view', 1, ...
                'version', 6, ...
                'core', 3);

            viewn = viewn + 1;
        end

        % Frame deadtime
        toppe.write2loop('delay', sys, ...
            'textra', deadtime_frameTR, ...
            'core', 4);

    end
end

% Write out the loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sys);
if ~isfile('./toppeN.entry')
    copyfile(which('toppeN.entry'), './');
end

cd(wd);

end

