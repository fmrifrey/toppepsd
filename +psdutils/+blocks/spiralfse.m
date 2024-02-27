classdef spiralfse

    properties
        sys
        dur_tipdown
        dur_refocus
        dur_readout
        dur_crusher
        tipdownfreq
        refocusfreq
        modulesfiletext
        coresfiletext
    end

    methods

        function obj = spiralfse(sys,varargin)
        % Input arguments:
        %   sys: toppe systemspecs object
        %
        %   tipdown_fa: tipdown flip angle (deg), default = 90
        %   tipdown_slthick: tipdown slab thickness (cm), default = 20
        %   tipdown_tbw: tipdown time bandwidth product, default = 4
        %   tipdown_dur: tipdown duration (us), default = 3200
        %   tipdown_nchop: number of dead samples before and after
        %       tipdown, default = [25, 0]
        %   tipdown_ftype: slr pulse type, default = 'min'
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
        %   refocus_fa: refocuser flip angle (deg), default = 180
        %   refocus_slthick: refocuser slab thickness (cm), default = 20
        %   refocus_tbw: refocuser time bandwidth product, default = 2
        %   refocus_dur: refocuser duration (us), default = 3200
        %   refocus_nchop: number of dead samples before and after
        %       refocuser, default = [25, 25]
        %   refocus_ftype: slr pulse type, default = 'ls'
        %   
        %   crush_ncycles: number of cycles of phase across slab
        %

            % Set default arguments
            arg.tipdown_fa = 90;
            arg.tipdown_slthick = 20;
            arg.tipdown_tbw = 4;
            arg.tipdown_dur = 3200;
            arg.tipdown_nchop = [25 0];
            arg.tipdown_ftype = 'min';
            arg.spiral_F = [24, 0, 0];
            arg.spiral_kxymax = 1.33;
            arg.spiral_kzmax = 0;
            arg.spiral_dir = 'io';
            arg.spiral_nshots = 1;
            arg.spiral_nnav = 100;
            arg.refocus_fa = 180;
            arg.refocus_slthick = 20;
            arg.refocus_tbw = 2;
            arg.refocus_dur = 3200;
            arg.refocus_nchop = [25, 25];
            arg.crush_ncycles = 16;
            arg.padtime = 1000;

            % Parse arguments
            arg = toppe.utils.vararg_pair(arg, varargin);

            % Set tipdown type
            if arg.tipdown_fa < 30
                tipdown_type = 'st';
            else
                tipdown_type = 'ex';
            end

            % Create SLR pulse
            [tipdown_rf, tipdown_gz, obj.tipdownfreq] = toppe.utils.rf.makeslr( ...
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
      
            % Create SLR pulse
            [refocus_rf, refocus_gz, obj.refocusfreq] = toppe.utils.rf.makeslr( ...
                arg.refocus_fa, ... % deg
                arg.refocus_slthick, ... % cm
                arg.refocus_tbw, ...
                arg.refocus_dur*1e-3, ... % ms
                arg.crush_ncycles, ...
                sys, ...
                'type', 'se', ...
                'writeModFile', false ...
                );
        
            % Write tipdown out to mod file
            toppe.writemod(sys, ...
                'rf', refocus_rf, ...
                'gz', refocus_gz, ...
                'nChop', arg.refocus_nchop, ...
                'ofname', 'refocus.mod');
      
            % Create initial spiral trajectory
            g_sp0 = psdutils.spiral.spgrad(sys, arg.spiral_dir, ...
                arg.spiral_F, arg.spiral_kxymax, arg.spiral_kzmax, ...
                arg.spiral_nshots, arg.spiral_nnav);
            
            % Write readout out to mod file
            toppe.writemod(sys, ...
                'gx', g_sp0(1,:)', ...
                'gy', g_sp0(2,:)', ...
                'gz', g_sp0(3,:)', ...
                'nChop', [10,10], ... % needs enough time for adc to catch up
                'ofname', 'readout.mod');

            % Calculate all module durations
            obj.dur_tipdown = length(toppe.readmod('tipdown.mod'))*sys.raster + arg.padtime; % us
            obj.dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + arg.padtime; % us
            obj.dur_refocus = length(toppe.readmod('refocus.mod'))*sys.raster + arg.padtime; % us

            % Set cores file text
            obj.coresfiletext = [ ...
                "2 1 0" % tipdown
                "2 2 0" % refocus
                "2 3 0" % readout
                ];

            % Set modules file text
            obj.modulesfiletext = [ ...
                sprintf("tipdown.mod %d 1 0 -1", obj.dur_tipdown)
                sprintf("readout.mod %d 0 1 -1", obj.dur_readout)
                sprintf("refocus.mod %d 1 0 -1", obj.dur_refocus)
                ];

            obj.sys = sys;
            
        end

        function viewn = play(obj,varargin)
        % Input arguments:
        %   view0: storage view index of first readout in train
        %   core0: core index of first core in block
        %   etl: number of TEs in echo train (including disdaqs),
        %       default = 1
        %   TE: echo time (ms), default = 'min'
        %   rotmatx: 3D transformation matrix or 3x3x3xetl list of 3D
        %       transformation matrices for each echo, default = eye(3)
        %   dabmode: turn dab 'on' or 'off'

            arg.view0 = 0;
            arg.core0 = 0;
            arg.etl = 1;
            arg.TE = 'min';
            arg.rotmatx = eye(3);
            arg.dabmode = 'on';
            arg = toppe.utils.vararg_pair(arg, varargin);

            % Set minimum TE and TR values
            minTE = 1e-3*(obj.dur_readout + obj.dur_refocus);
            if strcmpi(arg.TE, 'min')
                arg.TE = minTE;
            end

            % Translate single parameters into a list
            if size(arg.rotmatx,4) == 1
                arg.rotmatx = repmat(arg.rotmatx,[1,1,1,arg.etl]);
            end

            viewn = arg.view0 + 1;

            % Write tipdown to loop
            toppe.write2loop('tipdown.mod', obj.sys, ...
                'RFoffset', obj.tipdownfreq, ...
                'echo', 1, ...
                'slice', 1, ...
                'view', viewn, ...
                'RFspoil', 0, ...
                'version', 6, ...
                'trigout', 0, ...
                'core', arg.core0 + 1);

            % Write time delay to loop
            toppe.write2loop('delay', obj.sys, ...
                'textra', (arg.TE - obj.dur_refocus*1e-3)/2, ...
                'core', arg.core0 + 1);

            % Loop through echos in echo train
            for i = 1:arg.etl

                viewn = arg.view0 + i;

                % Write refocuser to loop
                toppe.write2loop('refocus.mod', obj.sys, ...
                    'RFoffset', obj.refocusfreq, ...
                    'echo', 1, ...
                    'slice', 1, ...
                    'view', viewn, ...
                    'RFspoil', 0, ...
                    'version', 6, ...
                    'trigout', 0, ...
                    'core', arg.core0 + 2);

                % Write time delay to loop
                toppe.write2loop('delay', obj.sys, ...
                    'textra', (arg.TE - minTE)/2, ...
                    'core', arg.core0 + 2);

                % Write readout to loop
                toppe.write2loop('readout.mod', obj.sys, ...
                    'echo', 1, ...
                    'slice', 1, ...
                    'view', viewn, ...
                    'RFspoil', 0, ...
                    'version', 6, ...
                    'trigout', 0, ...
                    'dabmode', arg.dabmode, ...
                    'rotmat', arg.rotmatx(:,:,i), ...
                    'core', arg.core0 + 3);

                % Write time delay to loop
                toppe.write2loop('delay', obj.sys, ...
                    'textra', (arg.TE - minTE)/2, ...
                    'core', arg.core0 + 3);

            end

        end
    end

end
