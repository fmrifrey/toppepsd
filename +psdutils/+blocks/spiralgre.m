classdef spiralgre

    properties
        sys
        dur_tipdown
        dur_readout
        dur_crusher
        tipdownfreq
        modulesfiletext
        coresfiletext
    end

    methods

        function obj = spiralgre(sys,varargin)
        % Input arguments:
        %   sys: toppe systemspecs object
        %
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

            % Set default arguments
            arg.tipdown_fa = 13;
            arg.tipdown_slthick = 20;
            arg.tipdown_tbw = 16;
            arg.tipdown_dur = 3200;
            arg.tipdown_nchop = [25 0];
            arg.tipdown_ftype = 'min';
            arg.spiral_F = [24, 0, 0];
            arg.spiral_kxymax = 1.33;
            arg.spiral_kzmax = 0;
            arg.spiral_dir = 'o';
            arg.spiral_nshots = 1;
            arg.spiral_nnav = 100;
            arg.crush_ncycles = 6;

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
        
            %% Generate crusher at the end of the readout
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

            % Set cores file text
            obj.coresfiletext = [ ...
                "2 1 0" % tipdown, delay
                "1 2" % readout
                "2 0 3" % delay, crush
                ];

            % Calculate all module durations
            obj.dur_tipdown = length(toppe.readmod('tipdown.mod'))*sys.raster + 20; % us
            obj.dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + 20; % us
            obj.dur_crusher = length(toppe.readmod('crusher.mod'))*sys.raster + 20; % us

            % Set modules file text
            obj.modulesfiletext = [ ...
                sprintf("tipdown.mod %d 1 0 -1", obj.dur_tipdown)
                sprintf("readout.mod %d 0 1 -1", obj.dur_readout)
                sprintf("crusher.mod %d 0 0 -1", obj.dur_crusher)
                ];

            obj.sys = sys;
            
        end

        function viewn = play(obj,varargin)
        % Input arguments:
        %   view0: storage view index of first readout in train
        %   core0: core index of first core in block
        %   etl: number of TRs in echo train (including disdaqs),
        %       default = 1
        %   ndd: number of disdaq TRs, default = 0
        %   TE: echo time value or an array of TE's for each echo (ms),
        %       default = 'min'
        %   TR: repitition time value or an array of TR's for each echo
        %       (ms), default = 'min'
        %   rotmatx: 3D transformation matrix or 3x3x3xetl list of 3D
        %       transformation matrices for each echo, default = eye(3)

            arg.view0 = 0;
            arg.core0 = 0;
            arg.etl = 1;
            arg.ndd = 0;
            arg.TE = 'min';
            arg.TR = 'min';
            arg.rotmatx = eye(3);
            arg = toppe.utils.vararg_pair(arg, varargin);

            minTR = 1e-3*(obj.dur_tipdown + obj.dur_readout + ...
                obj.dur_crusher + max(arg.TE)); % ms
            minTE = 0;

            % Set minimum TE and TR values
            if strcmpi(arg.TE, 'min')
                arg.TE = minTE;
            end

            if strcmpi(arg.TR, 'min')
                arg.TR = minTR;
            end

            % Translate single parameters into a list
            if length(arg.TE) == 1
                arg.TE = repmat(arg.TE,arg.etl,1);
            end
    
            if length(arg.TR) == 1
                arg.TR = repmat(arg.TR,arg.etl,1);
            end
    
            if size(arg.rotmatx,4) == 1
                arg.rotmatx = repmat(arg.rotmatx,[1,1,1,arg.etl]);
            end

            % Loop through echos in echo train
            for i = 1:arg.etl
                if i > arg.ndd % non-disdaq echoes
                    viewn = arg.view0 + i;
                    dabmode = 'off';
                else
                    viewn = 1;
                    dabmode = 'on';
                end

                % Write tipdown to loop
                toppe.write2loop('tipdown.mod', obj.sys, ...
                    'RFoffset', obj.tipdownfreq, ...
                    'echo', 1, ...
                    'slice', 1, ...
                    'view', viewn, ...
                    'RFspoil', 1, ...
                    'version', 6, ...
                    'trigout', 0, ...
                    'core', arg.core0 + 1);

                % Write echo time deadtime to loop
                toppe.write2loop('delay', obj.sys, ...
                    'textra', arg.TE(i), ... % ms
                    'core', arg.core0 + 1);

                % Write readout to loop
                toppe.write2loop('readout.mod', obj.sys, ...
                    'echo', 1, ...
                    'slice', 1, ...
                    'view', viewn, ...
                    'RFspoil', 1, ...
                    'version', 6, ...
                    'trigout', 0, ...
                    'dabmode', dabmode, ...
                    'rotmat', arg.rotmatx(:,:,:,i), ...
                    'core', arg.core0 + 2);
                 
                % Write repetition time deadtime to loop
                toppe.write2loop('delay', obj.sys, ...
                    'textra', arg.TR(i) - minTR, ... % ms
                    'core', arg.core0 + 3);

                % Crusher
                toppe.write2loop('crusher.mod', obj.sys, ...
                    'echo', 1, ...
                    'slice', 1, ...
                    'view', viewn, ...
                    'version', 6, ...
                    'core', arg.core0 + 3);

            end

        end
    end

end
