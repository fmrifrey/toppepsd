ntrains = 2;

% Set up sequence
sys = toppe.systemspecs('maxSlew',20);
save sys
toppe.writeentryfile('toppeN.entry');

% Add a spiral GRE readout block with all default parameters
readout_block = psdutils.blocks.spiralgre(sys);

% Make a deadtime block
deadtime_block.modulesfiletext = [];
deadtime_block.coresfiletext = "1 0";

% Write out the blocks
psdutils.filegen.writeblocks(readout_block, deadtime_block);

% Loop through echo trains and play sequence
toppe.write2loop('setup',sys,'version',6);
viewn = 0;
for trainn = 1:ntrains
    % Play spiral GRE readout
    viewn = readout_block.play('view0',viewn,'core0',0,'etl',20,'ndda',2,'TR',TRs);

    % Play deadtime
    toppe.write2loop('delay', sys, ...
        'textra', 500, ... % ms
        'core', 4);
end

% Finish loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry','seqstamp.txt',sys);