seq.nframes = 1;
seq.nshots = 1;
seq.nechoes = 16;
seq.N = 64;
seq.fov = 24;
save seq

% Set up sequence
sys = toppe.systemspecs('maxSlew',20);
save sys
toppe.writeentryfile('toppeN.entry');

% Add a spiral GRE readout block with all default parameters
readout_block = psdutils.blocks.spiralgre(sys,...
    'tipdown_slthick', 0.5, ...
    'spiral_kxymax', seq.N/seq.fov/2, ...
    'spiral_kzmax', 0, ...
    'spiral_nshots', seq.nechoes);

% Write out the blocks
psdutils.filegen.writeblocks(readout_block);

% Loop through echo trains and play sequence
toppe.write2loop('setup',sys,'version',6);
viewn = 0;

rotmats = zeros(3,3,seq.nechoes);
for echon = 1:seq.nechoes
    rotmats(:,:,echon) = eul2rotm(pi*(3 - sqrt(5))*echon*[1,0,0]);
end

for framen = 1:seq.nframes

    % Play spiral GRE readout echo train
    viewn = readout_block.play('view0',viewn,'rotmatx',rotmats,'core0',0,'etl',seq.nechoes,'TE',0,'TR',20);

end

% Finish loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry','seqstamp.txt',sys);