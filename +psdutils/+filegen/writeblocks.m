function writeblocks(varargin)

modulesfiletext = [];
coresfiletext = [];
nmods = length(varargin);
for imod = 1:nmods
    block = varargin{imod};
    modulesfiletext = [modulesfiletext; block.modulesfiletext];
    coresfiletext = [coresfiletext; block.coresfiletext];
end

psdutils.filegen.writemodulesfile(modulesfiletext);
psdutils.filegen.writecoresfile(coresfiletext);