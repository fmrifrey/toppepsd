function Wi = pipe_menon_dcf(G,itrmax)

    % Set default for itrmax
    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end
    
    % If G is a Gmri object, use its Gnufft object
    if isfield(G.arg,'Gnufft')
        G = G.Gnufft;
    end
    
    % Initialize weights to 1 (psf)
    Wi = ones(size(G,1),1);
    
    % Loop through iterations
    for itr = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( G.arg.st.interp_table(G.arg.st, ...
            G.arg.st.interp_table_adj(G.arg.st, Wi) ) );
        Wi = Wi ./ d;
        
    end
    
    % Normalize weights
    Wi = Wi / sum(abs(Wi));
    
end
