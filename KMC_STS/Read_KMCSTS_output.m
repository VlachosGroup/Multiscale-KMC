% Reads the output file from STS KMC simualtion

% Still need to make the script which does the averaging to compute
% sensitivities

function [t, N, N_int, W] = Read_KMCSTS_output(dir, N_record, n_params, n_specs)

cols = 1 + 2 * n_specs + n_params;                       % t, species, integral species, trajectory derivatives

fidread = fopen([dir '/KMC_STS_output.bin'],'r');
Y = fread(fidread,[N_record,cols],'double');
fclose(fidread);
Y(1,:) = [];        % Cut off t = 0 data

t = Y(:,1);
N = Y(:,2:1+n_specs);
N_int = Y(:,2+n_specs:1+2*n_specs);
W = Y(:,2+2*n_specs:end);

end