% Reads the output files from MSA

function [Y, Y2] = Read_KMCTTS_output(dir, N_record, n_params, n_specs)

cols = 1 + 2 * n_specs + n_params;                       % t, species, integral species, trajectory derivatives

fidread = fopen([dir '/MSA_output.bin'],'r');
Y = fread(fidread,[N_record,cols],'double');
fclose(fidread);

fidread = fopen([dir '/micro_derivs.bin'],'r');
Y2 = fread(fidread, n_params * n_specs * N_record, 'double');       % parameters * species * time points 
Y2 = reshape(Y2, [n_params, n_specs, N_record]);
fclose('all');

end