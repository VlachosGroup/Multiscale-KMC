% Reads the input file and passes system information to the main code
% Make sure to implement error handling in here

% Incomplete, not yet functional

function system_info = InputRead()

fid = fopen('MSA_input.txt');

fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
line = fgetl(fid);

species = strsplit(line,',')

system_info = [];

fclose('all');

end