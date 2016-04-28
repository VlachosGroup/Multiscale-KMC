function parts = strsplit(stringin,delimiter)

strln = length(stringin);
if (strln < 3)
    error('Input string consists of less than three characters.');
end

nwords = 0;

if (stringin(1) ~= delimiter) % record begins with a word
	nwords = nwords + 1;
    iw(1) = 1;
	inword = true;
end

for i = 2:strln

	if (stringin(i) == delimiter && ...
		stringin(i-1) ~= delimiter) % end of word
		inword = false;
		iw(2*nwords) = i-1;
    end

	if (stringin(i) ~= delimiter && ...
		stringin(i-1) == delimiter) % beginning of word
		inword = true;
		nwords = nwords + 1;
		iw(2*nwords-1) = i;
    end

end

if (inword)
	iw(2*nwords) = strln;
end

for i = 1:nwords
    parts{i} = stringin(iw(2*i-1):iw(2*i));
end

end
