function covar = cov_calc( data1, data2 )

matr = cov(data1,data2);
covar = matr(1,2);

end

