function out = blkdiagn(matrix, n)
out = matrix;
for ii = 2:n
    out = blkdiag(out,matrix); 
end
end