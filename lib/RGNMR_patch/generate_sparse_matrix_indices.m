function colind_A = generate_sparse_matrix_indices(omega,r,n2)
% Plain MATLAB version of c-utils/generate_sparse_matrix_indices.c (same output, no MEX file).
cols=[r*(omega(:,2)-1)+(1:r),r*n2+r*(omega(:,1)-1)+(1:r)];
colind_A=reshape(cols',1,[]);
end
