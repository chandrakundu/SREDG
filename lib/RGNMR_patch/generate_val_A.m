function val_A = generate_val_A(U,V,omega,~)
% Plain MATLAB version of c-utils/generate_val_A.c (same output, no MEX file).
values=[U(omega(:,1),:),V(omega(:,2),:)];
val_A=reshape(values',1,[]);
end
