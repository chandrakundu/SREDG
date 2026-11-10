function x=dcroot_zeroaware(w,a,Ind,d2,l,u)
% Original dcroot for positive distances; observed zero distances use preeedm_zero_scalar_min.
positive=Ind(d2(Ind)>0);x=dcroot(w,a,positive,d2,l,u);
zero=Ind(d2(Ind)==0);
x(zero)=preeedm_zero_scalar_min(w(zero),a(zero),l(zero),u(zero));
end
