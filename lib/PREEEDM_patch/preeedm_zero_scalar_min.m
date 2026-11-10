function x=preeedm_zero_scalar_min(w,a,l,u)
% Global minimizer of .5*(x-w)^2+a*sqrt(x), l<=x<=u, 0<=l<=u.
% t=sqrt(x): derivative 2*t^3-2*w*t+a. Compare endpoints and the
% larger nonnegative stationary root (the only possible interior minimum).
assert(isequal(size(w),size(a),size(l),size(u)));
assert(all(isfinite(w(:)))&&all(isfinite(a(:)))&&all(isfinite(l(:)))&&all(isfinite(u(:))));
assert(all(a(:)>=0)&&all(l(:)>=0)&&all(u(:)>=l(:)));
x=l; best=.5*(l-w).^2+a.*sqrt(l);
fu=.5*(u-w).^2+a.*sqrt(u); take=fu<best;x(take)=u(take);best(take)=fu(take);
no_weight=a==0;x(no_weight)=min(u(no_weight),max(l(no_weight),w(no_weight)));
indices=find(w>0&a>0&u>l);v=w(indices)/3;
c=a(indices)./(4*v.*sqrt(v));valid=c<=1;
indices=indices(valid);v=v(valid);c=c(valid);
candidate=4*v.*cos(acos(-min(1,max(0,c)))/3).^2;
inside=candidate>=l(indices)&candidate<=u(indices);indices=indices(inside);candidate=candidate(inside);
fc=.5*(candidate-w(indices)).^2+a(indices).*sqrt(candidate);
take=fc<best(indices);x(indices(take))=candidate(take);
end
