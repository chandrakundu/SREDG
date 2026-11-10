function SD=preeedm_masked_shortest_paths(D,W)
% Observed zero is a zero-weight edge; missing zero has no edge.
assert(isequal(D,D')&&isequal(W,W')&&all(diag(W)==0));
assert(all(isfinite(D(W)))&&all(D(W)>=0));
if all(D(W)>0)
 SD=graphallshortestpaths(sparse(D)); % no zero distances: same call as PREEEDM.m
else
 [i,j]=find(triu(W,1));weights=D(sub2ind(size(D),i,j));
 SD=distances(graph(i,j,weights,size(D,1)));
end
assert(all(isfinite(SD(:))),'The graph of observed pairs must be connected.');
end
