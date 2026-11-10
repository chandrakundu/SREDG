function S = graphallshortestpaths(A)
% Same output as the Bioinformatics Toolbox function of this name, for symmetric A, without the toolbox.
assert(isequal(A,A'), 'Expected an undirected symmetric adjacency matrix.');
S = distances(graph(A,'upper'));
end
