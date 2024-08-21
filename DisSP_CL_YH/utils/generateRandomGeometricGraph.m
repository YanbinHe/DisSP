function [Graph] = generateRandomGeometricGraph(n, r, width)
    % n is the number of nodes
    % r is the connection radius
    % width is the width of the square
    % this code assumes that the area is square, otherwise you need to
    % reshape with both width and height
    
    % Generate random positions within a square width * width
    positions = width * rand(n, 2);
    % Generate Euclidean distance matrix
    Dist = squareform(pdist(positions));

    % Generate binary adjacency matrix
    Adj = zeros(n,n);
    Adj(Dist < r) = 1;
    Adj = Adj - diag(diag(Adj)); % no self-loop
    % Generate the graph using adjancecy matrix
    G = graph(Adj);

    % Generate node neighbours
    N = size(G.Nodes,1);
    node_neigh = cell(N,1);

    for n = 1:N
        node_neigh{n,1} = neighbors(G,n)';
    end

    % Customize the graph
    Graph = struct;
    Graph.G = G;
    Graph.posi = positions;
    Graph.edge_list = table2array(G.Edges(:,1));
    Graph.node_neigh = node_neigh;
    Graph.num_node = N;
    Graph.D = diag(degree(G));
    Graph.Adj = Adj;
end
