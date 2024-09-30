function [metric]=RanDropDeAsyADMM(x_ini,G,Graph_param,param)

% x_ini : initial measurement by each sensor in the order of sesnor index
% G     : graph
% param : parameters
% option: options

% metric: encapsulate all results

%% settings
x_optimal     = mean(x_ini);
n             = Graph_param.num_node;
node_nei      = Graph_param.node_neigh;
x             = x_ini;
rho           = param.rho;

% define and initialize auxiliary and dual variables
z             = zeros(n,n);
v             = zeros(n,n);
degree_G      = diag(Graph_param.D);

% track values
error         = [];
err           = 1;
num_trans     = 0;
num_attempted = 0;
transnum      = [];
error_to_real = [];

max_drop      = 20;
num_drop      = 1;
while err > param.epsilon

    if num_attempted > 0 && mod(num_attempted,0.5*1e5) == 0 && num_drop <= max_drop % randomly drop one node
        n = Graph_param.num_node;
        i_drop = unidrnd(n);
        G_temp = rmnode(G,i_drop);
        G = G_temp;

        % regenerate graph params
        Graph_param.posi(i_drop,:) = [];
        Graph_param.edge_list = table2array(G.Edges(:,1));

        N = size(G.Nodes,1);
        node_neigh = cell(N,1);

        for num = 1:N
            node_neigh{num,1} = neighbors(G,num)';
        end
        Graph_param.node_neigh = node_neigh;
        Graph_param.num_node = N;
        Graph_param.D = diag(degree(G));

        % update
        n             = Graph_param.num_node;
        node_nei      = Graph_param.node_neigh;
        degree_G      = diag(Graph_param.D);

        z(i_drop,:)  = [];
        z(:,i_drop)  = [];
        v(i_drop,:)  = [];
        v(:,i_drop)  = [];
        x_ini(i_drop)= [];
        x(i_drop)    = [];
        x_optimal = mean(x_ini);
        num_drop = num_drop + 1;
    end
    
    % randomly choose one node and then randomly choose a neighbour of this
    % node to establish the connection
    i = unidrnd(n);
    temp_list = node_nei(i,1);
    temp_list = temp_list{1};
    j = temp_list(unidrnd(length(temp_list)));
    
    % we do not update the neighbours of i and j but to collect values from
    % them
    node_i_neighbour = temp_list;
    temp_list = node_nei(j,1);
    temp_list = temp_list{1};
    node_j_neighbour = temp_list;
    
    % temporary copy of auxiliary and dual variable
    z_temp = z;
    v_temp = v;

    % attempted transmission
    num_attempted = num_attempted + 2;

    if rand > param.P_transfail
        %% update of Node i,j
        sum_i = 0;
        
        for neighbour = node_i_neighbour
            sum_i = sum_i + rho * z_temp(i,neighbour) - v_temp(i,neighbour);
        end
        x(i) = (x_ini(i) + sum_i) / (1 + rho * degree_G(i));

        sum_j = 0;
        for neighbour=node_j_neighbour
            sum_j = sum_j + rho * z_temp(j,neighbour) - v_temp(j,neighbour);
        end

        x(j) = (x_ini(j) + sum_j) / (1 + rho * degree_G(j));
        %% update of zij
        z(i,j) = 0.5 * (x(i) + x(j));
        z(j,i) = z(i,j);

        %% update of vi|j and vj|i

        v(i,j) = v_temp(i,j) + rho * (x(i) - z(i,j));
        v(j,i) = v_temp(j,i) + rho * (x(j) - z(j,i));

        num_trans = num_trans + 2;

    end

    transnum = [transnum num_attempted];
    err = norm(x - x_optimal)^2; 
    error=[error err];
    error_to_real = [error_to_real norm(x - 1)^2];
end

metric = cell(3,1);
metric{1} = transnum;
metric{2} = error;
metric{3} = error_to_real;
metric{4} = num_trans;

end