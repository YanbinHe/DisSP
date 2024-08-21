clc
clear

addpath ./utils

% fix random seed to reproduce all the results
rng(exp(pi))

% flag
option.show_graph = 0;
option.test_admm = 1;
option.test_rangos = 1;
fontsizemean = 20;
%% generate the graph
Width   = 100;
Height  = 100;
radius  = 30; % radius
radiusn = radius/Width; % normalized to [0,1]
n       = 100:1000;
tmp     = sqrt(2 * log(n)./n);
n       = n(find(tmp < radiusn, 1)); % the number of nodes
prob    = 1 - 1/n^2; % the probability of connection

G = generateRandomGeometricGraph(n, radius, Width);

if option.show_graph
    figure(1)
    plot(G.G,'XData',G.posi(:,1),'YData',G.posi(:,2));
    xlim([0,100])
    ylim([0,100])
end

%% signal generation
noise = sqrt(0.01) * randn(n,1);
x_ini = ones(n,1) + noise;


%% randomized gossip
if option.test_rangos
    % proceed with optimized P
    Popt = OptP(n,G.node_neigh);
    % proceed with random P
    Pran = G.Adj .* rand(n);
    Pran = Pran ./ sum(Pran,2);
    % proceed with natural averaging P
    Pavg = G.Adj;
    Pavg = Pavg ./ sum(Pavg,2);

    RandomizedGossip(x_ini,G,option,transFail);
end


%% decentralized asychronized admm
if option.test_admm
    param.rho         = 0.4;
    param.epsilon     = 1e-12;

    %% study transmission failure
    for p = 0:0.1:0.9
        param.P_transfail = p;
        [result] = DeAsyADMM(x_ini,G,param);
        fname = ['./results/de_asy_admm_p_', num2str(p),'.mat'];
        save(fname,'result')
    end
    %% plot comparison
    error_collect = cell(10,2); % transmissions;error
    for i = 1:10
        p = (i-1) / 10;
        fname = ['./results/de_asy_admm_p_', num2str(p),'.mat'];
        load(fname,'result')

        error_collect{i,1} = result{1};
        error_collect{i,2} = result{2};
    end

    figure(2)
    for i = 1:10
        caption = ['p = ', num2str((i-1) / 10)];
        plot(error_collect{i,1},error_collect{i,2},'Display',caption)
        hold on
    end


    grid on
    set(gca, 'yscale', 'log');
    legend('boxoff')
    xlabel('Number of Transmission')
    ylabel('Error')
    ylim([1e-12,1e3])
end