clc
clear

addpath ./utils

% fix random seed to reproduce all the results
rng(1)

% flag
option.show_graph = 1;
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

[G,Graph_param] = generateRandomGeometricGraph(n, radius, Width);

if option.show_graph
    figure(1)
    plot(G,'XData',Graph_param.posi(:,1),'YData',Graph_param.posi(:,2));
    xlim([0,100])
    ylim([0,100])
end

%% signal generation
noise = sqrt(0.01) * randn(n,1);
x_ini = ones(n,1) + noise;


%% randomized gossip
% if option.test_rangos
%     % proceed with optimized P
%     Popt = OptP(n,Graph_param.node_neigh);
%     % proceed with random P
%     Pran = Graph_param.Adj .* rand(n);
%     Pran = Pran ./ sum(Pran,2);
%     % proceed with natural averaging P
%     Pavg = Graph_param.Adj;
%     Pavg = Pavg ./ sum(Pavg,2);
% 
%     RandomizedGossip(x_ini,G,option,transFail);
% end


%% decentralized asychronized admm
if option.test_admm
    param.rho         = 0.4;
    param.epsilon     = 1e-12;

    %% study transmission failure
    for p = 0:0.1:0.9
        param.P_transfail = p;
        [result] = DeAsyADMM(x_ini,Graph_param,param);
        fname = ['./results/de_asy_admm_p_', num2str(p),'.mat'];
        save(fname,'result')
    end

    %% study node drop
    for p = 0:0.1:0.9
        param.P_transfail = p;
        [result] = RanDropDeAsyADMM(x_ini,G,Graph_param,param);
        fname = ['./results/drop_de_asy_admm_p_', num2str(p),'.mat'];
        save(fname,'result')
    end
    %% load data for comparison
    error_collect = cell(10,6); % transmissions;error
    for i = 1:10
        p = (i-1) / 10;
        fname = ['./results/de_asy_admm_p_', num2str(p),'.mat'];
        load(fname,'result')

        error_collect{i,1} = result{1};
        error_collect{i,2} = result{2};
        error_collect{i,3} = result{3};

        fname = ['./results/drop_de_asy_admm_p_', num2str(p),'.mat'];
        load(fname,'result')

        error_collect{i,4} = result{1};
        error_collect{i,5} = result{2};
        error_collect{i,6} = result{3};
    end

    %% plot data for comparison
    figure(2)
    for i = 1:10
        caption = ['p = ', num2str((i-1) / 10)];
        plot(error_collect{i,1},error_collect{i,2},'Display',caption)
        hold on
        % plot(error_collect{i,1},error_collect{i,3},'Display',caption)
        % hold on
        % plot(error_collect{i,4},error_collect{i,5},':','Display',caption)
        % hold on
        % plot(error_collect{i,4},error_collect{i,6},':','Display',caption)
        % hold on
    end

    grid on
    set(gca, 'yscale', 'log');
    legend('boxoff')
    set(0,'DefaultLineLineWidth',3)
    set(0,'DefaultLineMarkerSize',14)
    set(0,'DefaultAxesFontWeight','bold')
    xlabel('Number of Transmission')
    ylabel('Error')
    ylim([1e-12,1e3])

    %%
    figure(3)
    for i = 3
        caption = ['p = ', num2str((i-1) / 10)];
        plot(error_collect{i,1},error_collect{i,2},'b-','Display',caption)
        hold on
        plot(error_collect{i,1},error_collect{i,3},':b','Display',caption)
        hold on
        plot(error_collect{i,4},error_collect{i,5},'r-','Display',caption)
        hold on
        plot(error_collect{i,4},error_collect{i,6},':r','Display',caption)
        hold on
    end


    grid on
    set(gca, 'yscale', 'log');
    legend('boxoff')
    set(0,'DefaultLineLineWidth',3)
    set(0,'DefaultLineMarkerSize',14)
    set(0,'DefaultAxesFontWeight','bold')
    xlabel('Number of Transmission')
    ylabel('Error')
    ylim([1e-12,1e3])
end