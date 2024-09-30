function [metric]=RandomizedGossip(x_ini,G,P,param)
n= G.num_node;
N      = G.node_neigh;
Z=x_ini;
error=[];
err=Inf;
avg=mean(x_ini);
transnum=[];
num=0;
while err>param.epsilon
    i=unidrnd(n);
    tmp=unifrnd(0,1);
    p=P(i,:);
    p=p(p~=0);
    p=cumsum(p);
    j_posi=find(p>tmp, 1 );
    j=N{i}(j_posi);
    num=num+2;
    if rand> param.P_transfail
%         num=num+2;
        avgz=(Z(i)+Z(j))/2;
        Z(i)=avgz;
        Z(j)=avgz;
    end
    err=sum((Z-avg).^2);
    transnum=[transnum num];
    error=[error err];
end
metric = cell(2,1);
metric{1} = transnum;
metric{2} = error;

end