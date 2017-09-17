clear all;
clc;
A=[0.678;1.1];
pi=[0.926, 0.074;0.023,0.977];
alpha=0.35;
delta=0.025;
beta=0.99;
sigma=2;
n=1000;
k_bar=A(2)*(alpha*A(2)/(1/beta-1+delta))^(alpha/(1-alpha))+(1-delta)*(alpha*A(2)/(1/beta-1+delta))^(1/(1-alpha));
K=linspace(0,k_bar,n+1);
K(1)=[];
V=zeros(n,2,2,2);
V(:,1,1,1)=K;
V(:,1,2,1)=K;
V(:,1,1,2)=K;
V(:,1,2,2)=K;
V(:,2,2,1)=10;
V(:,2,2,2)=20;
p=1;
t=1e-5;
l=1;
while [norm(V(:,:,end,1)-V(:,:,end-1,1));norm(V(:,:,end,2)-V(:,:,end-1,2))]>=t
    p=p+1;
    V(:,:,p,:)=V(:,:,p-1,:);
    m=zeros(n,n,2);
    for l = 1:2
        for i=1:n
            for j=1:n
            a=A(l)*K(i)^alpha+(1-delta)*K(i)-K(j);
                if a>0
                    m(i,j,l)=a^(1-sigma)/(1-sigma)+beta*(pi(l,1)*V(j,2,p-1,1)+pi(l,2)*V(j,2,p-1,2));
                else m(i,j,l)=-1000;
                end
            end
        end
    end
    [V(:,2,p,:),Ir]=max(m,[],2);

end
subplot(2,1,1);
plot(V(:,1,end,1),V(:,2,end,1));
title('Value Function, low');
xlabel('k');
ylabel('V(k)');
subplot(2,1,2);
plot(V(:,1,end,2),V(:,2,end,2));
title('Value Function, high');
xlabel('k');
ylabel('V(k)');
figure
k_prime=K(ind2sub(K,Ir));
plot(K,k_prime(:,1)); hold on;
plot(K,k_prime(:,2)); hold off;
legend('low state','high state','location','northwest');
title('Policy Function');