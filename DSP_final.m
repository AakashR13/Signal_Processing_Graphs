clear all;

p = 0.5;%freq lim
G = gsp_minnesota();
N = G.N;
x=rand(N,1);    
figure(1);
gsp_plot_signal(G,x);
title('original');
%display(x);
A = double(full(G.A));

%eigen
n = size(x,1);
[V,D] = eig(A);
[d,ind] = sort(diag(D),'descend');
D = D(ind,ind);
V = V(:,ind);
V = -1.*V;
%eigen

%scaling
scale = 1/D(1,1);
D = scale.*D;
d = scale.*d;
A = scale.*A;
%scaling
%{
%rounding
V = round(V,2);
d = round(d,2);
D = round(D,2);
%rounding
%}
%finding bandlimit/setting bandlimit
xf = V\x;
[B,I] = sort(xf);
cutoff = N*p;
I = I(1:cutoff);
xf(I) = 0;
%xf = round(xf,2);
k = nnz(xf);
x = V*xf;
%finding bandlimit/setting bandlimit

figure(2);
gsp_plot_signal(G,x);
title('bandlimited');

%constructing sampling operator
M = 1:k;%can be any k samples
K = 1:k;%must be first k columns
psi = zeros(k,n);
for x_ = 1:n
    for y_ = 1:k
        if x_ == M(y_)
            psi(y_,x_) = 1;
        end
    end
end
%constructing sampling operator

Vk(:,K) = V(:,K);
xm=psi*x;

%phi
phi = Vk*pinv(psi*Vk);
%phi = round(phi,2);
%phi

%recon
%x_ = phi*xm;
x_ = Vk*xf(1:k);
%x_ = round(x_,2);
figure(3);
gsp_plot_signal(G,x_);
title('reconstructed');
figure(4);
plot(1:N,abs(x-x_));
title('error(p=0.5)');
xlabel('coefficient');ylabel('magnitude');
MAE=mean(abs(x-x_));
display(MAE);
%recon
%1302.0708