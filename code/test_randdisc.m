%% 
% test for discrete number generator

addpath('toolbox/');

normalize = @(p)p./repmat(sum(p,1),[size(p,1) 1]);


N = 12;
Q = 1000;
p = normalize(rand(N,1));

err = [];
Qlist = round(10.^linspace(2,8,15));
for Q=Qlist
    v = randdisc(p,Q);
    p1 = hist(v,1:N)/Q; p1 = p1(:);
    err(end+1) = norm(p1-p)/norm(p);
end
clf
plot(log10(err));