function Ph=gpfitEB(x);

% Empirical Bayesian Method : Zhang 2010

% initial definitions
x=sort(x(:));
n=length(x);
m=20+round(n^0.5);
xn=x(n);

% define prior for GPD parameters (text around Equation 8)
p=0.3:0.1:0.9;
indsp=round(n*(1-p)+0.5);
indsq=round(n*(1-p.^2)+0.5);
xp=x(indsp)';
xq=x(indsq)';
khat=log(xq./xp-1)./log(p);
sighat=khat.*xp./(1-p.^khat);
sighat(khat==0)=(-xp(khat==0)./log(p(khat==0)));
kstar=-1;
sigstar=1./(2*median(sighat));

% estimates of theta^star (Equation 6)
j=1:m;
theta=(n-1)/((n+1)*xn)-(sigstar./kstar).*(1-((j-0.5)/m).^kstar);

% evaluate log-likelihood function for each value of theta^star
L=0*theta;
for i=1:m
    L(i)=n*loglik(theta(i),x);
end

% evaluate weighting terms (text under Equation 6)
w=theta;
for i=1:m
    w(i)=1/sum(exp(L-L(i)));
end

% calculate estimate of theta (Equation 5)
theta=sum(theta.*w);

% calculate parameter estimates
xi=mean(log(1-theta.*x));
sigma=-xi/theta;

Ph=[xi sigma];

function l =loglik(b,x)
% log-likelihood function (text above Equation 1)
% note that factor of n is applied when function is called
k=-mean(log(1-b*x));
if b==0
    l=k-1-log(mean(x));
else
    l=k-1+log(b/k);
end

end

end
