function P=AnnExcPrb(R,xi,sgm,NO);

P=1-exp(-NO.*(1-gpcdf(R,xi,sgm,0)));

return;