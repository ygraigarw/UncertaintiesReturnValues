function q=RtrVal(xi,sgm,NO);

q(xi==0)=sgm(xi==0)*log(NO);
q(xi~=0)=(sgm(xi~=0)./xi(xi~=0)).*(NO.^xi(xi~=0)-1);

return;