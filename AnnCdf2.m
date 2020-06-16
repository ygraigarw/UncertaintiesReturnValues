function y=AnnCdf2(X,xi,sgm,Rat,RP);

m=size(xi,1);
n=size(X,1);
A=(1+(xi./sgm)*X');
B=(-1./xi)*ones(1,n);
IsOK=A>0;
y=zeros(m,n);
y(IsOK)=exp(-Rat*RP*( A(IsOK).^B(IsOK) ));
y(IsOK==0)=1;
y=mean(y)';

return;