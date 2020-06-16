function Ph=gpfitMOM(x)

xn=max(x);
m=mean(x);
v=var(x);

% MOM estimates of parameters
sig=m*(1-0.5*(1-m^2/v));
xi=0.5*(1-m^2/v);
if xi<0 && xn>-sig/xi
    xi=-sig/xn;
end

Ph=[xi sig];

return;
