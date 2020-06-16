function Ph=gpfitPWM(x)

% PWM estimates of parameters
a0 = mean(x);
x = sort(x);
x = x(:);
n = length(x);
nn = n+0.35-(1:n);
a1 = sum(nn'.*x)/n/n;

sig = 2*a0*a1/(a0-2*a1);
xi = 2-a0/(a0-2*a1);
if xi<0 && x(n)>-sig/xi
    xi=-sig/x(n);
end

Ph=[xi sig];

return;