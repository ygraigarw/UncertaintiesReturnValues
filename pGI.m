function pGI(Enw0,Fmt);

if ~strcmp(get(gcf,'WindowStyle'),'docked')
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
end;

set(gcf,'Color',[1 1 1]);

Stn='jpeg';

if nargin==0;
    Enw0 = sprintf('GI_%s', datestr(now,30));
end;

if nargin<=1;
    Fmt=0;
end

if nargin==2;
    if Fmt==1;
        Stn='tiff';
    elseif Fmt==2;
        Stn='png';
    elseif Fmt==0;
        Stn='jpeg';
    elseif Fmt==3;
        Stn='jpeg';
    end;
end

tStr=sprintf('Enw=''%s.%s'';',Enw0,Stn);eval(tStr);

style = hgexport('factorystyle');
hgexport(gcf, Enw, style, 'Format', Stn);

return;