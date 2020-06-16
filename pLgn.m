function pLgn(Txt,Lct,FntSiz,LinWdt,LinStl,LinClr,Gap);

nT=size(Txt,1);

if 1;
    limits = objbounds(findall(gca));
    xl=limits(1); xh=limits(2); yl=limits(3); yh=limits(4); zl=limits(5); zh=limits(6);
    rx=xh-xl; ry=yh-yl;
else; %alternative
    t=get(gca,'xlim');xl=t(1);xh=t(2);
    t=get(gca,'ylim');yl=t(1);yh=t(2);
    rx=xh-xl; ry=yh-yl;
    xl=xl+rx*0.05;
    xh=xh-rx*0.05;
    yl=yl+ry*0.05;
    yh=yh-ry*0.05;
    rx=xh-xl; ry=yh-yl;
end;


if nargin==4;
    LinStl=cell(8,1);
    LinStl{1}='-'; %solid
    LinStl{2}='--';%dashed
    LinStl{3}=':';%dotted
    LinStl{4}='-.';%dash-dot
    LinStl{5}='-'; %solid
    LinStl{6}='--';%dashed
    LinStl{7}=':';%dotted
    LinStl{8}='-.';%dash-dot
end;

if nargin<6;
    LinClr=cell(8,1);
    LinClr{1}=zeros(3,1);
    LinClr{2}=zeros(3,1);
    LinClr{3}=zeros(3,1);
    LinClr{4}=zeros(3,1);
    LinClr{5}=ones(3,1)*0.7;
    LinClr{6}=ones(3,1)*0.7;
    LinClr{7}=ones(3,1)*0.7;
    LinClr{8}=ones(3,1)*0.7;
end;

if nargin<7;
    Gap=0.05;
end;

switch Lct;
    case 'NW';
        for j=1:nT;
            if j<=4;
                plot([xl xl+0.2*rx],[yh yh]-Gap*(j-0.5)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',LinClr{j});
                text(xl+0.25*rx,yh-Gap*(j-0.5)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','left','Interpreter','latex');
            elseif j<=8;
                plot([xl xl+0.2*rx],[yh yh]-Gap*(j-0.5)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',ones(3,1)*0.5,'color',LinClr{j});
                text(xl+0.25*rx,yh-Gap*(j-0.5)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','left','Interpreter','latex');
            end;                
        end;
    case 'NE';
        for j=1:nT;
           if j<=4;
               plot([xh-0.2*rx xh],[yh yh]-Gap*(j-0.5)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',LinClr{j});
               text(xh-0.25*rx,yh-Gap*(j-0.5)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','right','Interpreter','latex');
           elseif j<=8;
               plot([xh-0.2*rx xh],[yh yh]-Gap*(j-0.5)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',ones(3,1)*0.5,'color',LinClr{j});
               text(xh-0.25*rx,yh-Gap*(j-0.5)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','right','Interpreter','latex');
           end;
        end;
    case 'SW';
        for j=1:nT;
            if j<=4;
                plot([xl xl+0.2*rx],[yl yl]+Gap*(nT-j+1)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',LinClr{j});
                text(xl+0.25*rx,yl+Gap*(nT-j+1)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','left','Interpreter','latex');
           elseif j<=8;
               plot([xl xl+0.2*rx],[yl yl]+Gap*(nT-j+1)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',ones(3,1)*0.5,'color',LinClr{j});
               text(xl+0.25*rx,yl+Gap*(nT-j+1)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','left','Interpreter','latex');
            end;
        end;
    case 'SE';
        for j=1:nT;
           if j<=4;
               plot([xh-0.2*rx xh],[yl yl]+Gap*(nT-j+1)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',LinClr{j});
               text(xh-0.25*rx,yl+Gap*(nT-j+1)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','right','Interpreter','latex');
           elseif j<=8;
               plot([xh-0.2*rx xh],[yl yl]+Gap*(nT-j+1)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',ones(3,1)*0.5,'color',LinClr{j});
               text(xh-0.25*rx,yl+Gap*(nT-j+1)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','right','Interpreter','latex');
           end;
        end;
    case 'MW';
        ym=(yl+yh)/2+Gap*ry*nT/2;
        for j=1:nT;
            if j<=4;
                plot([xl xl+0.2*rx],ym+[0 0]-Gap*(j-0.5)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',LinClr{j});
                text(xl+0.25*rx,ym-Gap*(j-0.5)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','left','Interpreter','latex');
            elseif j<=8;
                plot([xl xl+0.2*rx],ym+[0 0]-Gap*(j-0.5)*ry,'color','k','linestyle',LinStl{j},'linewidth',LinWdt,'marker','o','color',ones(3,1)*0.5,'color',LinClr{j});
                text(xl+0.25*rx,ym-Gap*(j-0.5)*ry,Txt{j},'FontSize',FntSiz,'FontName','Garamond','HorizontalAlignment','left','Interpreter','latex');
            end;                
        end;
end;

return;