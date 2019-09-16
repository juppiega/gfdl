function plot_power_psi()

load streamfun_months.mat

 plot_scatter(P, psimax_e, 'P [W/m^2]', 'M^E [Tg/s]',false);
L = 2*pi*6371*L/360;
 xtick = plot_scatter(1./L, psimax_e, '1/L [1/km]', 'M^E [Tg/s]',false);

plot_scatter(P, -psimax_w, 'P [W/m^2]', '-M^W [Tg/s]', true);
plot_scatter(1./L, -psimax_w, '1/L [1/km]', '-M^W [Tg/s]', true, xtick);

end

function xtick = plot_scatter(x,y,xname,yname,monsoon_ignore, xtick)

mo = 1:12;

if nargin > 4 && monsoon_ignore
    mo = mo(4:end);
    x_mon = x(1:3);
    y_mon = y(1:3);
    x = x(4:end);
    y = y(4:end);
end
figure;
scatter(x,y,nan)
fs = 15;
set(gca, 'fontsize', fs)
xlabel(xname, 'fontsize', fs)
ylabel(yname, 'fontsize', fs)

c = cellstr(num2str(mo'));
text(x, y, c,'fontsize',13)
if monsoon_ignore
    c = cellstr(num2str((1:3)'));
    text(x_mon, y_mon, c,'fontsize',13,'color','r')
end
xlim([min(x)*0.9, max(x)*1.1])
if all(y < 0)
    ylim([min(y)*1.1, max(y)*0.7])
else
    ylim([min(y)*0.9, max(y)*1.1])
end

[c,p] = corr(x,y);
if c > 0
    dim = [.15 .1 .1 .8];
else
    dim = [.15 .1 .1 .25];
end

if monsoon_ignore
    str = ['w/o winter monsoon:',char(10),'$\mathrm{r = ',num2str(c,'%4.2f'),'}$ ',char(10),'$\mathrm{p = ',num2str(p,'%4.2E'),'}$'];
else
    str = ['$\mathrm{r = ',num2str(c,'%4.2f'),'}$ ',char(10),'$\mathrm{p = ',num2str(p,'%4.2E'),'}$'];
end
annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',fs,'interpreter','latex');
if strcmpi(xname, '1/L [1/km]')
    xlim([5.3,8.5]*1E-5);
    if nargin > 5 
        set(gca,'xtick',xtick)
    end
end
xtick = get(gca,'xtick');
annotation('textbox',[.85,.83,.1,.1],'String','b)','FitBoxToText','on','fontsize',fs,'edgecolor','none')

end