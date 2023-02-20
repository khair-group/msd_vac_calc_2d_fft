%%% 18-FEB-2023: See


clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
xlabel('$\tau$','FontSize',30,'Interpreter','latex');
ylabel('MSD','FontSize',30,'Interpreter','latex');
% title('Poisson, PBC, no coupling, steric repulsion','Interpreter','latex','FontSize',18);


axes1.XScale='log';
axes1.YScale='log';

xlimval=5e2;

xlim([8e-1 xlimval]);
% ylim([1e-5 1e2]);


dt=1.;
L=10;

msd_imp=importdata('msd_tot_from_fft.dat');


pl_1=plot(msd_imp(:,1),msd_imp(:,2),'-');
pl_1.LineWidth=3;
pl_1.Color='r';
pl_1.DisplayName='$N=100$';

vseries=importdata('vel_tseries_mu_0.10_del_0.05_beta_0.500_alph_0.300.dat');
[vel_mag] = ret_vel_mag(vseries(:,2),vseries(:,3));
alph=0.3;
bet=0.5;
sum_exp=alph+bet;
mean_vel=mean(vel_mag);
a=2.*mean_vel*mean_vel/(alph*alph);
c=2.*var(vel_mag)./(sum_exp*sum_exp);
f=@(x) a*(alph*x-1+exp(-alph*x))+c*(sum_exp*x-1+exp(-sum_exp*x));

fun1=fplot(f,[0. xlimval],'LineWidth',2,'LineStyle','-.','Color','k');
fun1.DisplayName='single particle result';

dim = [0.48 0.33 0.3 0.3];


str = {'$\alpha=0.3$, $\beta=0.5$', '$\mu=0.1$, $\delta=0.05$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'Interpreter','latex','FontSize',24,'EdgeColor','k',...
    'Color','k');

str = {'Phantom disks'};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'Interpreter','latex','FontSize',24,'EdgeColor','k',...
    'Color','k');

leg=legend;
leg.Location='southeast';
leg.Interpreter='latex';
leg.FontSize=22;
leg.NumColumns=1;