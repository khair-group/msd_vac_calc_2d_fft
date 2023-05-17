
clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
xlabel('$t$','FontSize',30,'Interpreter','latex');
ylabel('$\left<\mathbf{r}^2\right>$','FontSize',30,'Interpreter','latex');


axes1.XScale='log';
axes1.YScale='log';

xlimval=5e2;


xlim([2e-1 xlimval]);
ylim([1e-4 1e2]);


dt=0.1;
L=10;


msd_imp=importdata('msd_tot_from_fft.dat');
pl_3=plot(msd_imp(:,1),msd_imp(:,2),'o');
pl_3.LineWidth=3;
pl_3.Color=[0. 0.7 0.6];
pl_3.DisplayName='numerical data';


D_r=3e-1;
kappa=D_r;
bet=5e-1;
sum_exp=kappa+bet;
v0=0.1;
del=0.05;

mean_vel=v0;
var_vel_mag=del*del/3;

a=2.*mean_vel*mean_vel/(kappa*kappa);
c=2.*var_vel_mag./(sum_exp*sum_exp);
f=@(x) a*(kappa*x-1+exp(-kappa*x))+c*(sum_exp*x-1+exp(-sum_exp*x));

fun1=fplot(f,[0. xlimval],'LineWidth',2,'LineStyle','-','Color',[255./255, 142./255, 28./255]);
fun1.DisplayName='analytical result';


%%%%%%%%%% drawing slope line %%%%%%%%%%%%%%%%%%%

C0=1e-2;
expon=2;
f=@(x) C0*(x^(expon));
fun=fplot(f,[0.1 xlimval],'LineWidth',2,'LineStyle','-');
fun.HandleVisibility='off';
fun.Color='k';


sl1_lc='k';
sl1_words=text(0.25,0.0013,'slope 2','Interpreter','latex');
sl1_words.FontSize=22;
sl1_words.Rotation=37;
sl1_words.Color='k';


C0=9e-2;
expon=1;
f=@(x) C0*(x^(expon));
fun=fplot(f,[0.2 xlimval],'LineWidth',2,'LineStyle','-');
fun.HandleVisibility='off';
fun.Color='k';


sl2_lc='k';
sl2_words=text(101,20,'slope 1','Interpreter','latex');
sl2_words.FontSize=22;
sl2_words.Rotation=19;
sl2_words.Color='k';

dim = [0.48 0.33 0.3 0.3];


leg=legend;
leg.Location='southeast';
leg.Interpreter='latex';
leg.FontSize=22;
leg.NumColumns=1;
