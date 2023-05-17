clf;
axes1 = axes;

hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',2,'TickLength',[0.015 0.025]);
xlabel('$t$','FontSize',30,'Interpreter','latex');
ylabel('$\left<\mathbf{v}(t)\cdot\mathbf{v}(0)\right>/\left<v\right>^2$','FontSize',30,'Interpreter','latex');


% axes1.XScale='log';
% axes1.YScale='log';

xlimval=8;

xlim([0. xlimval]);



dt=0.1;
L=10;


v0=0.1;
del=0.05;
mean_vel=v0;
var_vel_mag=del*del/3;
D_r=5e-1;
kappa=D_r;
bet=2e-1;

vac_imp=importdata('vac_from_fft.dat');
pl_3=plot(vac_imp(:,1),vac_imp(:,2)./(mean_vel*mean_vel),'o');
pl_3.LineWidth=3;
pl_3.Color=[0. 0.7 0.6];
pl_3.DisplayName='numerical data';



a=mean_vel*mean_vel;
b=-kappa;
c=var_vel_mag;
d=-(bet+kappa);

nmlz=mean_vel*mean_vel;
a=a./nmlz;
c=c./nmlz;

f=@(x) a*exp(b*x)+c*exp(d*x);
fun1=fplot(f,[0. xlimval],'LineWidth',2,'LineStyle','-','Color',[255./255, 142./255, 28./255]);
fun1.DisplayName='analytical result';


 
% 
leg=legend;
leg.Location='northeast';
leg.Interpreter='latex';
leg.FontSize=22;
leg.NumColumns=1;