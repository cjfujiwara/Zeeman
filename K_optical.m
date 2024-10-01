function K_optical(out)

Blim=[120 130];
% Blim=[190 195];

B = out.B;
e_P32_A = out.Ep32(1,:);
e_S12_A = out.Es12(1,:);
dE_AA = (e_P32_A-e_S12_A)-(e_P32_A(1)-e_S12_A(1));

e_P32_L = out.Ep32(12,:);
e_S12_J = out.Es12(10,:);
dE_LJ = (e_P32_L-e_S12_J)-(e_P32_L(1)-e_S12_J(1));

hF=figure(1100);
hF.Color='w';
clf
co=get(gca,'colororder');
subplot(121);
set(gca,'PositionConstraint','outerposition','box','on','linewidth',1,'fontsize',12)
plot(B,dE_AA,'linewidth',2,'color',co(1,:)); hold on;
xlabel('field (G)');
ylabel('$\Delta f$ from B=0 (MHz)','interpreter','latex')
xlim(Blim);
legend({'$a\rightarrow a`$'},'interpreter','latex')
 
subplot(122);
set(gca,'PositionConstraint','outerposition','box','on','linewidth',1,'fontsize',12)
plot(B,dE_LJ,'linewidth',2,'color',co(2,:)); 
xlabel('field (G)');
ylabel('$\Delta f$ from B=0 (MHz)','interpreter','latex')
xlim(Blim);
legend({'$j\rightarrow l`$'},'interpreter','latex')



% plot(out.B,)
% 

% %% Optical High Field
% % Plot the high field imaging transitions
% Brange=[195 205];
% 
% hF4=figure(1006);
% hF4.Position(1:2)=[650 150];
% hF4.Position(3:4)=[450 600];
% set(gcf,'color','w');
% clf
% 
% % Plot each transitions for the |9,-9>-->|11,-11> and |9,-7>-->|11,-9>
% subplot(211);
% hold on
% clear ps
% co=get(gca,'colororder');
% dE1=Ep32(1,:)-Es12(1,:)-(Ep32(1,1)-Es12(1,1));
% dE2=Ep32(2,:)-Es12(2,:)-(Ep32(2,1)-Es12(2,1));
% 
% dE1=dE1';
% dE2=dE2';
% 
% ind1=find(B>Brange(1),1);
% ind2=find(B>Brange(2),1);
% 
% myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
% opt1=fitoptions(myfit);
% opt2=fitoptions(myfit);
% 
% 
% opt1.StartPoint=[-1.4 dE1(1)];
% opt2.StartPoint=[-1.4 dE2(1)];
% 
% fout1=fit(B(ind1:ind2)',dE1(ind1:ind2),myfit,opt1);
% fout2=fit(B(ind1:ind2)',dE2(ind1:ind2),myfit,opt2);
% 
% strFit1=['$' num2str(fout1.A,'%.3f') '~\mathrm{MHz/G}' ...
%     '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout1.b,'%.3f') '~\mathrm{MHz}$'];
% strFit2=['$' num2str(fout2.A,'%.3f') '~\mathrm{MHz/G}' ...
%     '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout2.b,'%.3f') '~\mathrm{MHz}$'];
% 
% disp(strFit1);
% disp(strFit2);
% 
% 
% 
% 
% 
% p2=plot(B,dE1,'-','linewidth',2);
% p3=plot(B,dE2,'-','linewidth',2);
% strs={['$|9/2,-9/2\rangle\rightarrow|11,-11/2\rangle,~(' num2str(fout1.A) '~\mathrm{MHz/G})$'],...
%     ['$|9/2,-7/2\rangle\rightarrow|11,-9/2\rangle,~(' num2str(fout2.A) '~\mathrm{MHz/G})$']};
% xlim([max([0 Brange(1)]) Brange(2)]);
% set(gca,'fontsize',12,'fontname','times','xgrid','on',...
%     'box','on','ygrid','on');
% xlabel('field (Gauss)');
% ylabel('energy relative to zero field (MHz)');
% legend([p2 p3],strs,'interpreter','latex','fontsize',8,'location','best');
% text(0.02,.02,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}\rightarrow^{40}\mathrm{K}~4\mathrm{P}_{3/2}$','interpreter','latex','units','normalized',...
%     'verticalalignment','bottom','fontsize',14,'horizontalalignment','left');
% % xlim([1);
% 
% 
% % Fitting
% ind1=find(B>Brange(1),1);
% ind2=find(B>Brange(2),1);
% 
% myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
% % opt.StartPoint=[.2 Es12(end,ind1)-Es12(1,ind1)];
% % keyboard
% fout=fit(B(ind1:ind2)',dE1(ind1:ind2)-dE2(ind1:ind2),myfit,opt);
% strFit=['$' num2str(fout.A,'%.3f') '~\mathrm{MHz/G}' ...
%     '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];
% 
% 
% % Plot the relative energy transitions for the |9,-9>-->|11,-11> and |9,-7>-->|11,-9>
% subplot(212);
% plot(B,dE1-dE2,'k-','linewidth',2);
% set(gca,'fontsize',12,'fontname','times','xgrid','on',...
%     'box','on','ygrid','on');
% 
% hold on
% 
% plot(B(ind1:ind2),feval(fout,B(ind1:ind2)),'r--');
% 
% xlim(Brange);
% xlabel('field (Gauss)');
% ylabel('relative energy (MHz)');
% 
% text(0.02,0.02,strFit,'units','normalized','interpreter','latex','verticalalignment','bottom');
end

