function Rb_transitions

[out,hF_S,hF_P1,hF_P3]=Rb_zeeman;


Es12=out.Es12;
B=out.B;

%% RF

hF1=figure(2003);
set(gcf,'color','w');
clf
axes
hold on

strs={'$|1,1\rangle\rightarrow |1,0\rangle$',...
    '$|1,0\rangle \rightarrow |1,-1\rangle$',...
    '$|2,2\rangle \rightarrow |2,1\rangle$',...
    '$|2,1\rangle \rightarrow |2,0\rangle$',...
    '$|2,0\rangle \rightarrow |2,-1\rangle$',...
    };

for kk=1:2
    plot(B,Es12(kk+1,:)-Es12(kk,:),'-','linewidth',2);
end

plot(B,Es12(end,:)-Es12(end-1,:),'-.','linewidth',2);
plot(B,Es12(end-1,:)-Es12(end-2,:),'-.','linewidth',2);
plot(B,Es12(end-2,:)-Es12(end-3,:),'-.','linewidth',2);

    
xlim([19 20]);

hF1.Position(3:4)=[600 300];

set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');

legend(strs,'interpreter','latex','fontsize',8,'location','northwest');

text(0.98,.02,'$^{87}\mathrm{Rb}~5\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','bottom','fontsize',18,'horizontalalignment','right');


%% uWave
hF2=figure(2004);

ind=find(B>30,1);

myfit=fittype('A*x+b','independent','x','coefficients',{'A','b'});

opt=fitoptions(myfit);
opt.StartPoint=[1.4 6840];
fout=fit(B(1:ind)',(Es12(end,1:ind)-Es12(1,1:ind))',myfit,opt);
str=['$' num2str(fout.A,'%.2f') '~\mathrm{MHz/G} + ' num2str(fout.b,'%.2f') '~\mathrm{MHz}$'];


set(gcf,'color','w');
clf

strs={'$|2,2\rangle\rightarrow|1,1\rangle$',str};



hold on
clear ps

co=get(gca,'colororder');

p1=plot(0:30,feval(fout,0:30),'-','linewidth',6,'color',[co(1,:) .4]);

% |2,2>-->|1,1>
p2=plot(B,Es12(end,:)-Es12(1,:),'k-','linewidth',2);
xlim([0 30]);

hF2.Position(3:4)=[600 300];

set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');

legend([p2 p1],strs,'interpreter','latex','fontsize',10,'location','northwest');
text(0.98,.02,'$^{87}\mathrm{Rb}~5\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','bottom','fontsize',18,'horizontalalignment','right');

%%
doSave=1;
if doSave
    fprintf('saving figures ...');
    
    % Save to png
    print(hF1,'Rb/Rb_RF.png','-dpng','-r400'); 
    print(hF2,'Rb/Rb_uWave.png','-dpng','-r400');  
    print(hF_S,'Rb/Rb_S12.png','-dpng','-r400'); 
    print(hF_P1,'Rb/Rb_P12.png','-dpng','-r400'); 
    print(hF_P3,'Rb/Rb_P32.png','-dpng','-r400'); 
    disp('done');
end
end
