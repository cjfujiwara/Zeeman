function K_transitions

[out,hF_S,hF_P1,hF_P3]=K_zeeman;

hF_S.Position(1:2)=[10 50];
hF_P1.Position(1:2)=[10 300];
hF_P3.Position(1:2)=[10 600];

Es12=out.Es12;
Ep32=out.Ep32;

B=out.B;
%% RF Transitions
hF1=figure(1003);
hF1.Position(1:2)=[650 50];
hF1.Position(3:4)=[600 300];
set(gcf,'color','w');
clf
axes
hold on
Brange=[19.7 20.3];
Brange=[9.7 10.3];

ind1=find(B>Brange(1),1);
ind2=find(B>Brange(2),1);
myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});

opt=fitoptions(myfit);
opt.StartPoint=[.2 Es12(2,ind1)-Es12(1,ind1)];
fout=fit(B(ind1:ind2)',(Es12(2,ind1:ind2)-Es12(1,ind1:ind2))',myfit,opt);
strFit=['$' num2str(1E3*fout.A,'%.2f') '~\mathrm{kHz/G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];

strs={'$-9/2\rightarrow-7/2$',
    '$-7/2\rightarrow-5/2$',
    '$-5/2\rightarrow-3/2$',
    '$-3/2\rightarrow-1/2$',
    '$-1/2\rightarrow+1/2$',
    '$+1/2\rightarrow+3/2$',
    '$+3/2\rightarrow+5/2$',
    '$+5/2\rightarrow+7/2$',
    '$+7/2\rightarrow+9/2$',
    };
co=get(gca,'colororder');

% p1=plot(Brange(1):Brange(2),feval(fout,Brange(1):Brange(2)),'-','linewidth',6,'color',[co(2,:) .4]);
hold on

inds=[1:9];
for kk=1:length(inds)
    plot(B,Es12(kk+1,:)-Es12(kk,:),'-','linewidth',2,'color',co(mod(kk-1,7)+1,:));
end
xlim(Brange);

set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
text(0.98,.02,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','bottom','fontsize',18,'horizontalalignment','right');
hold on


% legend([{strFit} strs(inds)],'interpreter','latex','fontsize',8,'location','northwest');
legend(strs(inds),'interpreter','latex','fontsize',8,'location','northwest');

%% RF High Field Transitions
hF1b=figure(1020);
hF1b.Position(1:2)=[650 50];
hF1b.Position(3:4)=[1000 300];
set(gcf,'color','w');
clf

Brange=[195 204];

ind1=find(B>Brange(1),1);
ind2=find(B>Brange(2),1);



strs={'$-9/2\rightarrow-7/2$',
    '$-7/2\rightarrow-5/2$',
    '$-5/2\rightarrow-3/2$',
    '$-3/2\rightarrow-1/2$',
    '$-1/2\rightarrow+1/2$',
    '$+1/2\rightarrow+3/2$',
    '$+3/2\rightarrow+5/2$',
    '$+5/2\rightarrow+7/2$',
    '$+7/2\rightarrow+9/2$',
    };
co=get(gca,'colororder');

% 

% -9-->-7
subplot(131);
inds=[1];

myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
opt=fitoptions(myfit);
opt.StartPoint=[.2 Es12(inds(1)+1,ind1)-Es12(inds(1),ind1)];
fout=fit(B(ind1:ind2)',(Es12(inds(1)+1,ind1:ind2)-Es12(inds(1),ind1:ind2))',myfit,opt);
strFit=['$' num2str(1E3*fout.A,'%.2f') '~\mathrm{kHz}/\mathrm{G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];


hold on
for kk=1:length(inds)
    plot(B,Es12(inds(kk)+1,:)-Es12(inds(kk),:),'-','linewidth',2,'color',co(mod(inds(kk)-1,7)+1,:));
    hold on
end
plot(Brange(1):Brange(2),feval(fout,Brange(1):Brange(2)),'--','linewidth',2,'color',[0 0 0 1]);

xlim(Brange);
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
% legend({strs{inds},strFit},'interpreter','latex','fontsize',8,'location','southeast');
legend(strs(inds),'interpreter','latex','fontsize',8,'location','northwest');

text(.05,.01,strFit,'units','normalized','fontsize',8,'interpreter','latex',...
    'verticalalignment','bottom')


% -9-->-7
subplot(132);
inds=[2];

myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
opt=fitoptions(myfit);
opt.StartPoint=[.2 Es12(inds(1)+1,ind1)-Es12(inds(1),ind1)];
fout=fit(B(ind1:ind2)',(Es12(inds(1)+1,ind1:ind2)-Es12(inds(1),ind1:ind2))',myfit,opt);
strFit=['$' num2str(1E3*fout.A,'%.2f') '~\mathrm{kHz}/\mathrm{G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];

for kk=1:length(inds)
    plot(B,Es12(inds(kk)+1,:)-Es12(inds(kk),:),'-','linewidth',2,'color',co(mod(inds(kk)-1,7)+1,:));
    hold on
end
plot(Brange(1):Brange(2),feval(fout,Brange(1):Brange(2)),'--','linewidth',2,'color',[0 0 0 1]);

xlim(Brange);
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
legend(strs(inds),'interpreter','latex','fontsize',7,'location','northwest');
text(.05,.01,strFit,'units','normalized','fontsize',8,'interpreter','latex',...
    'verticalalignment','bottom')

% -7-->-5
subplot(133);
inds=[3];
myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
opt=fitoptions(myfit);
opt.StartPoint=[.2 Es12(inds(1)+1,ind1)-Es12(inds(1),ind1)];
fout=fit(B(ind1:ind2)',(Es12(inds(1)+1,ind1:ind2)-Es12(inds(1),ind1:ind2))',myfit,opt);
strFit=['$' num2str(1E3*fout.A,'%.2f') '~\mathrm{kHz}/\mathrm{G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];
for kk=1:length(inds)
    plot(B,Es12(inds(kk)+1,:)-Es12(inds(kk),:),'-','linewidth',2,'color',co(mod(inds(kk)-1,7)+1,:));
    hold on
end
plot(Brange(1):Brange(2),feval(fout,Brange(1):Brange(2)),'--','linewidth',2,'color',[0 0 0 1]);
xlim(Brange);
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
legend(strs(inds),'interpreter','latex','fontsize',8,'location','northwest');
text(.05,.01,strFit,'units','normalized','fontsize',8,'interpreter','latex',...
    'verticalalignment','bottom')


% legend([{strFit} strs(inds)],'interpreter','latex','fontsize',8,'location','northwest');

%% RF High Field Transitions
hF1b=figure(1120);
hF1b.Position(1:2)=[650 50];
hF1b.Position(3:4)=[500 500];
set(gcf,'color','w');
clf

Brange=124+[-.5 .5];

co=get(gca,'colororder');

subplot(121);

% -9-->-7
inds=[1 2 3 4 5 6 7 8 9 17 16 15 14 13 12 11];
come=jet(length(inds));

hold on
legStr={};
for kk=1:length(inds)
    
    if inds(kk)>9
        linespec='--';
    else
        linespec='-';
    end
    plot(B,Es12(inds(kk)+1,:)-Es12(inds(kk),:),'linewidth',2,'color',come(kk,:),...
        'linestyle',linespec);
    hold on
    legStr{kk}=['$' char(96+inds(kk)) '\longleftrightarrow ' char(96+inds(kk)+1) '$'];
end

xlim(Brange);
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
% legend({strs{inds},strFit},'interpreter','latex','fontsize',8,'location','southeast');
legend(legStr,'interpreter','latex','fontsize',8,'location','northwest');


subplot(122);
inds=[1 18;%-9/2-->-7/2
    2 18;%-7/2-->-7/2
    2 17;
    3 17;
    ];%-7/2-->-5/2
% come=jet(size(inds,1));
hold on
come=get(gca,'colororder');
legStr={};
for kk=1:length(inds)
    plot(B,Es12(inds(kk,2),:)-Es12(inds(kk,1),:),'linewidth',2,'color',co(kk,:),...
        'linestyle','-');
    hold on
    legStr{kk}=['$' char(96+inds(kk,1)) '\longleftrightarrow ' char(96+inds(kk,2)) '$'];
end

xlim(Brange);
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');
% legend({strs{inds},strFit},'interpreter','latex','fontsize',8,'location','southeast');
legend(legStr,'interpreter','latex','fontsize',8,'location','best');
%% uWave Low Field

Brange=[-1 30];

ind1=find(B>Brange(1),1);
ind2=find(B>Brange(2),1);

myfit=fittype('A*x+b','independent','x','coefficients',{'A','b'});

opt=fitoptions(myfit);
opt.StartPoint=[2.51 1285+2.51*min(Brange)];
fout=fit(B(ind1:ind2)',(Es12(end,ind1:ind2)-Es12(1,ind1:ind2))',myfit,opt);
str=['$' num2str(fout.A,'%.2f') '~\mathrm{MHz/G} + ' num2str(fout.b,'%.1f') '~\mathrm{MHz}$'];

hF2=figure(1004);
hF2.Position(1:2)=[600 100];
hF2.Position(3:4)=[600 300];

set(gcf,'color','w');
clf

strs1={'$|9/2,-9/2\rangle\rightarrow|7/2,-7/2\rangle$',str};



hold on
clear ps

co=get(gca,'colororder');

p1=plot(Brange(1):Brange(2),feval(fout,Brange(1):Brange(2)),'-','linewidth',6,'color',[co(1,:) .4]);

% |9,-9>-->|7,-7>
p2=plot(B,Es12(end,:)-Es12(1,:),'k-','linewidth',2);
xlim([max([0 Brange(1)]) Brange(2)]);


set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');

legend([p2 p1],strs1,'interpreter','latex','fontsize',10,'location','northwest');


text(0.98,.02,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','bottom','fontsize',18,'horizontalalignment','right');
%% uWave Plane Selection

Brange=[120 135];

ind1=find(B>Brange(1),1);
ind2=find(B>Brange(2),1);

% myfit=fittype('A*x+b','independent','x','coefficients',{'A','b'});
% 
% opt=fitoptions(myfit);
% opt.StartPoint=[2.51 1285+2.51*min(Brange)];
% fout=fit(B(ind1:ind2)',(Es12(end,ind1:ind2)-Es12(1,ind1:ind2))',myfit,opt);
% str=['$' num2str(fout.A,'%.2f') '~\mathrm{MHz/G} + ' num2str(fout.b,'%.1f') '~\mathrm{MHz}$'];


myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
opt=fitoptions(myfit);
opt.StartPoint=[.2 Es12(end,ind1)-Es12(1,ind1)];
fout=fit(B(ind1:ind2)',(Es12(end,ind1:ind2)-Es12(1,ind1:ind2))',myfit,opt);
strFit=['$' num2str(fout.A,'%.3f') '~\mathrm{MHz/G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];




hF3=figure(1005);
set(gcf,'color','w');
clf

strs1={'$|9/2,-9/2\rangle\rightarrow|7/2,-7/2\rangle$',strFit};



hold on
clear ps

co=get(gca,'colororder');

p1=plot(Brange(1):Brange(2),feval(fout,Brange(1):Brange(2)),'-','linewidth',6,'color',[co(1,:) .4]);

% |9,-9>-->|7,-7>
p2=plot(B,Es12(end,:)-Es12(1,:),'k-','linewidth',2);
xlim([max([0 Brange(1)]) Brange(2)]);

hF3.Position(3:4)=[600 300];

set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');

legend([p2 p1],strs1,'interpreter','latex','fontsize',10,'location','northwest');


text(0.98,.02,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','bottom','fontsize',18,'horizontalalignment','right');

%% uWave Plane Selection

Brange=[120 135];

ind1=find(B>Brange(1),1);
ind2=find(B>Brange(2),1);

% myfit=fittype('A*x+b','independent','x','coefficients',{'A','b'});
% 
% opt=fitoptions(myfit);
% opt.StartPoint=[2.51 1285+2.51*min(Brange)];
% fout=fit(B(ind1:ind2)',(Es12(end,ind1:ind2)-Es12(1,ind1:ind2))',myfit,opt);
% str=['$' num2str(fout.A,'%.2f') '~\mathrm{MHz/G} + ' num2str(fout.b,'%.1f') '~\mathrm{MHz}$'];


myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
opt=fitoptions(myfit);
opt.StartPoint=[.2 Es12(end-1,ind1)-Es12(2,ind1)];
fout=fit(B(ind1:ind2)',(Es12(end-1,ind1:ind2)-Es12(2,ind1:ind2))',myfit,opt);
strFit=['$' num2str(fout.A,'%.3f') '~\mathrm{MHz/G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];




hF5=figure(1007);
set(gcf,'color','w');
clf

strs1={'$|9/2,-7/2\rangle\rightarrow|7/2,-5/2\rangle$',strFit};



hold on
clear ps

co=get(gca,'colororder');

p1=plot(Brange(1):Brange(2),feval(fout,Brange(1):Brange(2)),'-','linewidth',6,'color',[co(1,:) .4]);

% |9,-9>-->|7,-7>
p2=plot(B,Es12(end-1,:)-Es12(2,:),'k-','linewidth',2);
xlim([max([0 Brange(1)]) Brange(2)]);

hF5.Position(3:4)=[600 300];

set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy (MHz)');

legend([p2 p1],strs1,'interpreter','latex','fontsize',10,'location','northwest');


text(0.98,.02,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','bottom','fontsize',18,'horizontalalignment','right');
%% Optical High Field
% Plot the high field imaging transitions
Brange=[195 205];

hF4=figure(1006);
hF4.Position(1:2)=[650 150];
hF4.Position(3:4)=[450 600];
set(gcf,'color','w');
clf

% Plot each transitions for the |9,-9>-->|11,-11> and |9,-7>-->|11,-9>
subplot(211);
hold on
clear ps
co=get(gca,'colororder');
dE1=Ep32(1,:)-Es12(1,:)-(Ep32(1,1)-Es12(1,1));
dE2=Ep32(2,:)-Es12(2,:)-(Ep32(2,1)-Es12(2,1));

dE1=dE1';
dE2=dE2';

ind1=find(B>Brange(1),1);
ind2=find(B>Brange(2),1);

myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
opt1=fitoptions(myfit);
opt2=fitoptions(myfit);


opt1.StartPoint=[-1.4 dE1(1)];
opt2.StartPoint=[-1.4 dE2(1)];

fout1=fit(B(ind1:ind2)',dE1(ind1:ind2),myfit,opt1);
fout2=fit(B(ind1:ind2)',dE2(ind1:ind2),myfit,opt2);

strFit1=['$' num2str(fout1.A,'%.3f') '~\mathrm{MHz/G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout1.b,'%.3f') '~\mathrm{MHz}$'];
strFit2=['$' num2str(fout2.A,'%.3f') '~\mathrm{MHz/G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout2.b,'%.3f') '~\mathrm{MHz}$'];

disp(strFit1);
disp(strFit2);





p2=plot(B,dE1,'-','linewidth',2);
p3=plot(B,dE2,'-','linewidth',2);
strs={['$|9/2,-9/2\rangle\rightarrow|11,-11/2\rangle,~(' num2str(fout1.A) '~\mathrm{MHz/G})$'],...
    ['$|9/2,-7/2\rangle\rightarrow|11,-9/2\rangle,~(' num2str(fout2.A) '~\mathrm{MHz/G})$']};
xlim([max([0 Brange(1)]) Brange(2)]);
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');
xlabel('field (Gauss)');
ylabel('energy relative to zero field (MHz)');
legend([p2 p3],strs,'interpreter','latex','fontsize',8,'location','best');
text(0.02,.02,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}\rightarrow^{40}\mathrm{K}~4\mathrm{P}_{3/2}$','interpreter','latex','units','normalized',...
    'verticalalignment','bottom','fontsize',14,'horizontalalignment','left');
% xlim([1);


% Fitting
ind1=find(B>Brange(1),1);
ind2=find(B>Brange(2),1);

myfit=fittype(@(A,b,x) A*(x-Brange(1))+b,'independent','x','coefficients',{'A','b'});
% opt.StartPoint=[.2 Es12(end,ind1)-Es12(1,ind1)];
% keyboard
fout=fit(B(ind1:ind2)',dE1(ind1:ind2)-dE2(ind1:ind2),myfit,opt);
strFit=['$' num2str(fout.A,'%.3f') '~\mathrm{MHz/G}' ...
    '(B-' num2str(Brange(1)) '~\mathrm{G}) + ' num2str(fout.b,'%.3f') '~\mathrm{MHz}$'];


% Plot the relative energy transitions for the |9,-9>-->|11,-11> and |9,-7>-->|11,-9>
subplot(212);
plot(B,dE1-dE2,'k-','linewidth',2);
set(gca,'fontsize',12,'fontname','times','xgrid','on',...
    'box','on','ygrid','on');

hold on

plot(B(ind1:ind2),feval(fout,B(ind1:ind2)),'r--');

xlim(Brange);
xlabel('field (Gauss)');
ylabel('relative energy (MHz)');

text(0.02,0.02,strFit,'units','normalized','interpreter','latex','verticalalignment','bottom');
%%
doSave=0;
if doSave
    fprintf('saving figures ...');
    
    % Save to png
    print(hF1,'K/K_RF.png','-dpng','-r400'); 
    print(hF1b,'K/K_RF_high_field.png','-dpng','-r800'); 

    print(hF2,'K/K_uWave.png','-dpng','-r400');  
    print(hF3,'K/K_uWave_high.png','-dpng','-r400');  
    print(hF4,'K/K_optical_high.png','-dpng','-r400');  

    print(hF_S,'K/K_S12.png','-dpng','-r400'); 
    print(hF_P1,'K/K_P12.png','-dpng','-r400'); 
    print(hF_P3,'K/K_P32.png','-dpng','-r400'); 
    disp('done');
end
end

