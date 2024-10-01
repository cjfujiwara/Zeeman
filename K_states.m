function K_states(out)


% Get the magnetic field
B=out.B;

%% 4S1/2 State Tracking

    function [amplitudes,basis_mImJ,str]=getState_S12(eigenstate_number)
        v = out.vs12(:,eigenstate_number,:);    % The Vector
        v_last = out.vs12(:,eigenstate_number,end); % Last Vector (should be paschen-back)
        
        % Make string |mI,mJ> for last magnetic field
        str = ['$a\rightarrow |' ...
            num2str(out.mImJ_S12(eigenstate_number,1)) ',' ...
            num2str(out.mImJ_S12(eigenstate_number,2)) '\rangle$'];        
        basis_inds = find(abs(v_last)~=0,2);
        
        amplitudes=zeros(numel(out.B),numel(basis_inds));
        basis_mImJ=zeros(numel(basis_inds),2);
        for kk=1:length(basis_inds)
            this_amplitude = out.vs12(basis_inds(kk),eigenstate_number,:);
            
            this_amplitude = reshape(this_amplitude,[numel(B) 1]);
            amplitudes(:,kk) = this_amplitude;
            basis_mImJ(kk,:)=out.basis_mImJ_S12(basis_inds(kk),:);
        end      
    end

    function [amplitudes,basis_mImJ,str]=getState_P12(eigenstate_number)
        v = out.vp12(:,eigenstate_number,:);    % The Vector
        v_last = out.vp12(:,eigenstate_number,end); % Last Vector (should be paschen-back)
        
        % Make string |mI,mJ> for last magnetic field
        str = ['$a\rightarrow |' ...
            num2str(out.mImJ_P12(eigenstate_number,1)) ',' ...
            num2str(out.mImJ_P12(eigenstate_number,2)) '\rangle$'];        
        basis_inds = find(abs(v_last)~=0,2);
        
        amplitudes=zeros(numel(out.B),numel(basis_inds));
        basis_mImJ=zeros(numel(basis_inds),2);
        for kk=1:length(basis_inds)
            this_amplitude = out.vp12(basis_inds(kk),eigenstate_number,:);
            
            this_amplitude = reshape(this_amplitude,[numel(B) 1]);
            amplitudes(:,kk) = this_amplitude;
            basis_mImJ(kk,:)=out.basis_mImJ_P12(basis_inds(kk),:);
        end
    end  

    function [amplitudes,basis_mImJ,str]=getState_P32(eigenstate_number,numstates)
        
        if nargin ==1
            numstates=2;
        end
        v = out.vp32(:,eigenstate_number,:);    % The Vector
        v_last = out.vp32(:,eigenstate_number,end); % Last Vector (should be paschen-back)
        
        % Make string |mI,mJ> for last magnetic field
        str = ['$a\rightarrow |' ...
            num2str(out.mImJ_P32(eigenstate_number,1)) ',' ...
            num2str(out.mImJ_P32(eigenstate_number,2)) '\rangle$'];        
        
        [~,bb] = sort(abs(v_last),'descend');
        bb=bb(1:numstates);
%         bb = sort(bb(1:numstates));
        basis_inds = bb;
        amplitudes=zeros(numel(out.B),numel(basis_inds));
        basis_mImJ=zeros(numel(basis_inds),2);
        
        for kk=1:length(basis_inds)
            this_amplitude = out.vp32(basis_inds(kk),eigenstate_number,:);
            
            this_amplitude = reshape(this_amplitude,[numel(B) 1]);
            amplitudes(:,kk) = this_amplitude;
            basis_mImJ(kk,:)=out.basis_mImJ_P32(basis_inds(kk),:);
        end      
    end
%%

f_s12=figure(10001);
f_s12.Color='w';
f_s12.Position=[50 50 1100 200];
clf
inds = [1 10 2 18];
legStr={};
co=get(gca,'colororder');
linespecs={'-','--','-.'};
for mm=1:length(inds)
    legStr={};
    subplot(1,length(inds),mm)
    [amplitudes,basis_mImJ,str]=getState_S12(inds(mm));
        set(gca,'PositionConstraint','outerposition');

    for rr=1:size(basis_mImJ)
        str=['$|\langle' num2str(basis_mImJ(rr,1), '%+i') ',' num2str(basis_mImJ(rr,2),'%+.1f') '|{\bf ' num2str(inds(mm),'%02i') '}\rangle|'];
        plot(B,abs(amplitudes(:,rr)),'color',co(mm,:),'linestyle',linespecs{rr},'linewidth',2);
        hold on;
        
        if inds(mm)<=26
            str = [str '~{\bf ' char(96+inds(mm)) '}$'];
        else
            str=[str '$'];
        end
        
        legStr{end+1} = str;



    end
    legend(legStr,'interpreter','latex','location','northeast','fontsize',8);
    set(gca,'box','on','linewidth',1,'fontsize',10,'fontname','times');
    xlabel('field (G)');
    ylabel('amplitude');
    xlim([0 230]);
    ylim([0 1]);
    
            if mm==1
            text(.02,.02,'$^{40}\mathrm{K}~4 S_{1/2}$','interpreter','latex',...
                'units','normalized','verticalalignment','bottom',...
                'horizontalalignment','left');
            end
end
%%

f_p12=figure(10002);
f_p12.Color='w';
f_p12.Position=[50 250 1100 200];
clf
inds = [1 10 2 18];
legStr={};
co=get(gca,'colororder');
linespecs={'-','--','-.'};
for mm=1:length(inds)
    legStr={};
    subplot(1,length(inds),mm)
    set(gca,'PositionConstraint','outerposition');

    [amplitudes,basis_mImJ,str]=getState_P12(inds(mm));
    
    for rr=1:size(basis_mImJ)
        str=['$|\langle' num2str(basis_mImJ(rr,1), '%+i') ',' num2str(basis_mImJ(rr,2),'%+.1f') '|{\bf ' num2str(inds(mm),'%02i') '}\rangle|'];
        plot(B,abs(amplitudes(:,rr)),'color',co(mm,:),'linestyle',linespecs{rr},'linewidth',2);
        hold on;        
        if inds(mm)<=26
            str = [str '~{\bf ' char(96+inds(mm)) '}$'];
        else
            str=[str '$'];
        end        
        legStr{end+1} = str;
    end
    legend(legStr,'interpreter','latex','location','northeast','fontsize',8);
    set(gca,'box','on','linewidth',1,'fontsize',10,'fontname','times');
    xlabel('field (G)');
    ylabel('amplitude');
    xlim([0 230]);
    ylim([0 1]);

  if mm==1
        text(.02,.02,'$^{40}\mathrm{K}~4 P_{1/2}$','interpreter','latex',...
            'units','normalized','verticalalignment','bottom',...
            'horizontalalignment','left');
    end
end

%%

f_p32=figure(10003);
f_p32.Color='w';
f_p32.Position=[50 350 1100 200];
clf
inds = [1 12 2 36]; % Track a,b,q,r states
numstates=[1 1 2 4];
legStr={};
co=get(gca,'colororder');
linespecs={'-','--','-.',':'};
for mm=1:length(inds)
    legStr={};
    subplot(1,length(inds),mm)
    set(gca,'PositionConstraint','outerposition');

    N=numstates(mm); 
    [amplitudes,basis_mImJ,str]=getState_P32(inds(mm),N);    
    for rr=1:size(basis_mImJ)
        str=['$|\langle' num2str(basis_mImJ(rr,1), '%+i') ',' num2str(basis_mImJ(rr,2),'%+.1f') '|{\bf ' num2str(inds(mm),'%02i') '}\rangle|'];
        plot(B,abs(amplitudes(:,rr)),'color',co(mm,:),'linestyle',linespecs{rr},'linewidth',2);
        hold on;        
        if inds(mm)<=26
            str = [str '~{\bf ' char(96+inds(mm)) '}$'];
        else
            str=[str '$'];
        end
        ylim([0 1]);
        legStr{end+1} = str;
    end
    legend(legStr,'interpreter','latex','location','northeast','fontsize',8);
    set(gca,'box','on','linewidth',1,'fontsize',10,'fontname','times');
    xlabel('field (G)');
    ylabel('amplitude');
    xlim([0 230]);
    ylim([0 1]);

          if mm==1
            text(.02,.02,'$^{40}\mathrm{K}~4 P_{3/2}$','interpreter','latex',...
                'units','normalized','verticalalignment','bottom',...
                'horizontalalignment','left');
        end
end
set(gca,'PositionConstraint','outerposition');


% %%
% % b = state #2 (9/2,-7/2)
% % r = state #18 (7/2,-7/2)
% % q = state #17 (7/5,-5/2)
% 
% % Get the |9,-9> state
% e_99=out.Es12(1,:);
% 
% % |mJ=-1/2>
% v_99_a=out.vs12(18,1,:);
% v_99_a=reshape(v_99_a,[length(B) 1]);
% 
% % Get the |9,-7> state
% e_97=out.Es12(2,:);
% 
% % |mJ=-1/2>
% v_97_a=out.vs12(17,2,:);    
% v_97_a=reshape(v_97_a,[length(B) 1]);
% 
% % |mJ=+1/2>
% v_97_b=out.vs12(9,2,:);                 
% v_97_b=reshape(v_97_b,[length(B) 1]);
% 
% % get the |r> state = |7,-7> state
% v_77_a = out.vs12(9,end,:);%v_77_a=v_77_a(:);
% v_77_b = out.vs12(17,end,:);%v_77_b=v_77_b(:);
% 
% 
% % keyboard
% %%
% % Get the |11,-11> state
% e_1111=out.Ep32(1,:);
% 
% 
% 
% % |11,-11> --> |mJ=-3/2>
% v_1111_a=out.vp32(36,1,:);
% v_1111_a=reshape(v_1111_a,[length(B) 1]);
% 
% % |11,-11> --> |mJ=-3/2>
% v_119_a=out.vp32(35,2,:);
% v_119_a=reshape(v_119_a,[length(B) 1]);
% 
% % |11,-11> --> |mJ=-1/2>
% v_119_b=out.vp32(27,2,:);
% v_119_b=reshape(v_119_b,[length(B) 1]);
% 
% hF=figure;
% hF.Color='w';
% hF.Position=[100 100 450 600];
% 
% str1='$|a\rangle = |9/2,-9/2\rangle\rightarrow |-4,-1/2\rangle$';
% str2='$|b\rangle =|9/2,-7/2\rangle\rightarrow {\bf a}|-3,-1/2\rangle+b|-4,+1/2\rangle$';
% str3='$|b\rangle= |9/2,-7/2\rangle\rightarrow a|-3,-1/2\rangle+{\bf b}|-4,+1/2\rangle$';
% 
% subplot(211);
% plot(B(2:end),v_99_a(2:end),'-','linewidth',1);
% hold on
% plot(B(2:end),abs(v_97_a(2:end)),'linewidth',2);
% plot(B(2:end),abs(v_97_b(2:end)),'linewidth',2);
% xlim([0 210]);
% ylim([0 1]);
% xlabel('magnetic field (gauss)');
% set(gca,'fontsize',12);
% 
% legend({str1,str2,str3},'interpreter','latex','fontsize',8);
% 
% str1='$|a`\rangle =|11/2,-11/2\rangle\rightarrow |-4,-3/2\rangle$';
% str2='$|b`\rangle=|11/2,-9/2\rangle\rightarrow {\bf a}|-3,-3/2\rangle+b|-4,-1/2\rangle$';
% str3='$|b`\rangle=|11/2,-9/2\rangle\rightarrow a|-3,-3/2\rangle+{\bf b}|-4,-1/2\rangle$';
% 
% text(0.02,.02,'$^{40}\mathrm{K}~4\mathrm{S}_{1/2}$','interpreter','latex','units','normalized',...
%     'verticalalignment','bottom','fontsize',18);
% 
% % str3='$\langle m_I=-3,m_J=-3/2|F=11/2,m_F-9/2\rangle$';
% 
% ylabel('amplitude');
% set(gca,'fontsize',12,'fontname','times','xgrid','on',...
%     'box','on','ygrid','on');
% 
% subplot(212);
% plot(B(2:end),v_1111_a(2:end),'-','linewidth',1);
% hold on
% plot(B(2:end),abs(v_119_a(2:end)),'linewidth',2);
% plot(B(2:end),abs(v_119_b(2:end)),'linewidth',2);
% xlim([0 210]);
% ylim([0 1]);
% xlabel('magnetic field (gauss)');
% set(gca,'fontsize',12);
% legend({str1,str2,str3},'interpreter','latex','fontsize',8);
% ylabel('amplitude');
% set(gca,'fontsize',12,'fontname','times','xgrid','on',...
%     'box','on','ygrid','on');
% 
% % keyboard
% % plot(B,v_99(end,
% text(0.02,.02,'$^{40}\mathrm{K}~4\mathrm{P}_{3/2}$','interpreter','latex','units','normalized',...
%     'verticalalignment','bottom','fontsize',18);
%%
doSave=0;
if doSave
    fprintf('saving figures ...');
    
    % Save to png
    print(hF,'K/K_high_field_CG.png','-dpng','-r400'); 

    disp('done');
end


end

