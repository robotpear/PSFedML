clear all
load('sonar.mat')

para = [];
para.rho = 0.2;
% beta = [0.05,0.1,0.15,0.2,0.25];
% lambda = [0.05,0.1,0.15,0.2,0.25];
para.Cum = 20;
para.max = 100;
para.maxc = 1;
para.maxs = 1;
para.var = 0.0;
rng(9)

%% sonar
fprintf('load dataset sonar');
para.beta = 1;
para.lambda = 0.1;
[Hp, Vres, V, Err, Errs, Bp,Gp] = PSFedML(data, labels,para);

% for i=1:5
%     rng(6)
%     para.beta = beta(i);
%     for j=1:5
%         para.lambda = lambda(j);
%         [Hp, Vres, V,Errs,Bp,Gp] = PSFedML(data, labels,para);
%         Vs = [beta(i);lambda(j);Vres(1);Vres(2);Vres(3)];
%         Vac(i,j)=Vres(1);
%         Vp(i,j)=Vres(2);
%         Vnmi(i,j)=Vres(3);
%     end
% end
% figure;
% set(gcf,'unit','normalized','position',[0.1,0.25,0.8,0.5])
% subplot(1,3,1)
% bar3(Vac)
% set(gca,'xticklabel',[0.001,0.1,0.2,0.5,1])
% set(gca,'yticklabel',[0.001,0.1,0.2,0.5,1])
% xlabel('Value of \lambda','Rotation',22);
% ylabel('Value of \beta','Rotation',-35);
% zlabel('ACC');
% 
% 
% subplot(1,3,2)
% bar3(Vp)
% set(gca,'xticklabel',[0.001,0.1,0.2,0.5,1])
% set(gca,'yticklabel',[0.001,0.1,0.2,0.5,1])
% xlabel('Value of \lambda','Rotation',22);
% ylabel('Value of \beta','Rotation',-35);
% zlabel('Purity');
% 
% subplot(1,3,3)
% bar3(Vnmi)
% set(gca,'xticklabel',[0.001,0.1,0.2,0.5,1])
% set(gca,'yticklabel',[0.001,0.1,0.2,0.5,1])
% xlabel('Value of \lambda','Rotation',22);
% ylabel('Value of \beta','Rotation',-35);
% zlabel('NMI');