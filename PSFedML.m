%
function [ Hp, Vres, V, Err, Errs, Bp ,Gp] = PSFedML(Data, labels,para)
%
% Partially Shared Federated Multiview Learning (PSFedML)



%% initial 
% Cnum = length(unique(labels));
Cnum = para.Cum;
rho = para.rho;
maxIter = para.max;
maxIterc = para.maxc;
maxIters = para.maxs;
thresh = 1e-4;
beta = para.beta;
lambda = para.lambda;
var = para.var;



%num of views 
nView = length(Data);
%num of samples(assumption: all clients hold same num of samples)

%num of features
nFea = zeros(nView,1);
nSamp = size(Data{1}, 2);
sumf = 0;
for n = 1:nView
    nFea(n) = size(Data{n}, 1);
    sumf = nFea(n)+sumf;
end
%% clients initialize matrix B^(p), H^(p)

Xp = cell(nView,1);
Bp = cell(nView,1);
Hp = cell(nView,1);
for n = 1:nView
%     for  i = 1:nSamp
%         Data{n}(:,i) = (Data{n}(:,i) - mean(Data{n}(:,i)))/std(Data{n}(:,i));
%     end

%     FlattenedData = Data{n}(:); 
%     
%     MappedFlattened = mapminmax(FlattenedData, 0, 1); 
%     
%     Xp{n} = reshape(MappedFlattened, size(Data{n})); 

    Xp{n}=scaleSVM(Data{n},0,1);
%     Xp{n}=Data{n};
    

    Bp{n} = rand(nFea(n),Cnum);
    Hp{n} = rand(Cnum,nSamp);
end

%% initial graph regular
for n=1:nView
    fea{n} = Xp{n}';
    options = [];
    options.Metric = 'Euclidean';
    options.NeighborMode = 'KNN';
    options.k = 10;
    options.WeightMode = 'HeatKernel';
    options.t = 1;
    W{n} = constructW(fea{n},options);
    D{n} = diag(sum(W{n},2));
end






%% serve initialize mateix H, G, V

K = rho*Cnum;
Gp = cell(nView,1);
V = rand(nSamp,K);
for n=1:nView
    Gp{n} = rand(Cnum,K);
end


%%  update

for iter = 1:maxIter


%     if mod(iter, 1) == 0
% %        fprintf('numOfOutliers = %d, ratio = %f\n', length(Idx),ratio);
% %        fprintf('%dth iteration, obj = %f \n', it, obj);
%        fprintf('processing iteration %d...\n', iter);
%     end


%% client update   
    for m = 1:maxIterc
    % update B^{n},loop over clients (in parallel)
        
        for n = 1:nView 
            Hpt = Hp{n}';      
            Bp{n} = Bp{n}.*((Xp{n}*Hpt)./max((Bp{n}*Hp{n}*Hpt),1e-5));
        end
        clear Hpt
    
    
    % update Hp^{n}
        for n = 1:nView
            Bpt = Bp{n}'; 
            Hp{n} = Hp{n}.*((Bpt*Xp{n}+lambda*Hp{n}*W{n})./max((Bpt*Bp{n}*Hp{n}+lambda*Hp{n}*D{n}),1e-5));
        end
        clear Bpt
    


        if (m > 0) 
            loss = 0; 
            for n = 1:nView
    %             loss1 = sum(sum(sum((Xp{n} - Bp{n}*Hp{n}).^2)));
    %             loss = sum(loss1);
                loss1 = (norm((Xp{n}-Bp{n}*Hp{n}),'fro'))^2+trace(Hp{n}*(D{n}-W{n})*Hp{n}');
                loss = loss1 + loss;
                
            end
            loss = loss/nView;
            Err(m+(iter-1)*maxIterc) = loss;
            clear loss1
        end

    % convergence or not
        if(m > 1)
            diffc = abs(Err(m-1) - Err(m));
            if(diffc < thresh)
                Err(m+1:maxIterc)=Err(m);
                break;
            end
        end
    

    end


a=randsrc(nView,1,[0 1; var 1-var]);
M = [1:nView];
nViewG = a'.*M;
nViewG = find(nViewG~=0);



 %% server update 

    for m = 1:maxIters
            

            
        % update G
        for n = nViewG
            dp{n} = 0.5./sqrt(sum(Gp{n}.*Gp{n},2)+1e-5);
            Dp{n} = diag(dp{n});
            Gp{n}=inv((Hp{n}*Hp{n}'+beta*Dp{n}))*Hp{n}*V;
        end

%         % sort Gp
        for n = nViewG
            for i = 1:Cnum
               gp{n}(i) = norm(Gp{n}(i,:),2);  
            end
            [~,Gp_index{n}] = sort(gp{n},'descend');  
            Gp_index{n} = Gp_index{n}(1:K);  
            Gp_id{n} = zeros(Cnum,K);  
            for i = 1:K
                Gp_id{n}(Gp_index{n}(i),i)=1;
            end
            Gp{n} = Gp_id{n};
        end
    

        %% update V

          sumhg =zeros(nSamp,K);
          for n =nViewG
              avgw(n) = nFea(n)/sumf; 
              sumhg = avgw(n)*Hp{n}'*Gp{n}+sumhg;
          end
          V = sumhg;

        
%         % sort Gp
%         for n = 1:nView
%             for i = 1:Cnum
%                gp{n}(i) = norm(Gp{n}(i,:),2);  
%             end
%             [~,Gp_index{n}] = sort(gp{n},'descend');  
%             Gp_index{n} = Gp_index{n}(1:K);  
%                     Gp{n}=Gp_index{n};
%         end

        
    
        
        % server loss function value

        if (m > 0)
            loss = 0;
            for n = nViewG
                loss = (norm((Hp{n}'*Gp{n}-V),'fro')^2 + beta*trace(Gp{n}'*Dp{n}*Gp{n}))+loss; 
%                 loss = norm((Hp{n}'*Gp{n}-V),'fro')^2 +loss;
            end
    %         loss = sum(sum((hg-V).^2));
            Errs(m+(iter-1)*maxIters) = loss;
%             ERRS(iter,m) = Errs(m);
%             plot(Errs)
%             hold on 
        end



        
    
%         % convergence or not
        if(m > 1)
            diffs = abs(Errs(m-1) - Errs(m));
            if(diffs < thresh)
                Errs(m:maxIters)=Errs(m);
                break;
            end
        end


    end
    
    
    %% transform V to clients
    for n= nViewG
        Hpt = Hp{n}';
        [row,col]=find(Gp_id{n}~=0);
        for i =1:length(row)
            Hpt(:,row(i))=V(:,col(i));
        end
        Hp{n}=Hpt';
%         Hp{n}=(Hp{n}'*Gp_id{n})';
    end



end
%% Clustering and plot


%     res1 = kmeans(Hp{2}', length(unique(labels)),'emptyaction', 'singleton');
%     Hpres = ClusteringMeasure(labels, res1);
%     HAcc = Hpres(1)
%     HPurity = Hpres(2)
%     HNMI = Hpres(3)

    %
    [res2,C] = kmeans(V, length(unique(labels)),'emptyaction', 'singleton');
    Vres = ClusteringMeasure(labels, res2);
    VAcc = Vres(1)
    VPurity = Vres(2)
    VNMI = Vres(3)
    Gp_index;

%     figure
%     scatter(Data{1}(:,1),'LineWidth',10);
%     scatter(V,res2)
%     gscatter(V(:,1),V(:,2),res2,'bgmrk')
%     hold on
%     plot(C(:,1),C(:,2),'kx')
%     legend('C1','C2','C3','C4','Centroid')


%% loss value plot
%     figure('Name','client loss :Err')
%     plot(Err)
%     xlabel('iter');

    figure('Name','server loss :Errs')
    plot(Errs,'-bs',...
    'LineWidth',1,...
    'MarkerIndices',[1:10:100],...
    'MarkerSize',6,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5,0.5,0.5])
    ylabel('Objective Function Value')
    xlabel('Iterations');
    grid on


%     
    %

end
%