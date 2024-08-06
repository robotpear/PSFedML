%kmeans 
clear all
load('Caltech101-7.mat')

%100leaves
% label = length(unique(Y))
% 
% for n = 1:3
% 
%    for i = 1:10
%        res = kmeans(X{n}, label,'emptyaction', 'singleton');
%        Gres(i,:) = ClusteringMeasure(Y, res);
%    end
%    fprintf('view %d\n',n)
%    mean(Gres)
% %     std(Gres)
% end

%sonar
% cluster = length(unique(labels))
% 
% for n = 1:3
%     data{n}=data{n}';
% 
%    for i = 1:10
%        res = kmeans(data{n}, cluster,'emptyaction', 'singleton');
%        Gres(i,:) = ClusteringMeasure(labels, res);
%    end
%    fprintf('view %d\n',n)
%    mean(Gres)
% %     std(Gres)
% end

%uci-digts
% cluster = length(unique(truth))
% data{1} = mfeat_kar;
% data{2} = mfeat_fou;
% data{3} = mfeat_fac;
% 
% for n = 1:3
%    for i = 1:10
%        res = kmeans(data{n}, cluster,'emptyaction', 'singleton');
%        Gres(i,:) = ClusteringMeasure(truth, res);
%    end
%    fprintf('view %d\n',n)
%    mean(Gres)
%     std(Gres)
% end

%MSRC
% cluster = length(unique(truth))
% data{1} = msr1;
% data{2} = msr2;
% data{3} = msr3;
% data{4} = msr4;
% data{5} = msr5;
% 
% for n = 1:5
%    for i = 1:10
%        res = kmeans(data{n}, cluster,'emptyaction', 'singleton');
%        Gres(i,:) = ClusteringMeasure(truth, res);
%    end
%    fprintf('view %d\n',n)
%    mean(Gres)
% %     std(Gres)
% end

%Caltech101

% cluster = length(unique(Y))
% 
% for n = 1:6
%    for i = 1:10
%        res = kmeans(X{n}, cluster,'emptyaction', 'singleton');
%        Gres(i,:) = ClusteringMeasure(Y, res);
%    end
%    fprintf('view %d\n',n)
%    mean(Gres)
% %     std(Gres)
% end