clc; clear; % clear，清除
mousetype = '';

%path for C57
if mousetype == 'C57'
   path1 = '';


%%path for Cko 路径
elseif mousetype == 'CKO'
    path1 = ''; 


%%path for ANT
elseif mousetype == 'ANT'
    path1 = ''; 
end

path = path1;
% path = [path11;path12;path13;path14;path15;path16;path17;path18;path19;path20];
% for d = 1: size(path,1)
 
[data, text, ~] = xlsread([path, 'behavior']);% behavior文件读取

%% define zero site || 将食物位置定义成零点；

food_position = data(4, 3:4);
realx= (data(1,3)-data(3,3))*2;  %宽度对应像素数
realy=data(1,4)-data(2,4);       %高度对应像素数

data = data(6:end,[2 3 4]);      %截取轨迹数据及帧数
data(:,2) = (data(:,2) - food_position(1))/realx;
data(:,3) = (data(:,3) - food_position(2))/realy;

data_trial = {};
hand = 'lef';% define forelimb direction


jj = 1;
trial_start = 1;
for ii = 1:(length(data)-1)
    if isnan(data(ii,1))
        trial_end = ii-1;
        data_trial{jj} = data(trial_start:trial_end, :);
        trial_start = trial_end + 2;
        jj = jj + 1;
    end   
end

for ii = 1:length(data_trial)
    if strcmp(hand,'left')
        X_min = min(data_trial{ii}(:,2));
        turn_point = find(data_trial{ii}(:,2) == X_min,1);
        data_trial{ii} = data_trial{ii}(1:turn_point,:);
    else
        X_max = max(data_trial{ii}(:,2));
        turn_point = find(data_trial{ii}(:,2) == X_max,1);
        data_trial{ii} = data_trial{ii}(1:turn_point,:);
    end
end

% figure  
% 
% for ii = 1:length(data_trial)
%    plot(data_trial{ii}(:,2)*-1, data_trial{ii}(:,3) * -1); hold on;
% end
% hold on;
%% hausdorff distance 
dhf = []; D1 = [];D2=[];Dmin1=[];Dmin2=[];a = 1;
for aa = 1:(length (data_trial)-1)
   for bb = aa+1:length (data_trial)
       D2=[];D1 = [];Dmin1=[];Dmin2=[];
       for k= 1:length(data_trial{aa})
           for m=1:length(data_trial{bb})
               dx = data_trial{aa}(k,2)-data_trial{bb}(m,2);
               dy = data_trial{aa}(k,3)-data_trial{bb}(m,3); 
               D1(m) = sqrt(dx.^2 + dy.^2);
           end
           Dmin1(k) = min(D1); % 单点到另外一条线的最短距离
       end
       DMin1 = max(Dmin1);     % 所有点到另外一条线最短距离的最大值
       for j = 1:length(data_trial{bb})
           for n = 1:length(data_trial{aa})
               dx = data_trial{bb}(j,2)-data_trial{aa}(n,2);
               dy = data_trial{bb}(j,3)-data_trial{aa}(n,3); 
               D2(n) = sqrt(dx.^2 + dy.^2);
           end
           Dmin2(j) = min(D2);   
       end
       DMin2 = max(Dmin2);
       dhf(a) = max(DMin1,DMin2);
       a = a + 1;
    end
end
Av_dhf  = mean(dhf,'all');

%%

% end



%%
ShowWin = [-1, 1];
FrameRate = 30;
mousetype = 'C57';

load([path,'matlab.mat'],'neuron');
[num, text] = xlsread([path,'timeStamps1.csv']);
neuron_dff = [];totalBase =[];
p1=[];
Tneuron_dffz = [];  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
neuronC1= neuron.S;
neuronC1(:,:)=0;
neuronC = neuronC1;
n1 = FrameRate./10;
n2= FrameRate./5;

for k = 1: size (neuron.C,1)
    for z =n1:size (neuron.C,2)-10
        if neuron.S(k,z)>0
            N_max = max(neuron.C_raw(k,z-n1:z+10));
            X_max = find(neuron.C_raw(k,z-n1:z+10)==N_max,2);
            X_max = X_max-n1-1;
            neuronC(k,z-n1:z+X_max) = neuron.C_raw (k,z-n1:z+X_max);
        end
    end
end

%%%%%%%%%%%%%%%%%%
    
neuron_dff = neuronC;
% neuron_dff = normalize(neuronC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

success = [];miss = [];nograsp = [];drop = [];grasp=[];fail = [];total=[];
for ii = 1:length(text)
    if strcmp(text{ii,4},'s')
        success = cat(1,success, ii-2);
    elseif strcmp(text{ii,4},'m')
        miss = cat (1,miss, ii-2);
    elseif strcmp(text{ii,4},'n')
        nograsp = cat (1,nograsp, ii-2);
    elseif strcmp(text{ii,4},'d')
        drop = cat (1,drop, ii-2);
    end
end
total = [success;miss;nograsp;drop];
total = sort(total);

totalBase = neuron_dff;totalNeuron1=zeros(size(neuron_dff,1), 0.5*FrameRate+1);
totalNeuron = zeros(size(neuron_dff,1), 0.5*FrameRate+1);totalNeuron2 = [];
for ii = 1:length(total)
    Window(1) = total(ii) - FrameRate*0.1;
    Window(2) = total(ii) + FrameRate*0.4;
    totalNeuron = totalNeuron + neuron_dff(:,Window(1):Window(2));
    totalBase (:,Window(1):Window(2))= -1; %将所有实验窗中的数据归零
    totalNeuron1 = neuron_dff(:,Window(1):Window(2));
    if ii ==1
        totalNeuron2 = totalNeuron1;
    else
        totalNeuron2 = [totalNeuron2 totalNeuron1];
    end
end
totalNeuron = totalNeuron./length(total);  %平均值
totalBase (:,all(totalBase<-0.5,1)) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%去掉不反应的细胞
for ii = 1: size(totalBase,1) 
    p1(ii,:) = ranksum(totalBase(ii,:),totalNeuron(ii,:));% 判断运动时间内是否与背景显著差异；ranksum test；
end

aa =1; bb = 1;Tneuron_dff1=[];aver_total=[];aver_Base=[];
for jj = 1:size(totalNeuron,1)
    aver_total (jj,:) = mean(totalNeuron(jj,:));
    aver_Base(jj,:) = nanmean (totalBase(jj,:));
        if  p1(jj,:)<0.05 & aver_total (jj,:)> aver_Base(jj,:)%同时满足显著差异且高于背景；
          Tneuron_dff1(bb,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞
          bb = bb + 1;
        end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tneuron = mean(Tneuron_dff1,1);totalNeurons=[];
aa = 1;
for ii = 1:length(total)
    Window(1) = total(ii) - FrameRate*0.1;
    Window(2) = total(ii) + FrameRate*0.4; 
%     if mean(Tneuron(:,Window(1):Window(2)),'all')~= 0
       totalNeurons(aa,:) = Tneuron(:,Window(1):Window(2));
       aa = aa+1;
%     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% population neurons between


total_corr_trials =[];a = 1;
for ii = 1: size(totalNeurons,1)-1
    for jj =ii+1:size(totalNeurons,1)
    total_corr_trialm = corrcoef(totalNeurons(ii,:),totalNeurons(jj,:));
    total_corr_trials(a) = total_corr_trialm (1,2); 
    a = a + 1;
    end
end

Mean_total_corr  = mean (total_corr_trials,'all');

dhf_corr(:,1) = dhf(1,:);
dhf_corr(:,2) = total_corr_trials(1,:);
dhf_corr2 = dhf_corr(all(~isnan(dhf_corr),2),:); 
dhf_corr3 = dhf_corr2(all(dhf_corr2(:,2) > 0, 2),:);
figure 
scatter (dhf_corr2(:,1),dhf_corr2(:,2));
figure
scatter (dhf_corr3(:,1),dhf_corr3(:,2));


