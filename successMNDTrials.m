clc;clear;clear all; clear global;
ShowWin = [-1, 1];
FrameRate = 30;
mousetype = '';

%%%%%%%%%%% C57 %%%%%%%%%%%%%
if mousetype == 'C57'
filename1= '';

filename = [filename1];

%%%%%%%%%%%%%%% CKO %%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif mousetype == 'CKO'
    
filename1= '';

filename = [filename1];

elseif mousetype =='ANT'
    
filename1 = '';

filename = [filename1];
end


for cc = 1:size(filename,1)
load([filename(cc,:),'\matlab.mat'],'neuron');
[num, text] = xlsread([filename(cc,:),'\timeStamps1.csv']);
neuron_dff = [];totalBase =[];
totalResp=[];graspResp=[];failResp=[];totalResp2=[];graspResp=[];failResp=[];
p1=[];p2=[];p3=[];

%%%%%%%%%%%%%%%%%%%%%%%%%
neuronC1= neuron.S;
neuronC1(:,:)=0;
neuronC = neuronC1;
n1 = FrameRate./10;
n2= FrameRate./5;

for k = 1: size (neuron.C,1)
    for z =1:size (neuron.C,2)-10
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

success = [];miss = [];nograsp = [];drop = [];grasp=[];fail=[];
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
grasp =[success;drop];
fail= [miss;nograsp];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalBase = neuron_dff;totalNeuron =[];totalNeuron1=zeros(size(neuron_dff,1), 0.5*FrameRate+1);
totalNeuron = zeros(size(neuron_dff,1), 0.5*FrameRate+1);totalNeuron2 = [];
for ii = 1:length(total)
    Window(1) = total(ii) - FrameRate*0.1;
    Window(2) = total(ii) + FrameRate*0.4;
    totalNeuron = totalNeuron + neuron_dff(:,Window(1):Window(2));
    totalBase (:,Window(1):Window(2))= -1; %将所有实验窗中的数据归-1
    totalNeuron1 = neuron_dff(:,Window(1):Window(2));
    if ii ==1
        totalNeuron2 = totalNeuron1;
    else
        totalNeuron2 = [totalNeuron2 totalNeuron1];
    end
end

totalNeuron = totalNeuron./length(total);  %平均值
totalBase (:,all(totalBase<-0.5,1)) = [];




% Base = totalBase(:,20:100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%total
for ii = 1: size(totalBase,1) 
    p1(ii,:) = ranksum(totalBase(ii,:),totalNeuron2(ii,:));% 判断运动时间内是否与背景显著差异；ranksum test；
end

bb = 1;Tneuron_dff=[];aver_total=[];aver_Base=[];SD=[];
for jj = 1:size(totalNeuron,1)
    SD = std(neuron_dff(jj,:));
    aver_total (jj,:) = mean(totalNeuron(jj,:));
    aver_Base(jj,:) = mean (totalBase(jj,:));
    s(jj,:) = totalNeuron(jj,:)-aver_Base(jj,:);
        if p1(jj)<0.05 & aver_total (jj,:)> aver_Base(jj,:) %同时满足显著差异且高于背景；
          Tneuron_dff(bb,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞 
          bb =bb+ 1;
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tneuron = mean (Tneuron_dff,1); totalResp1 = [];
for ii = 1:length(total)
    Window = total(ii) + ShowWin.*30;
    totalResp1(ii,:) = Tneuron(:,Window(1):Window(2));
end
[~, peak_time] = max(totalResp1(:,:),[],2);
[~, I_t1] = sort(peak_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

totalResp = zeros(size(Tneuron_dff,1), (ShowWin(2)-ShowWin(1)).*30+1);%筛选总的细胞
for ii = 1:length(total)
    Window = total(ii) + ShowWin.*30;
    totalResp = totalResp + Tneuron_dff(:,Window(1):Window(2)); 
end
 totalResp = totalResp./length(total);  %平均值

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Neu_cor =[];
for ii =1 : size(totalResp,1)-1
    for jj = ii+1:size(totalResp,1)
        Neu_corm = corrcoef (totalResp(ii,:),totalResp(jj,:));
        Neu_cor(ii,jj) = Neu_corm(1,2);
    end
end
Neu_cor_Mean(cc) = 2 * mean(Neu_cor,'all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~, peak_time1] = max(totalResp(:,:),[],2);
[~, I_t] = sort(peak_time1);

% bb = 1;Tneuron_dff2=[];
% for jj = 1:size(totalNeuron,1)
%     aver_total (jj,:) = mean(totalNeuron(jj,:));
%     aver_Base(jj,:) = mean (totalBase(jj,:));
%         if aver_total(jj,:)<aver_Base(jj,:) & p1(jj,:)<0.05 %同时满足显著差异且低于背景；
%           Tneuron_dff2(bb,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞 
%           bb =bb+ 1;
%         end
% end

% totalResp2 = zeros(size(Tneuron_dff2,1), (ShowWin(2)-ShowWin(1)).*30+1);%筛选总的细胞
% for ii = 1:length(total)
%     Window = total(ii) + ShowWin.*30;
%     totalResp2 = totalResp2 + Tneuron_dff2(:,Window(1):Window(2)); 
% end
%  totalResp2 = totalResp2./length(total);  %平均值
% 
% 
% [~, peak_time] = min(totalResp2(:,:),[],2);
% [~, I_t2] = sort(peak_time);
% %%%%%%%%%%%%%%%%%%grasp
% grasp_trial=[];
% grasp_trial = zeros(size(neuron_dff,1), 0.5*FrameRate+1);
% for ii = 1:length(grasp)
%      Window(1) = grasp(ii) - FrameRate*0.1;
%     Window(2) = grasp(ii) + FrameRate*0.4;
%     grasp_trial = grasp_trial + neuron_dff(:,Window(1):Window(2));
% end
% grasp_trial = grasp_trial./length(grasp);
% 
% for ii = 1: size(totalBase,1) 
%     p2(ii,:) = ranksum(totalBase(ii,:),grasp_trial(ii,:));% 判断运动时间内是否与背景显著差异；ranksum test；
% end
% 
% bb = 1;aa = 1;aver_grasp = [];Gneuron_dff=[];Gneuron_dff2=[];
% for jj = 1:size(totalNeuron,1)
%     aver_grasp (jj,:) = mean(grasp_trial(jj,:));
%     aver_Base(jj,:) = mean (totalBase(jj,:));
%         if  p2(jj,:)<0.05  & aver_grasp (jj,:)> aver_Base(jj,:) %同时满足显著差异且高于背景；
%           Gneuron_dff(bb,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞 
%           averG(bb) = aver_grasp(jj);
%           bb =bb+ 1;
% %         elseif aver_grasp(jj,:)<aver_Base(jj,:) & p2(jj,:)<0.05 %同时满足显著差异且高于背景；
% %           Gneuron_dff2(aa,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞 
% %           aa =aa+ 1;
%         end
% end
% 
% graspResp = zeros(size(Gneuron_dff,1), (ShowWin(2)-ShowWin(1)).*30+1);
% for ii =1:length(grasp)
%     Window = grasp(ii)+ ShowWin.*30;
%     graspResp = graspResp + Gneuron_dff(:,Window(1):Window(2));
% end
% graspResp = graspResp./length(grasp);
% 
% [~, peak_time2] = max(graspResp(:,:),[],2);
% [~, I_g] = sort(peak_time2);


% graspResp2 = zeros(size(Gneuron_dff2,1), (ShowWin(2)-ShowWin(1)).*30+1);%筛选总的细胞
% for ii = 1:length(grasp)
%     Window = grasp(ii) + ShowWin.*30;
%     graspResp2 = graspResp2 + Gneuron_dff2(:,Window(1):Window(2)); 
% end
%  graspResp2 = graspResp2./length(grasp);  %平均值
% 
% 
% [~, peak_time] = min(graspResp2(:,:),[],2);
% [~, I_g2] = sort(peak_time);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% miss
% fail_trial=[];
% fail_trial = zeros(size(neuron_dff,1), 0.5*FrameRate+1);
% for ii = 1:length(fail)
%      Window(1) = fail(ii) - FrameRate*0.1;
%     Window(2) = fail(ii) + FrameRate*0.4;
%     fail_trial = fail_trial + neuron_dff(:,Window(1):Window(2));
% end
% fail_trial = fail_trial./length(fail);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%去掉不反应的细胞
% 
% for ii = 1: size(totalBase,1) 
%     p3(ii,:) = ranksum(totalBase(ii,:),fail_trial(ii,:));% 判断运动时间内是否与背景显著差异；ranksum test；
% end
% 
% bb = 1;aa = 1;aver_fail=[];Mneuron_dff=[];Mneuron_dff2=[];
% for jj = 1:size(totalNeuron,1)
%     aver_fail (jj,:) = mean(fail_trial(jj,:));
%     aver_Base(jj,:) = mean (totalBase(jj,:));
%         if p3(jj,:)<0.05  & aver_fail (jj,:)> aver_Base(jj,:)%同时满足显著差异且高于背景；
%           Mneuron_dff(bb,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞 
%           bb =bb+ 1;
% %         elseif aver_fail(jj,:)<aver_Base(jj,:) & p3(jj,:)<0.05 %inhibited neurons
% %           Mneuron_dff2(aa,:) = neuron_dff(jj,:);  
% %           aa= aa+1;
%         end
% end
% 
% failResp = zeros(size(Mneuron_dff,1), (ShowWin(2)-ShowWin(1)).*30+1);
% for ii = 1:length(fail)
%     Window = fail(ii) + ShowWin.*30;
%     failResp = failResp + Mneuron_dff(:,Window(1):Window(2));
% end
% failResp = failResp./length(fail);
% 
% [~, peak_time3] = max(failResp(:,:),[],2);
% [~, I_f] = sort(peak_time3);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bb = 1;Tneu=[];
% for jj = 1:size(totalNeuron,1)
%     aver_fail (jj,:) = mean(fail_trial(jj,:));
%     aver_Base(jj,:) = mean(totalBase(jj,:));
%     aver_grasp (jj,:) = mean(grasp_trial(jj,:));
%         if (p2(jj,:)<0.05 & aver_grasp (jj,:)> aver_Base(jj,:)) || ...
%                 (p3(jj,:)<0.05 & aver_fail (jj,:)> aver_Base(jj,:)) %同时满足显著差异且高于背景；
%           Tneu(bb,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞 
%           bb = bb + 1;
%         end
% end

% TneuResp = [];
% TneuResp = zeros(size(Tneu,1), (ShowWin(2)-ShowWin(1)).*30+1);
% for ii = 1:length(fail)
%     Window = total(ii) + ShowWin.*30;
%     TneuResp = TneuResp + Tneu(:,Window(1):Window(2));
% end
% TneuResp = TneuResp./length(total);
% 
% [~, peak_time4] = max(TneuResp(:,:),[],2);
% [~, I_tT] = sort(peak_time4);

% failResp2 = zeros(size(Mneuron_dff2,1), (ShowWin(2)-ShowWin(1)).*30+1);
% if length(fail)>10
%     g = randperm(length(fail),10);
%   for ii = 1:10
%     Window = fail(g(ii)) + ShowWin.*30;
%     failResp2 = failResp2 + Mneuron_dff2(:,Window(1):Window(2));
%   end
%  failResp2 = failResp2./10;
% else
%     for ii =1:length(fail)
%         Window = fail(ii)+ ShowWin.*30;
%          failResp2 = failResp2 + Mneuron_dff2(:,Window(1):Window(2));
%     end
%     failResp2 = failResp2./length(fail);
% end
% failResp2 = zeros(size(Mneuron_dff2,1), (ShowWin(2)-ShowWin(1)).*30+1);
% for ii = 1:length(fail)
%     Window = fail(ii) + ShowWin.*30;
%     failResp2 = failResp2 + Mneuron_dff2(:,Window(1):Window(2));
% end
% failResp2 = failResp2./length(fail);
% 
% [~, peak_time] = min(failResp2(:,:),[],2);
% [~, I_f2] = sort(peak_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bb = 1;aa = 1;aver_fail=[];NTneuron_dff=[];NTneuron_dff2=[];
% for jj = 1:size(totalNeuron,1)
%     aver_fail (jj,:) = mean(fail_trial(jj,:));
%     aver_Base(jj,:) = mean (totalBase(jj,:));
%     aver_grasp (jj,:) = mean(grasp_trial(jj,:));
%         if (aver_fail(jj,:)>aver_Base(jj,:)||aver_grasp (jj,:)>aver_Base(jj,:)) & ...
%                 (p3(jj,:)<0.05||p2(jj,:)<0.05) %同时满足显著差异且高于背景；
%           NTneuron_dff(bb,:) = neuron_dff(jj,:);% new neuron_dff 去掉不反应的细胞 
%           bb =bb+ 1;
%         elseif (aver_fail(jj,:)<aver_Base(jj,:)||aver_grasp (jj,:)<aver_Base(jj,:)) & ... 
%                 (p3(jj,:)<0.05||p2(jj,:)<0.05) %inhibited neurons
%           NTneuron_dff2(aa,:) = neuron_dff(jj,:);  
%           aa= aa+1;
%         end
% end
% 
% NtotalResp = zeros(size(NTneuron_dff,1), (ShowWin(2)-ShowWin(1)).*30+1);
% for ii = 1:length(total)
%     Window = total(ii) + ShowWin.*30;
%     NtotalResp = NtotalResp + NTneuron_dff(:,Window(1):Window(2));
% end
% NtotalResp = NtotalResp./length(total);
% 
% [~, peak_time] = max(NtotalResp(:,:),[],2);
% [~, I_NT] = sort(peak_time);
% 
% % NtotalResp2 = zeros(size(NTneuron_dff2,1), (ShowWin(2)-ShowWin(1)).*30+1);
% % for ii = 1:length(total)
% %     Window = total(ii) + ShowWin.*30;
% %     NtotalResp2 = NtotalResp2 + NTneuron_dff2(:,Window(1):Window(2));
% % end
% % NtotalResp2 = NtotalResp2./length(total);
% % 
% % [~, peak_time] = max(NtotalResp2(:,:),[],2);
% % [~, I_NT2] = sort(peak_time);
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AllNeuron_total1=[]; 
AllNeuron_total1 = mean(totalResp);
   



X = [-1*FrameRate:1*FrameRate];% X axis
X = X./FrameRate;% 转化为时间
Y1 = [1:size(graspResp,1)];%Y axis
Y2 = [1:size(failResp,1)];%Y axis
Y3 = [1:size(totalResp1,1)];%Y axis
% Y4 = [1:size(graspResp2,1)];%Y axis
% Y5 = [1:size(failResp2,1)];%Y axis
% Y6 = [1:size(totalResp2,1)];%Y axis


figure
subplot(1,3,1)
imagesc(X,Y1,graspResp(I_g,:),[0,0.1]);xlabel('Time/s');ylabel('Neurons');
subplot(1,3,2)
imagesc(X,Y2,failResp(I_f,:),[0,0.1]);xlabel('Time/s');ylabel('Neurons');
subplot(1,3,3)
imagesc(X,Y3,totalResp1(I_t1,:),[0,0.1]);xlabel('Time/s');ylabel('Neurons');
% subplot(2,3,4)
% imagesc(X,Y4,graspResp2(I_g2,:),[-0.02,0.05]);xlabel('Time/s');ylabel('Neurons');
% subplot(2,3,5)
% imagesc(X,Y5,failResp2(I_f2,:),[-0.02,0.05]);xlabel('Time/s');ylabel('Neurons');
% subplot(2,3,6)
% imagesc(X,Y6,totalResp2(I_t2,:),[-0.02,0.05]);xlabel('Time/s');ylabel('Neurons');
A(cc)=(a-1)/(a+c-2);
B(cc)=(b-1)/(b+c-2);

cd(filename(cc,:));
save data3.mat;
end




