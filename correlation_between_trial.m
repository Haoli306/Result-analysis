clc;clear;clear all; clear global;
ShowWin = [-1, 1];
FrameRate = 30;
mousetype = 'C57';

%%%%%%%%%%% C57 %%%%%%%%%%%%%
if mousetype == 'C57'
filename1= '';

filename = [filename1];


%%%%%%%%%%%%%%% CKO %%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif mousetype == 'CKO'
    
filename1= '';

filename = [filename1];

elseif mousetype == 'ANT';
filename1 = '';


filename = [filename1];
end

for cc = 1:size(filename,1)
load([filename(cc,:),'\matlab.mat'],'neuron');
[num, text] = xlsread([filename(cc,:),'\timeStamps1.csv']);
neuron_dff = [];totalBase =[];
p1=[];p2=[];p3=[];
Tneuron_dffz = []; Gneuron_dffz = []; Fneuron_dffz = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% neuron_dff = normalize(neuronC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% for ii = 1:size(neuron.C,1)
%     neuron_dff(ii,:) = (neuron.C_raw(ii,:) - mean(neuron.C_raw(ii,:)))/std(neuron.C_raw(ii,:));% z-score of each neuron
% end
success = [];miss = [];nograsp = [];drop = [];grasp=[];fail = [];
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

totalBase = neuron_dff;totalNeuron1=zeros(size(neuron_dff,1), 0.5*FrameRate+1);
totalNeuron = zeros(size(neuron_dff,1), 0.5*FrameRate+1);totalNeuron2 = [];
for ii = 1:length(total)
    Window(1) = total(ii) - FrameRate*0.1;
    Window(2) = total(ii) + FrameRate*0.4;
    totalNeuron = totalNeuron + neuron_dff(:,Window(1):Window(2));
    totalBase (:,Window(1):Window(2))= -1; 
    totalNeuron1 = neuron_dff(:,Window(1):Window(2));
    if ii ==1
        totalNeuron2 = totalNeuron1;
    else
        totalNeuron2 = [totalNeuron2 totalNeuron1];
    end
end
totalNeuron = totalNeuron./length(total);  
totalBase (:,all(totalBase<-0.5,1)) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1: size(totalBase,1) 
    p1(ii,:) = ranksum(totalBase(ii,:),totalNeuron(ii,:));% ranksum test；
end

aa =1; bb = 1;Tneuron_dff1=[];Tneuron_dff2 =[];Tneuron_dff=[];
for jj = 1:size(totalNeuron,1)
    aver_total (jj,:) = mean(totalNeuron(jj,:));
    aver_Base(jj,:) = nanmean (totalBase(jj,:));
        if  p1(jj,:)<0.05 & aver_total (jj,:)> aver_Base(jj,:)
          Tneuron_dff1(bb,:) = neuron_dff(jj,:); 
          Tneuron_dff(bb,:) = neuron.S(jj,:);
          bb = bb + 1;
        end 
end
act_rate(cc) = size(Tneuron_dff1,1)/size(totalNeuron,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tneuron = mean(Tneuron_dff,1);totalNeurons=[];Timing = [];
aa = 1;Tneuroncon = []; totalcon = [];bb = 1; d =1; totalReac = [];
for ii = 1:length(total)
     Window(1) = total(ii) - FrameRate*0.1;
     Window(2) = total(ii) + FrameRate*0.4; 
     Window(3) = total(ii) - FrameRate*0.6;
     Window(4) = total(ii) - FrameRate*0.1;
     totalReac (ii,:) = mean(Tneuron(:,Window(1):Window(2)) - Tneuron(:,Window(3):Window(4))); 
     totalBs (ii,:) = Tneuron(:,Window(3):Window(4)); %%baseline的movement  
    
     if mean(Tneuron(:,Window(1):Window(2)),'all')~= 0
        totalNeurons(aa,:) = Tneuron(:,Window(1):Window(2));
        for k = 1:size(totalNeurons,2)
            if totalNeurons(aa,k)~=0
                Timing(d) =(k-1)/FrameRate;
                d = d+1;
            end
        end     
        aa = aa+1;
     end
     reoccurance1(cc) = (aa-1)/length(total);
     Tneuroncon = Tneuron_dff(:,Window(1):Window(2));
     m = mean(Tneuroncon,'all');
     if m ~= 0
       totalcon(:,bb) = mean (Tneuroncon,2);
       bb = bb + 1;
     end
end
TimingSD(cc) = std(Timing);
%%%求和
totalBsSum = sum(totalBs);
totalBsm = totalBs(:)';
totalSD(cc) = std(totalBsm);
totalReac_mean(cc) = mean(totalReac);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%population neurons between trials
totalcon_corr = [];
for ii = 1:size(totalcon,2)
    for jj = (ii + 1):size(totalcon,2)
        totalcon_corrm = corrcoef(totalcon(:,ii),totalcon(:,jj));
        totalcon_corr (ii,jj) = totalcon_corrm(1,2);
    end
end
Meantotalcon_corr(cc) = nanmean(totalcon_corr,'all').*2;

total_corr_trials =[];ptotal_dis=[];
for ii = 1: size(totalNeurons,1)-1
    for jj =ii+1:size(totalNeurons,1)
    total_corr_trialm = corrcoef(totalNeurons(ii,:),totalNeurons(jj,:));
    total_corr_trials(ii,jj) = total_corr_trialm (1,2); 
    
    end
end
Trials{cc} = totalNeurons;
Mean_total_corr (cc) = mean (total_corr_trials,'all').*2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% totalResp = zeros(size(Tneuron_dff,1), FrameRate*2+1);%筛选总的细胞
% for ii = 1:length(total)
%     Window(1) = total(ii) - FrameRate;
%     Window(2) = total(ii) + FrameRate;
%     totalResp = totalResp + Tneuron_dff(:,Window(1):Window(2)); 
% end
%  totalResp = totalResp./length(total);  %平均值
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timings=[]; t = 1;Ws=[];
%  for ii = 1: size(totalResp,1)
%      for jj = 1:size(totalResp,2)
%          if totalResp(ii,jj) ~= 0
%             Timings(t) = (jj-1)/FrameRate;
%             t = t+1;
%          end
%      end
%  end
% TimingSDs(cc) = std(Timings);
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%correlation between trials
% total_aver_corr = [];total_corr_trialm=[];total_averNeuron=[];total_corr_trial=[];
% total_trial = zeros(size(neuron_dff,1), 0.5*FrameRate+1);%%计算整体不同trial的神经元平均活性
% total_dism=[];total_disS=[];Sg_reoccurance=[];TimeSDm=[];m=1;
% for jj = 1: size (Tneuron_dff,1) 
%     d = 1;h=1;Time =[];
%    for ii = 1:length(total)
%     Window(1) = total(ii) - FrameRate;
%     Window(2) = total(ii) + FrameRate;
%       if any(Tneuron_dff(jj,Window(1):Window(2))>0)
%          total_averNeuron(d,:) = Tneuron_dff(jj,Window(1):Window(2));
%          d =d + 1;
%       end
%    end
%    for k = 1: size (total_averNeuron,1)
%      for z = 1:size (total_averNeuron,2)
%          if total_averNeuron(k,z) ~= 0  % 单神经元calcium event对应的时间
%             Time(h) = (z-1)./FrameRate;
%             h = h + 1;
%          end
%      end
%    end
%    if length(Time) > 2
%      TimeSDm(m) = std (Time);
%      m = m+1;
%    end
%   Sg_reoccurance(jj) = (d-1)/length(total);% 计算单个神经元激活率
%   for aa =1:size(total_averNeuron,1)
%      for bb = aa+1:size(total_averNeuron,1)
%          total_corr_trialm(:,:) = corrcoef(total_averNeuron(aa,:),total_averNeuron(bb,:)) ;
%          total_corr_trial(aa,bb) = total_corr_trialm(1,2);
%          total_dism(aa,bb) = (sum((total_averNeuron(aa,:)-total_averNeuron(bb,:)).^2)).^0.5; 
%      end
%   end
%    total_disS(jj,:)= 2.*nanmean(total_dism,'all');
%    total_aver_corr (jj,:) = 2.*nanmean(total_corr_trial,'all');
% end
% TimeSD(cc) = nanmean(TimeSDm);
%  
% 
% total_aver_corr(isnan(total_aver_corr))=[];
% total_aver_corrN1 = mean(total_aver_corr);
% total_disSMean = nanmean(total_disS);
% reoccurance2(cc) = mean(Sg_reoccurance);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(filename(cc,:));
save corr1.mat;
end



