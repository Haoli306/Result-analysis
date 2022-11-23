clc;%close all;
clear;clear all; clear global;
path = 'D:\work\tracking\';
filename = 'KO-B1.xlsx';
% filename = 'c57-x3.xlsx';
sheet1= 'ko-b1-t1'; sheet2= 'ko-b1-t3'; sheet3= 'ko-b1-t6';
% sheet1 = 'x3-c57-d1';sheet2 = 'x3-c57-d3';sheet3 = 'x3-c57-d6';
sheet = [sheet1;sheet2;sheet3];
for d = 1: size(sheet,1)
[data{d}, text, ~] = xlsread([path, filename], sheet(d,:));
food_position = data{d}(4, 3:4);
realx= data{d}(1,11);
realy=data{d}(1,12);
data{d} = data{d}(6:end,2:4);
text = text(7:end, 6);
data{d}(:,2) = (data{d}(:,2) - food_position(1))/realx;
data{d}(:,3) = (data{d}(:,3) - food_position(2))/realy;
data_trial = {};
hand = 'left';
end
figure 
plot1 = plot(1*data{1}(:,2),-1*data{1}(:,3),'linewidth',1,'color','b');
axis([-3 4 -1.5 2]);plot1.Color(4)=0.2;
hold on;
plot2 = plot(1*data{2}(:,2)+1.5,-1*data{2}(:,3),'linewidth',1,'color','r');plot2.Color(4)=0.2;
hold on ;
plot3 = plot(1*data{3}(:,2)+3,-1*data{3}(:,3),'linewidth',1,'color','m'); plot3.Color(4)=0.2;

hold off;