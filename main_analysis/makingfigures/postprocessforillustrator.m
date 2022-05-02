%% fig_1_3_licks:
set(gcf,'Color','w');
fig=gcf;
% fig.Units = 'inches';
% fig.Position = [0 0 6 2.4];
fig.Children.YLabel.String={};
fig.Children.XLabel.String={};
fig.Children.XTickLabel={};
fig.Children.YTickLabel={};
print2eps('fig1_3_licks',gcf,'-painters')
%export_fig(gcf,sprintf('%s',['sfig1_2_v1_1','.eps']))
%% v1_1
set(gcf,'Color','w');
fig=gcf;
fig.Children.Title.String={};
print2eps('sfig1_2_v1_1',gcf,'-painters')
export_fig(gcf,sprintf('%s',['sfig1_2_v1_1','.tiff']))
%% LM_0
set(gcf,'Color','w');
fig=gcf;
fig.Children.Title.String={};
print2eps('sfig1_2_LM_0',gcf,'-painters')
export_fig(gcf,sprintf('%s',['sfig1_2_LM_0','.tiff']))
%% fig_1_4_V1imagesc
set(gcf,'Color','w');
fig=gcf;
print2eps('fig1_4_V1_imagesc',gcf,'-painters')
%% fig_1_5
set(gcf,'Color','w');
fig=gcf;
% fig.Children.Title.String={};
% fig.Children.YLabel.String={};
% fig.Children.XLabel.String={};
% fig.Children.XTickLabel={};
% fig.Children.YTickLabel={};
print2eps('fig1_5_V1diff',gcf,'-depsc2','-painters')
%% fig1_6_LMbs
set(gcf,'Color','w');
fig=gcf;
% fig.Children.Title.String={};
% fig.Children.YLabel.String={};
% fig.Children.XLabel.String={};
% fig.Children.XTickLabel={};
% fig.Children.YTickLabel={};
print2eps('fig1_6_LMbs',gcf,'-depsc2','-painters')
%% fig1_6_V1bs
set(gcf,'Color','w');
fig=gcf;
% fig.Children.Title.String={};
% fig.Children.YLabel.String={};
% fig.Children.XLabel.String={};
% fig.Children.XTickLabel={};
% fig.Children.YTickLabel={};
print2eps('fig1_6_V1bs',gcf,'-depsc2','-painters')
%% fig1_7_V1perf
set(gcf,'Color','w');
fig=gcf;
% fig.Children.Title.String={};
% fig.Children.YLabel.String={};
% fig.Children.XLabel.String={};
% fig.Children.XTickLabel={};
% fig.Children.YTickLabel={};
print2eps('fig1_7_V1perf',gcf,'-depsc2','-painters')
%% fig1_sups
%a=dir('/mnt/data/Mitra/figs/P2_L/plots/v7/figs/');
a=dir('/mnt/data/Mitra/figs/P2_L/bothdircombined/2019_07_27/');
for i=3:length(a)
   % uiopen(['/mnt/data/Mitra/figs/P2_L/plots/v7/figs/',a(i).name],1)
    uiopen(['/mnt/data/Mitra/figs/P2_L/bothdircombined/2019_07_27/',a(i).name],1)
    f=gcf;
    set(f,'Color','w')
    box('off')
    %cd('/mnt/data/Mitra/figs/P2_L/plots/v7')
    cd('/mnt/data/Mitra/figs/P2_L/bothdircombined/2019_07_27')
    %export_fig(f,sprintf('%s',[a(i).name,'.png']),'-m2','-painters')
    print2eps(sprintf('%s',[a(i).name]),gcf,'-depsc2','-painters')
end
%% fig 2 - rasters of go nogo with/out silencing

% rsters
subplot(1,2,1);xlim([-100 800])
box('off')
subplot(1,2,2);xlim([-100 800])
box('off')
%traces
xlim([-100 800]);ylim([0 18])
box('off')
axx=gca;
axx.XTick=0:100:1000;
axx.XTickLabel={'0','','200','','400','','600','','800','','1000'};
axx.YTick=0:5:20;
axx.YTickLabel={'0','5','10','15','20','25'};
set(gcf,'Color','w')
