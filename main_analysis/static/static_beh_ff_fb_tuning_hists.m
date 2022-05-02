%%% first clear all, then run the first few sections of
%%% staticBehaviorpropertiesVSsilencing.m, once for FF and once for FB.
%%% Then the FF and FB cells will be used here. The post plots will depebd
%%% on include_phys (set in  staticBehaviorpropertiesVSsilencing.m). p
%%% values are ranksum

%%% this is pre cleanign for nans
figure;
subplot(1,4,1);histogram(FF.pre.all_lick_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.pre.all_lick_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','b'); ylabel('lick tuning')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.pre.all_lick_tun,FB.pre.all_lick_tun)))

subplot(1,4,2);histogram(FF.pre.all_earlystim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.pre.all_earlystim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','b'); ylabel('early stim tuning')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.pre.all_earlystim_tun,FB.pre.all_earlystim_tun)))

subplot(1,4,3);histogram(FF.pre.all_earlystim_tun_sign,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.pre.all_earlystim_tun_sign,[-1:0.1:1],'Normalization','pdf','FaceColor','b'); ylabel('early stim tuning sign')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.pre.all_earlystim_tun_sign,FB.pre.all_earlystim_tun_sign)))

subplot(1,4,4);histogram(FF.pre.all_stim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.pre.all_stim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','b'); ylabel('stim tuning')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.pre.all_stim_tun,FB.pre.all_stim_tun)))

%%% this is post cleanign for nans (what is actually used for models)
figure;
subplot(1,4,1);histogram(FF.post.all_lick_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.post.all_lick_tun,[-1:0.1 :1],'Normalization','pdf','FaceColor','b'); ylabel('lick tuning')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.post.all_lick_tun,FB.post.all_lick_tun)))
hold on; line([nanmedian(FF.post.all_lick_tun),nanmedian(FF.post.all_lick_tun)],[0,3],'Color','r')
hold on; text(-1,1.8,sprintf('ff median = %d',nanmedian(FF.post.all_lick_tun)))
hold on; line([nanmedian(FB.post.all_lick_tun),nanmedian(FB.post.all_lick_tun)],[0,3],'Color','b')
hold on; text(-1,1.6,sprintf('fb median = %d',nanmedian(FB.post.all_lick_tun)))

subplot(1,4,2);histogram(FF.post.all_earlystim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.post.all_earlystim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','b'); ylabel('early stim tuning')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.post.all_earlystim_tun,FB.post.all_earlystim_tun)))
hold on; line([nanmedian(FF.post.all_earlystim_tun),nanmedian(FF.post.all_earlystim_tun)],[0,3],'Color','r')
hold on; text(-1,1.8,sprintf('ff median = %d',nanmedian(FF.post.all_earlystim_tun)))
hold on; line([nanmedian(FB.post.all_earlystim_tun),nanmedian(FB.post.all_earlystim_tun)],[0,3],'Color','b')
hold on; text(-1,1.6,sprintf('fb median = %d',nanmedian(FB.post.all_earlystim_tun)))

subplot(1,4,3);histogram(FF.post.all_earlystim_tun_sign,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.post.all_earlystim_tun_sign,[-1:0.1:1],'Normalization','pdf','FaceColor','b'); ylabel('early stim tuning sign')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.post.all_earlystim_tun_sign,FB.post.all_earlystim_tun_sign)))

subplot(1,4,4);histogram(FF.post.all_stim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','r')
hold on;histogram(FB.post.all_stim_tun,[-1:0.1:1],'Normalization','pdf','FaceColor','b'); ylabel('stim tuning')
legend('ff','fb')
text(-1,2,num2str(ranksum(FF.post.all_stim_tun,FB.post.all_stim_tun)))
