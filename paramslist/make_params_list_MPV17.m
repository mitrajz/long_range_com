

% assigning channels:
n_ephys_channles = 128;
n_aux_channels = 6;
photodiode_channel = 2;
laser_channel = 1;
lick_channel = 6; % analog lick channel
valve_channel = 5;
running_channel = 4;
% beh: lick threshold fo grooming: (stationary mice, no speed threshold)
grooming_lick_threshold = 0.5;
% the freq for filtering pd signal:
pdfiltfreq = 100;

% threshold for detecting onset of pd and laser (on derivative):
DiffTh=5;
DiffThL = 2000;
lasershape = 'ramp';

% recording file params:
% where you expect SST messages to arrive, you can check the table beforehand
% to make sure the node id is correct.

Netevent_nodeID = 103;
rec_number = '0';
processor_number = '101';
animalname='MPV17';
%expname={'/stimlaserL_2018-01-17_19-08-45'};
rotate = 1;
if ~rotate
expname={'/task_2019-01-30_21-11-58'};
oemessagesfilename = '30-01-19-21-12-46.dat';
else
expname={'/taskn_2019-01-30_22-23-19'};
oemessagesfilename = '30-01-19-22-23-49.dat';
end
% spike extraction
units_to_extract = 'gmu'; % could be gmu: good and mu, or 'all': all except noise

% window for PAllOn
WindowSec = 1;
samplingf = 30000;

% laserdelay binning params
% L
Laserdelaybinnum=2000;
Laserdelaybinwindow=3;

% set stimulus and behavior analysis preferances

% L
if ~rotate
    tStart = 1*10^6;
    tEnd =  122588192;
else
    tStart = 1*10^6;
    tEnd =  118402288;
end
%%%%%% params for detecting licks:
if rotate
lickthprctile = 95;
lickcontforms = 20;
thr = 0.2; 

else
lickthprctile = 99.5;
lickcontforms = 40;
thr = 0.1;
end

% saving directory
savingdirroot = '/mnt/data/Mitra/figs/MPV17/'; 
