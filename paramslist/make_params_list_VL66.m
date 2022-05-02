% assigning channels:
n_ephys_channles = 128;
n_aux_channels = 6;
photodiode_channel = 2;
laser_channel = 1;

% the freq for filtering pd signal:
pdfiltfreq = 100;

% threshold for detectinf onset of pd and laser (on derivative):
DiffTh=1;
DiffThL = 2000;
lasershape = 'ramp';

% recording file params:
% where you expect SST messages to arrive, you can check the table beforehand
% to make sure the node id is correct.

Netevent_nodeID = 103;
rec_number = '0';
processor_number = '101';
animalname='VL66';
%expname={'/stimlaserL_2018-01-16_18-59-01'};
expname={'/stimlaserS_2018-01-16_20-16-27'};

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

do_parse_ssmessages = 1;
% in case of skipping ParseSSmessages, tStart and tend must be given as a
% heuristic for start and end of pd signal, also in cutting experiments,
% tEnd must be set manually to cutting threshold


% L
tStart_p = 0.6*10^6; 
tEnd_p =   129749424;

tEnd = tEnd_p;
tStart = tStart_p;

cut_rec = 0;
tEnd_cut =tEnd_p;


do_behavior_analysis = 1;
lick_channel = 4;
running_channel = 3;

grooming_speed_threshold = 4;
grooming_lick_threshold = 0.1;


% saving directory
savingdirroot = '/mnt/data/Mitra/figs/P1_S/LED/FB/VL66/'; 
