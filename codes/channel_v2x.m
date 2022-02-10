%%channel V2X
clear all
close all

%% set tx
s = qd_simulation_parameters;                           % New simulation parameters
s.center_frequency = 2.4e9;                             % 2.4 GHz center frequency
s.use_absolute_delays = 1;                              % Include delay of the LOS path
s.sample_density = 2;                                 % Minimum possible sample density

l = qd_layout(s);
l.no_tx = 1; 
l.tx_position(:,:) = [ -200 ; 0 ; 25 ]; 
l.tx_array = qd_arrayant('multi', 8, 0.5, 12);
l.tx_array.no_elements = 10;
l.tx_array.rotate_pattern(45,'z');
l.tx_array.center_frequency = l.simpar.center_frequency;
l.tx_name = {'BS'};
%% set rx
l.no_rx = 2;
l.rx_array(1,1) = qd_arrayant('multi');
l.rx_array(1,1).center_frequency = l.simpar.center_frequency;
l.rx_array(1,2) = qd_arrayant('multi');
l.rx_array(1,2).center_frequency = l.simpar.center_frequency;
l.rx_array(1,1).no_elements = 5;
l.rx_array(1,2).no_elements = 5;

%% t1
t1 = qd_track('linear',200,3*pi/4);                    
t1.name = 'Car1';                                   
t1.initial_position(:,1) = [ 100 ; -200 ; 2 ];        
% t1.set_speed( 80/3.6 ); 
% t1.interpolate('time',10e-3,[],[],1);              
t1.movement_profile = [ 0,0 ; 3,60 ; 6,150 ; 8,get_length(t1)]';
dist1 = t1.interpolate_movement( 1e-3 );  
t1.interpolate('time',10e-3,[],[],1);
t1.calc_orientation;                                   
%% t2 
t2 = qd_track('linear',150,3*pi/4);                      
t2.name = 'Car2';                                   
t2.initial_position(:,1) = [ 120 ; -250 ; 2 ];          
% t2.set_speed( 100/3.6 ); 
% t2.interpolate('time',10e-3,[],[],1);                    % Interpolate to to 10 ms grid
t2.movement_profile = [ 0,0 ; 3,30 ; 6,120 ; 8,get_length(t2)]';
dist2 = t2.interpolate_movement( 1e-3 );  
t2.interpolate('time',10e-3,[],[],1);
t2.calc_orientation;                                    
%%
No_seg = 8;
Snaps_per_seg_1 = round(t1.no_snapshots/10);
Snaps_per_seg_2 = round(t2.no_snapshots/10);

t1.segment_index = [0:1:No_seg-1]*Snaps_per_seg_1+1;
t2.segment_index = [0:1:No_seg-1]*Snaps_per_seg_2+1;

S1 = '3GPP_38.901_UMa_NLOS';
S2 = '3GPP_38.901_UMa_LOS';
S3 = '3GPP_38.901_UMa_NLOS';
S4 = '3GPP_38.901_UMa_LOS';
S5 = '3GPP_38.901_UMa_LOS';
S6 = '3GPP_38.901_UMa_NLOS';
S7 = '3GPP_38.901_UMa_LOS';
S8 = '3GPP_38.901_UMa_LOS';


t1.scenario = {S1, S2, S3, S4, S5, S6, S7, S8};
t1.interpolate_positions( s.samples_per_meter );

t2.scenario = {S1, S2, S3, S4, S5, S6, S7, S8};
t2.interpolate_positions( s.samples_per_meter );
%% rx_rack

l.rx_track(1,1) = t1;
l.rx_track(1,2) = t2;
l.rx_track.correct_overlap;
%% channel generation
p = l.init_builder;
gen_parameters(p);
c = get_channels(p);
cn = c.merge;


%%
c1 = cn(1).interpolate( dist1 );
c2 = cn(2).interpolate( dist2 );
h1 = c1.fr( 10e6,64 );
h2 = c2.fr( 10e6,64 );
h = cat(1,h1,h2);

%% plt h
plt_h_t = squeeze(h(1,1,:,:));
mesh(abs(plt_h_t));

%% power
pow = 10*log10(reshape(sum(abs(h(: ,: ,: ,:)).^2,3),[],1));
time = (0:c1.no_snap-1) *0.001;

%% plot
% figure;
% plot (time , pow);
% xlabel('Time[s]');
% ylabel('Power[dB]');
% title('Path Gains');

% p1 = pow(1:100*r);
% p2 = pow(100*r+1 : 150*r);
% p3 = pow(150*r+1 : 240*r);
% p4 = pow(240*r+1 : 300*r);
% p5 = pow(300*r+1 : 400*r);
% p6 = pow(400*r+1 : end);
% 
% stat_p = [p1;p2;p3;p4;p5;p6];
% g  = [1*ones(length(p1),1);2*ones(length(p2),1);3*ones(length(p3),1);...
%    4*ones(length(p4),1);5*ones(length(p5),1);6*ones(length(p6),1)];
% % path_gains(:,:) = c.coeff(1,1,:,:);
% 
% figure
% boxplot(stat_p, g);
% xlabel('Segments');
% ylabel('Averaging path gains[dB]');
% title('Statistics');
% % figure
% % boxplot(10*log10(abs(path_gains))');

set(0,'defaultTextFontSize', 18)                      	% Default Font Size
set(0,'defaultAxesFontSize', 18)                     	% Default Font Size
set(0,'defaultAxesFontName','Times')               	    % Default Font Type
set(0,'defaultTextFontName','Times')                 	% Default Font Type
set(0,'defaultFigurePaperPositionMode','auto')       	% Default Plot position
set(0,'DefaultFigurePaperType','<custom>')             	% Default Paper Type
set(0,'DefaultFigurePaperSize',[14.5 7.8])            	% Default Paper Size

[ map,x_coords,y_coords] = l.power_map( '3GPP_38.901_UMa_LOS','quick',5,-500,500,-500,500,1.5 );
P = 10*log10( sum(sum(cat(3,map{:}),3),4));                    % Total received power

l.visualize;                                   % Show BS and MT positions on the map
hold on
imagesc( x_coords, y_coords, P );                       % Plot the received power
hold off
axis([-300 300 -300 300])                               % Plot size
caxis( max(P(:)) + [-20 0] )                            % Color range 
colmap = colormap;
colormap( colmap*0.5 + 0.5 );                           % Adjust colors to be "lighter"
set(gca,'layer','top')                                  % Show grid on top of the map
