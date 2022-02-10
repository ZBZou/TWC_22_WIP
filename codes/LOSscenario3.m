%% Date: Jan 24, 2016 6:00 pm. This is the third code from Quadriga (with some modification)


%% Initial Parameters

close all
clear all

%% General simlation parameters
Nt = 10;  %numbe of trials
%  v1 = 15;
% v2 = 25;
LL = zeros(10,Nt);
hDB = cell(10,1);                       % DB for Channel Gains

for ii = 1 : Nt
s = qd_simulation_parameters;          % Basic simulation parameters
s.center_frequency      = 5.9e9;    % Center Frequency is 5.9 GHz
s.sample_density        = 2.5;        % 2 samples / half wave length
s.use_absolute_delays = 1;          % Delays are "on"

%%

% Velocity = round((v2-v1)*rand(1) )+ v1;
    Velocity = 20;       % Vehicle velocity measured in m/s.
%     lambda = (3e8/s.center_frequency);
%     d_coh = lambda/(2*pi);
%     UpdateTime = d_coh/Velocity;       % The channel coefficents are generated each T_coherent
    UpdateTime = 0.25;
    UpdateDistance = Velocity * UpdateTime; % Distance (in meter) between the
                                            % adjacent snapshots (i.e. the channel
                                            % coefficents are calculated each UpdateDistance).
                                            % UpdateDistance = d_coh
    
    SegmentLength = 10;  % Segment length (in meter), which must be less than
    %or equal the decorelation distance of the LSPs
    
    SnapshotsPerSegment = SegmentLength/UpdateDistance; % Number of snapshots per segment
    if SnapshotsPerSegment - floor(SnapshotsPerSegment) < 0.5; % This is to make sure the the number of snapshots is an integer.
        SnapshotsPerSegment = floor(SnapshotsPerSegment);
    else SnapshotsPerSegment = ceil(SnapshotsPerSegment);
    end
    
    TrackLength = 200;   % Track length (in m) (i.e. length of the road stretch) 
    NmberOfSegments = TrackLength/SegmentLength;
    if rem(NmberOfSegments,1)~=0  % This is to make sure the the number of segmants is an integer.
        error('Number of segments must be an integer')
    end
    
    NumberOfSnapshots = SnapshotsPerSegment * NmberOfSegments; %No. of snapshots in the track
    
    Y = 1;
    for i = 1 : NmberOfSegments
        X(i) = Y;                % X is avector containing the initial snapshot index for each segment
        Y = Y + SnapshotsPerSegment;
    end
    
    
    %% Define a track
    
    t = qd_track('linear',TrackLength,pi);     %Track length is 200m west (<---)
    t.initial_position = [100;1.5;0];       % Starting point of the trck
    t.interpolate_positions( s.samples_per_meter );
    t.segment_index = X;
    S1 = 'WINNER_UMi_B1_LOS1';
    S2 = 'WINNER_UMi_B1_LOS2';
    S3 = 'WINNER_UMi_B1_LOS3';
    S4 = 'WINNER_UMi_B1_LOS4';
    S5 = 'WINNER_UMi_B1_LOS5';
    S6 = 'WINNER_UMi_B1_LOS6';
    S7 = 'WINNER_UMi_B1_LOS7';
    S8 = 'WINNER_UMi_B1_LOS8';
    S9 = 'WINNER_UMi_B1_LOS9';
    S10 = 'WINNER_UMi_B1_LOS10';
    
    
    t.scenario = {S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S10, S9, S8, S7, S6, S5, S4, S3, S2, S1};
    t.interpolate_positions( s.samples_per_meter );
    
    
    %% Define network layout
    
    l=qd_layout(s);
    l.tx_position = [0;4.5;2.5];                 % 2.5 m tx height
    l.rx_position = [100;1.5;1.5];               % 1.5 m tx height
    l.tx_array.generate( 'omni' );
    l.rx_array = l.tx_array;
    
    % l.visualize;
    % % view(-33, 60);
    % axis([-110 110 -7 7 0 3 ]);
    
    l.track = t;
    
    % l.set_pairing('sf',10,26.989)
    % [map, pox, poy] = power_map('WINNER_UMi_B1_LOS', 'quick',10, -100,0,100,0, 27);
    
    %% LSPs generations
    
    % n = [1 1 1 2 3 4  2 3 4 2];          % Number of fixed scatterers (only on the
    % side of the road) in each segment
    n = [ 2 3 4  4 4 5 5];
    
    
%     Nt = 2;
    
    % jj = 1;
    % kk = 1;
%     LL = zeros(10,Nt);
%     for ii = 1 : Nt
        L=[];
        n2=[1 1 1] + round(3*rand(1));
        K = 15;                % Max number of additional scatteres (on the street)
        for i = 1 : length(n)
            L(i) = n(i) + round(rand(1)*K); %Total no. of scatteres
        end
        L=[n2 L];
        LL(:,ii) = L;     %matrix of L for each run
        
        k = 10;
        for i = 1:10
            L2(i) = L(k);
            k = k-1;
        end
        L_total = [L L2];
        
        % h_hor = l.rx_position(1); %initial TX-RX horizontal distance
        % h_RSU = l.tx_position(3);    % RSU height
        % h_Vehicle = l.rx_position(3); % Vehicle height
        % h = h_RSU - h_Vehicle;      % Absolute hight
        % r = sqrt( h^2 + h_hor^2);   %Max Comm range
        % MaxSegPerEllipse = sqrt(r)/2; %Max number of segments per Ellipse
        % Prob = 0.5;
        % MinSegPerEllipse = round(Prob * MaxSegPerEllipse); %minimum number of segments with no scattering
        % d = MinSegPerEllipse * SegmentLength;  %Length of the minimum number of accounted segments
        
        p = l.init_builder;
        gen_parameters(p);
        
        % MinSegPerEllipse = 3;
        % n2 = [1 1 2];
        % for i = 1 : MinSegPerEllipse
        %     p(i).scenpar.NumClusters = n2(i); %Randomly change the total no. of clusters
        %     p(i).scenpar.xpr_mu      = 100;              % Disable XPR
        %     p(i).scenpar.xpr_sigma   = 0;
        %     p(i).update_parameters;
        % end
        for i = 1:10
            p(i).scenpar.NumClusters = L(i);
%             p(i).scenpar.xpr_mu      = 100;              % Disable XPR
%             p(i).scenpar.xpr_sigma   = 0;
%             p(i).update_parameters;
        end
        % L_Total = [1 1 1 L 1 1 1];
        
        
        % p(4).scenpar.NumClusters = 1;
        % p.update_parameters;
        % p(5).scenpar.NumClusters = 1;
        % p.update_parameters;
        % p(6).scenpar.NumClusters = 1;
        % p.update_parameters;
        
        % p.scenpar.xpr_mu      = 100;              % Disable XPR
        % p.scenpar.xpr_sigma   = 0;
        
        c = get_channels(p);
        cn = c.merge;
        
        t.set_speed( Velocity );       %Car speed is 30km/h =8.3m/s
        dist = t.interpolate_movement( UpdateTime); % Specifies the number of snapshops used in pdp
        %If we use value of 1s here, then
        % d=v*t = 8.3m/s*1s=8.3m --> for 200m
        %track length, we have d=200/8.3=25 snapshots
        %(i,e. positions). Also, we need 25
        %seconds to pass the 200m track. Thus,
        %if segment length (typically for
        %LOS)=20m, then 25/20=1 snapshop per
        %segment.
        ci = cn.interpolate( dist , 'linear' );
        
        %% Plotting
        %     pwr_orig = 10*log10(squeeze(abs(cn.coeff(1,1,1,:))).^2);
        %     nsnap = cn.no_snap;
        %     dist_orig = (0:nsnap-1) * t.get_length/(nsnap-1);
        %     pwr_int = 10*log10(squeeze(abs(ci.coeff(1,1,1,:))).^2);
        %
        %     figure
        %     plot( dist_orig,pwr_orig , 'r','Linewidth',2 )
        %     hold on
        %     plot( dist,pwr_int ,'b','Linewidth',2 )
        %     hold off
        %     axis(  [ min(dist) , max(dist) ,...
        %         min( pwr_orig( pwr_orig>-Inf ) ) , ...
        %         max( pwr_orig( pwr_orig>-Inf ) )+10 ] )
        %
        %     xlabel('Distance from start point [m]');
        %     ylabel('Power [dB]');
        %     title('RX Power before and after interpolation')
        %     grid on
        
        
        %% RX Power distribution along the track when the RSU is placed at the center of the track
        %     dist2=dist-100; % Making the center of the distance vector zero
        %     pwr_int2 =(squeeze(abs(ci.coeff(1,1,1,:))).^2);  % RX power in Watt.
        %     pwrNormalized = pwr_int2/max(pwr_int2); % Power normalizeation
        %     figure
        %     plot( dist2,pwrNormalized ,'b','Linewidth',2 )
        %     xlabel('Distance from the RSU [m]');
        %     ylabel('Power [Watt]');
        %     title('Normalized RX power distribution along the road (RSU@center)')
        %     grid on
        %
        %     figure
        %     hist(pwrNormalized,sqrt( length(pwrNormalized)))
        %     hold on
        %     hist(-pwrNormalized,sqrt( length(pwrNormalized)),'r')
        %     xlabel('Normalized RX power')
        %     ylabel('Number of samples')
        %     title('RX power distribution of the propagation paths')
        %     grid on
        
        %% Histogram Plot of the Laplace pdf
        % N = length(pwrNormalized);
        % mu = mean(pwrNormalized);
        % sigma = std(pwrNormalized);
        % pwrNormalized = mu+pwrNormalized;
        % fx = 1/sqrt(2*(sigma^2))*exp(-sqrt(2/(sigma^2))*abs(pwrNormalized-mu)); % calculate the Laplace PDF
        %
        % bins = sqrt(length(pwrNormalized));
        % range = max(pwrNormalized)-min(pwrNormalized);
        % width = range/bins;
        % area = N*width;
        % fx = area*fx;
        % figure
        % hist(pwrNormalized,bins)
        % hold on
        % plot(pwrNormalized,fx,'-r')
        %% PDP Plot
        h = ci.fr( 10e6,64 );
        h = squeeze(h);
        pdp = 10*log10(abs(ifft(h,[],1).').^2);
        
        H_FirstLocation(:,ii) = h(:,1);  %for the first location only
        H_MiddleLocation(:,ii) = h(:,floor((size(h,2))/2) );  %for the middle location only
        H_LastLocation(:,ii) = h(:,end);  %for the last location only
%         H_ArbitraryLocation(:,ii) = h(:,10);  %for the 55th location only
        %    H(jj:ii*size(h,1), kk:ii*size(h,2)) = h; %for all locations:this brought out-of-memory Matlab message
        %    jj = jj + size(h,1);
        %    kk = kk + size(h,2);
        
        for iii = 1:10
            hDB{iii}=[hDB{iii},h(:,2*iii)];
        end
        counter =  ii
    end
    
    H_v1_FirstLocation = H_FirstLocation;
    H_v1_MiddleLocation = H_MiddleLocation;
    H_v1_LastLocation = H_LastLocation;
%     H_v1_ArbitraryLocation = H_ArbitraryLocation;
    ZZ=(abs(H_MiddleLocation.')).^2;
    Z =surf(10*log10(ZZ));
    
  
figure
imagesc(pdp(:,:));
caxis([ max(max(pdp))-50 max(max(pdp))-5 ])
colorbar

cm = colormap('hot');
colormap(cm(end:-1:1,:));

% set(gca,'XTick',1:4:64);
% set(gca,'XTickLabel',(0:8:64)/10e6*1e6);
xlabel('Delay [\mus]');
% set(gca,'YTick',1:ci.no_snap/8:ci.no_snap);
% set(gca,'YTickLabel', (0:ci.no_snap/8:ci.no_snap)/ci.no_snap * 20 );
ylabel('Distance [m]');
title('Pwer-delay profile')

%% Power in the first delaty (corresponding to the first OFDM sub-carrier)
figure
plot(pdp(:,1));
xlabel('Distance (i.e. snapshot) [m]');
ylabel('Power [dB]');
title('Power of the first delay tap')
grid on
%% 3D plot of pdp a
figure
pdp_3D = surf(pdp);
xlabel('Delay [\mus]');
ylabel('Distance [m]');
zlabel('Channel gain [dB]')
title('3D plot of pdp')

%% Channel selectivity plot
pfp = 10*log10(abs((h).').^2);       % power-frequency profile

figure
pfp_3D = surf(pfp);
xlabel('Frequency (OFDM subcarriers) [MHz]');
ylabel('Distance(channel snapshots) [m]');
zlabel('Channel gain [dB]')
title('3D plot of pfp')



%% Plotting the distribution of the scatterers as a discrete RV
figure
stem(1:20,L_total)
% hold on
% plot(1:20,L_total,'r')
xlabel('Segment number');
ylabel('Number of scatterers')
title('Distribution of the scatteres along the segments of the track')
grid on
% figure
% stem(1:20,L_total/max(L_total))
% hold on
% plot(1:20,L_total/max(L_total),'r')
% xlabel('Number of scatterers');
% ylabel('pdf')
% title('pdf of the number of scatteres in each segment')


%% Other plots
% figure
% plot(pdp)
%
% figure
% plot(mean(pdp))
%
% figure
% plot(pfp)
%
% figure
% plot(mean(pfp))












