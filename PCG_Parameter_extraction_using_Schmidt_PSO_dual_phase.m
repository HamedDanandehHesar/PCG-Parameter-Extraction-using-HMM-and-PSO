addpath('Springer-Segmentation');

%%
close all;
clear 

%% Load Data 
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load(file);
PCG = data.PCG;
fs = data.fs;
ECG = data.ECG;
%% Load the default options:
% These options control options such as the original sampling frequency of
% the data, the sampling frequency for the derived features and whether the
% mex code should be used for the Viterbi decoding:
springer_options = default_Springer_HSMM_options;
load('HMM_init_parameters_Springer.mat')
springer_options.audio_Fs = fs;

%% Run the HMM on an unseen test recording:
% And display the resulting segmentation
pcg = PCG';
[assigned_states] = runSpringerSegmentationAlgorithm(pcg, springer_options.audio_Fs, B_matrix, pi_vector, total_obs_distribution, true);

PCG_normalized = (pcg-min(pcg))/(max(pcg)-min(pcg))-0.5;

S1_segments = (assigned_states==1);
S2_segments = (assigned_states==3);

ECG_normalized = 1*(ECG-min(ECG))/(max(ECG)-min(ECG))-0.5;
figure(2),plot(PCG_normalized),hold on,plot(S1_segments,'r'),plot(ECG_normalized-.25,'k'),hold off
ylim([-1 1.3])


%% Get S1 peaks

sig_diff = diff([0; S1_segments]);
[index1,~,~] = find(sig_diff==1);
[index2,~,~] = find(sig_diff==-1);
S1_peaks = [];
for i=1:min(length(index1),length(index2))
    %%%% this to ensure that we find max peaks
    ind = max(index1(i)-0.05*fs,1):min(index2(i)+0.05*fs,length(pcg));
    temp =  (pcg(ind));
    [~,idx] = max(temp);
    S1_peaks = [S1_peaks ind(1)+idx];

end
% %% remove the last S1 peak
% S1_peaks(end) = [];
% figure(4),plot(PCG_normalized),hold on,plot(ECG_normalized-.25,'k'),plot(S1_peaks,abs(PCG_normalized(S1_peaks)),'*r'),plot(S1_segments,'r'),hold off
% ylim([-1.5 1.5])
% title('S1 Peaks')


%% Get S2 peaks

sig_diff = diff([0; S2_segments]);
[index1,~,~] = find(sig_diff==1);
[index2,~,~] = find(sig_diff==-1);
S2_peaks = [];
for i=1:min(length(index1),length(index2))
    %%%% this to ensure that we find max peaks
    ind = max(index1(i)-0.05*fs,1):min(index2(i)+0.05*fs,length(pcg));
    temp =  (pcg(ind));
    [~,idx] = max(temp);
    S2_peaks = [S2_peaks ind(1)+idx];

end
%% remove the last S1 peak
% S2_peaks(end) = [];
% figure(5),plot(PCG_normalized),hold on,plot(ECG_normalized-.25,'k'),plot(S2_peaks,abs(PCG_normalized(S2_peaks)),'*r'),plot(S2_segments,'r'),hold off
% ylim([-1.5 1.5])
% title('S2 Peaks')

%% Phase wrapping (New Method)




bins = 0.5*fs;
PCG = pcg';
length_sig = length(PCG);
[Phase_PC, Omega] = calculate_dual_phase_S1_S2_custom(S1_peaks, S2_peaks, length_sig, fs);%% goooooood
[Phase_DTW, Omega_DTW] = calculate_dual_phase_S1_S2_DTW( ...
    pcg, S1_peaks, S2_peaks, length_sig, fs);

figure(7),plot(6*PCG_normalized),hold on,plot(Phase_DTW,'g'),plot(ECG_normalized-.25,'k'),plot(S2_peaks,abs(6*PCG_normalized(S2_peaks)),'*r'),hold off

[~, pcg_mean_dual] = pcgsd_extractor_simple_ver_MAD_1(PCG,Phase_PC, bins);

figure(8),plot(pcg_mean_dual)

%%
pcg_mean = pcg_mean_dual;

%% new lines
pcg_mean = pcg_mean_dual;
bins = length(pcg_mean);
phase = -pi:2*pi/(bins-1):pi;

L_num_of_kernels = 35;

pcg_mean_temp = zeros(1,bins);

ai = [];
bi = [];
tetai  = [];
fi = [];
phi = [];

%% پیدا کردن پالس های قوی

[pks,locs] = findpeaks(abs(pcg_mean),'MinPeakHeight',0.01*max(abs(pcg_mean)));

[~,idx] = sort(pks,'descend');
locs = locs(idx);

theta_candidates = phase(locs);
E_total = norm(pcg_mean.^2);

for i = 1:min(L_num_of_kernels,length(theta_candidates))

disp(num2str(i))

theta_i = theta_candidates(i);

pcg_mean_temp1 = pcg_mean - pcg_mean_temp;

%% حدود پارامترها

lb = [ -1.2*max(abs(pcg_mean))   0.001   -pi   -pi ];
ub = [  1.2*max(abs(pcg_mean))   1       pi    pi ];

%% تابع هزینه

myfun = @(x) norm( pcg_mean_temp1 - ...
x(1).*exp(-(rem(phase-theta_i+pi,2*pi)-pi).^2./(2*x(2).^2)) ...
.*cos(x(3).*phase-x(4)) ,2);

options = optimoptions('particleswarm','SwarmSize',200,'MaxIter',80);

OptimumParams = particleswarm(myfun,4,lb,ub,options);

%% پارامترها

ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
fi_1 = OptimumParams(3);
phi_1 = OptimumParams(4);

tetai_1 = theta_i;

ai = [ai ai_1];
bi = [bi bi_1];
tetai = [tetai tetai_1];
fi = [fi fi_1];
phi = [phi phi_1];

%% ساخت کرنل

dt = rem(phase - tetai_1 + pi,2*pi)-pi;

kernel = ai_1 .* exp(-dt.^2./(2*bi_1.^2)) .* cos(fi_1.*phase-phi_1);

pcg_mean_temp = pcg_mean_temp + kernel;

figure(12)
plot(pcg_mean_temp)
drawnow

E_model = norm(pcg_mean_temp.^2);
    %% شرط پایان (۹۹٪ انرژی)
    % if (E_model / E_total) >= 0.999
    %     disp('*** 99% of signal energy reconstructed → stopping ***');
    %     break
    % end
end


figure(11),plot(pcg_mean),hold on,plot(pcg_mean_temp),hold off
PCG_synth = zeros(1,size(PCG,2));
Phase = Phase_DTW;
for kk = 1:length(ai)
    dtetai = rem(Phase - tetai(kk) + pi,2*pi)-pi;
    PCG_synth = PCG_synth + ai(kk) .* exp(-dtetai .^2 ./ (2*bi(kk) .^ 2)).*cos(fi(kk).*Phase-phi(kk));
end
figure(13),plot(PCG),hold on,plot(PCG_synth),hold off






function [pcgsd,pcg_mean] = pcgsd_extractor_simple_ver_MAD_1(pcg,phase,bins)

x1 = pcg(1,:);

meanPhase = zeros(1,bins);
ECGmean   = zeros(1,bins);
ECGsd     = zeros(1,bins);

% -------- تابع محاسبه ميانگين وزن‌دار robust --------
function wm = robust_gauss_mean(vals)

    m = median(vals);

    MAD = median(abs(vals - m));
    h = 1.4826 * MAD;

    if h < eps
        wm = m;
        return
    end

    w = exp(-((vals - m).^2) / (h^2));
    wm = sum(w .* vals) / sum(w);

end

% --------- bin اول ----------
I = find( phase >= (pi - pi/bins) | phase < (-pi + pi/bins) );

if ~isempty(I)

    vals = x1(I);

    meanPhase(1) = -pi;
    ECGmean(1) = robust_gauss_mean(vals);
    ECGsd(1) = std(vals);

else

    meanPhase(1) = 0;
    ECGmean(1) = 0;
    ECGsd(1) = -1;

end

% -------- بقيه bin ها --------
for i = 1:bins-1

    I = find( phase >= 2*pi*(i-0.5)/bins - pi & ...
              phase <=  2*pi*(i+0.5)/bins - pi );

    if ~isempty(I)

        vals = x1(I);

        meanPhase(i+1) = mean(phase(I));
        ECGmean(i+1) = robust_gauss_mean(vals);
        ECGsd(i+1) = std(vals);

    else

        meanPhase(i+1) = 0;
        ECGmean(i+1) = 0;
        ECGsd(i+1) = -1;

    end

end

% -------- پر کردن bin های خالی --------
K = find(ECGsd == -1);

for i = 1:length(K)

    k = K(i);

    switch k

        case 1
            meanPhase(k) = -pi;
            ECGmean(k) = ECGmean(k+1);
            ECGsd(k) = ECGsd(k+1);

        case bins
            meanPhase(k) = pi;
            ECGmean(k) = ECGmean(k-1);
            ECGsd(k) = ECGsd(k-1);

        otherwise
            meanPhase(k) = mean([meanPhase(k-1), meanPhase(k+1)]);
            ECGmean(k) = mean([ECGmean(k-1), ECGmean(k+1)]);
            ECGsd(k) = mean([ECGsd(k-1), ECGsd(k+1)]);

    end

end

pcg_mean = ECGmean;
pcgsd = ECGsd;

end





function [Phase, Omega] = calculate_dual_phase_S1_S2_custom(S1_peaks, S2_peaks, length_sig, fs)

% -------------------------------------------------------------------------
% Custom dual-phase mapping for PCG:
%   S1 = 0
%   S2 = pi/2
%   Next S1 = 2*pi (which will wrap to 0)
%
% This produces a **stable cardiac phase** for mean PCG extraction.
% -------------------------------------------------------------------------

Phase = zeros(1, length_sig);
Omega = zeros(1, length_sig);

% ----------------------------
% 1) Combine and sort peaks
% ----------------------------
all_peaks  = [S1_peaks(:); S2_peaks(:)];
labels     = [ones(numel(S1_peaks),1); 2*ones(numel(S2_peaks),1)];
[pk, idx]  = sort(all_peaks);
lb         = labels(idx);

% ----------------------------
% 2) Ensure starting with S1
% ----------------------------
if lb(1) == 2   % starts with S2 → remove it
    pk(1) = []; 
    lb(1) = [];
end

% ----------------------------
% 3) Enforce exact alternation: S1 S2 S1 S2 ...
% ----------------------------
i = 1;
while i < numel(lb)
    if lb(i) == lb(i+1)
        pk(i+1) = [];
        lb(i+1) = [];
    else
        i = i + 1;
    end
end

% ----------------------------
% 4) Assign custom target phases
%     S1(k) = 2*pi*n
%     S2(k) = pi/2 + 2*pi*n
%     where n = floor((k-1)/2)
% ----------------------------
target_phase = zeros(size(pk));

for k = 1:numel(pk)
    n = floor((k-1)/2);  % cycle index: S1,S2 share same n
    if lb(k) == 1
        target_phase(k) = 2*pi*n;
    else
        target_phase(k) = (pi/2) + 2*pi*n;
    end
end

% ----------------------------
% 5) Linear interpolation between consecutive peaks
% ----------------------------
for k = 1:numel(pk)-1
    
    i1 = pk(k);
    i2 = pk(k+1);
    
    th1 = target_phase(k);
    th2 = target_phase(k+1);

    bins = i2 - i1;

    dtheta = (th2 - th1) / bins;
    omega_val = fs * dtheta;

    theta = th1;
    Phase(i1) = wrapToPi(theta);
    Omega(i1) = omega_val;

    for j = i1+1:i2
        theta = theta + dtheta;
        Phase(j) = wrapToPi(theta);
        Omega(j) = omega_val;
    end
end

% ----------------------------
% 6) After last peak — continue with average cycle
% ----------------------------
if numel(pk) > 1
    RR = mean(diff(pk));
else
    RR = fs; % fallback
end

dtheta = 2*pi / RR;      % same as S1→next S1
omega_val = fs * dtheta;
theta = target_phase(end);

for j = pk(end)+1:length_sig
    theta = theta + dtheta;
    Phase(j) = wrapToPi(theta);
    Omega(j) = omega_val;
end
%% 8) Before first peak  (NEW FIX)

theta = target_phase(1);

for j = pk(1)-1:-1:1
    theta = theta - dtheta;
    Phase(j) = wrapToPi(theta);
    Omega(j) = omega_val;
end

end


function [Phase, Omega] = calculate_dual_phase_S1_S2_DTW( ...
    pcg, S1_peaks, S2_peaks, length_sig, fs)

% -------------------------------------------------------------------------
% Dual-phase mapping for PCG using DTW between peaks
%
%   S1 = 0
%   S2 = pi/2
%   Next S1 = 2*pi (wraps to 0)
%
% Same structure as the linear version, but interpolation between
% peaks is replaced with DTW-based nonlinear warping.
% -------------------------------------------------------------------------

Phase = zeros(1,length_sig);
Omega = zeros(1,length_sig);

env = abs(hilbert(pcg));

%% 1) Combine and sort peaks
all_peaks  = [S1_peaks(:); S2_peaks(:)];
labels     = [ones(numel(S1_peaks),1); 2*ones(numel(S2_peaks),1)];
[pk, idx]  = sort(all_peaks);
lb         = labels(idx);

%% 2) Ensure starting with S1
if lb(1)==2
    pk(1)=[];
    lb(1)=[];
end

%% 3) Enforce alternation
i=1;
while i<numel(lb)
    if lb(i)==lb(i+1)
        pk(i+1)=[];
        lb(i+1)=[];
    else
        i=i+1;
    end
end

%% 4) Target phases
target_phase=zeros(size(pk));

for k=1:numel(pk)

    n=floor((k-1)/2);

    if lb(k)==1
        target_phase(k)=2*pi*n;
    else
        target_phase(k)=pi/2+2*pi*n;
    end

end

%% 5) DTW interpolation between peaks

for k=1:numel(pk)-1

    i1=pk(k);
    i2=pk(k+1);

    th1=target_phase(k);
    th2=target_phase(k+1);

    seg = env(i1:i2);
    N = length(seg);

    ref = linspace(0,1,N);          % reference time
    phase_ref = linspace(th1,th2,N);

    sig = seg(:)';
    refsig = resample(sig,N,length(sig));

    [~,ix,iy] = dtw(refsig,sig);

    phase_seg = zeros(1,N);

    for n=1:length(ix)
        phase_seg(iy(n)) = phase_ref(ix(n));
    end

    phase_seg = fillmissing(phase_seg,'linear');

    Phase(i1:i2)=wrapToPi(phase_seg);

    dtheta = mean(diff(phase_seg));
    Omega(i1:i2)=fs*dtheta;

end

%% 6) RR estimation

if numel(pk)>1
    RR=mean(diff(pk));
else
    RR=fs;
end

dtheta=2*pi/RR;
omega_val=fs*dtheta;

%% 7) After last peak

theta=target_phase(end);

for j=pk(end)+1:length_sig

    theta=theta+dtheta;
    Phase(j)=wrapToPi(theta);
    Omega(j)=omega_val;

end

%% 8) Before first peak

theta=target_phase(1);

for j=pk(1)-1:-1:1

    theta=theta-dtheta;
    Phase(j)=wrapToPi(theta);
    Omega(j)=omega_val;

end

end



