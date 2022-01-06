% This script is an example of pulse adjustment (8° pulse) on the 1st subject of the control database
% Matlab version : R2020b

clear
close all
inputPath = pwd;

%% Add utils

addpath(fullfile(inputPath, 'utils'))

%% Load data

% Load reference B1+ map (in Hz/V) and reference brain mask
load(fullfile(inputPath, '..', 'Data', 'B1_ref.mat'));

% Load B1+ maps of the 1st subject of the control database (in Hz/V)
load(fullfile(inputPath, '..', 'Data', 'ControlDatabase.mat'));

% Load UP
load(fullfile(inputPath, '..', 'Data', 'Pulses', 'UP.mat'));

figure; 
UP.visualizeMP();
drawnow;

[UP_E_per_ch, UP_E_tot] = UP.getEnergyPerChannel(); % V^2.s
UP_E_per_ch = max(UP_E_per_ch);
UP_Pmax = max(UP.getMaxRFPowerPerChannel());% V^2

% Load SUP (8° GRAPE pulse)
load(fullfile(inputPath, '..', 'Data', 'Pulses', 'SUP.mat'));
figure; 
UP.visualizeMP();
drawnow;

[SUP_E_per_ch, SUP_E_tot] = SUP.getEnergyPerChannel(); % V^2.s
SUP_E_per_ch = max(SUP_E_per_ch);
SUP_Pmax = max(SUP.getMaxRFPowerPerChannel());% V^2

fprintf('\n                          UP         SUP\n')
fprintf('E per channel (mJ)   %10.2f %10.2f\n', 1e3*UP_E_per_ch/50, 1e3*SUP_E_per_ch/50);
fprintf('E total       (mJ)   %10.2f %10.2f\n', 1e3*UP_E_tot/50, 1e3*SUP_E_tot/50);
fprintf('Peak power    (W)    %10.2f %10.2f\n\n', UP_Pmax/50, SUP_Pmax/50);


%% Pulse adjustment

adjtype = 'full'; % 'diagonal'
do_resample = true; % false

% Extract 3 slices to simulate a 3-slice B1+ map scan
% taking only one slice is not enough, taking more than 3 slices does not
% bring much (and costs scan time in practice)
Mask_ref_3slc = Mask_ref_40slc; 
Mask_ref_3slc(:,:,~get_slice_subset(3, size(Mask_ref_40slc, 3), 6:32)) = 0;

clear adjSUP
figure;
L = zeros(Nc,Nc,numel(datav));
for i = 1:numel(datav)
    fprintf('progress SUP_adjustment = %d/%d\n', i, numel(datav));
    data = datav(i);
    
    switch adjtype
        case 'full'
            [adjSUP(i),L(:,:,i)] = SUP_adjustment(SUP, data.b1, data.MaskSignal, b1_ref, Mask_ref_3slc);
        case 'diagonal'
            [adjSUP(i),L(:,:,i)] = SUP_adjustment_diag(SUP, data.b1, data.MaskSignal, b1_ref, Mask_ref_3slc);
        otherwise
            error('unknown adjustment type %s; expecting full|diagonal', adjtype);
    end
    
    [adjSUP_E_per_ch, adjSUP_E_tot] = adjSUP(i).getEnergyPerChannel(); % V^2.s
    adjSUP_E_per_ch = max(adjSUP_E_per_ch);
    adjSUP_Pmax = max(adjSUP(i).getMaxRFPowerPerChannel());% V^2

    lambda = max([1, adjSUP_E_per_ch/SUP_E_per_ch, adjSUP_E_tot/SUP_E_tot, adjSUP_Pmax/SUP_Pmax]); 
    
    % relax Energy and peak power constraints
    adjSUP(i) = adjSUP(i).scaletime(lambda);
    
    % resample to 10 us 
    if (do_resample)
        adjSUP(i) = adjSUP(i).sampled(1e-5);
    end
    
    subplot(1,3,1);
    imagesc(real(L(:,:,i)), [-2,2]); axis equal;
    colorbar;
    title('Re(L)')
    subplot(1,3,2);
    imagesc(imag(L(:,:,i)), [-0.5,0.5]); axis equal;
    colorbar;
    title('Im(L)')
    subplot(1,3,3);
    imagesc(abs(L(:,:,i)), [0,2]); axis equal;
    colorbar;
    title('abs(L)')
    drawnow;
end

%% Comparison to UP

FANRMSE = zeros(numel(datav), 2);
FA_UP = zeros([dim, numel(datav)]);
FA_adjSUP = zeros([dim, numel(datav)]);

figure;

for i = 1:numel(datav)
        fprintf('progress FA simulation = %d/%d\n', i, numel(datav));

    data = datav(i);
    b1 = 2*pi*reshape(data.b1, Npos, Nc).'; % rad/s/V
    
    fa = UP.stasim(pos, b1, 0); 
    fa = 180/pi*reshape(fa, dim);
    FANRMSE(i,1) = calcNRMSE(fa(data.MaskBrain>0), UP.FA, 1);
    FA_UP(:,:,:,i) = applyMask(data.MaskBrain, fa);

    fa = adjSUP(i).stasim(pos, b1, 0);
    fa = 180/pi*reshape(fa, dim);
    FANRMSE(i,2) = calcNRMSE(fa(data.MaskBrain>0), SUP.FA, 1);
    FA_adjSUP(:,:,:,i) = applyMask(data.MaskBrain, fa);

    
    subplot(1,2,1);
    view3dp(1, FA_UP(:,:,:,i), 'cscale', [0,10]);
    colormap(jet);
    title(['avg FA-NRMSE UP = ', num2str((FANRMSE(i,1))*100,2), '%'])
    colorbar
    subplot(1,2,2);
    view3dp(1, FA_adjSUP(:,:,:,i), 'cscale', [0,10]);
    colormap(jet);
    title(['avg FA-NRMSE adjSUP = ', num2str((FANRMSE(i,2))*100,2), '%']);
    colorbar
    drawnow;
end

%% FA-NRMSE comparison (scatter plot)

figure;
title ('FA-NRMSE of UP versus SUP')
plot(FANRMSE(:,1)*100, FANRMSE(:,2)*100, 'o');
hold on;
plot([0,15],[0,15], 'r');
xlabel('FANRMSE UP (%)')
ylabel('FANRMSE SUP (%)')
xlim([0,15]);
ylim([0,15]);
grid on
legend({'Subjects 1-15', 'identity line'})

%% FA maps (sagittal)

figure;
subplot(1,2,1);
view3dp(1, FA_UP, 'cscale', [0,10]);
colormap(jet);
title(['avg FA-NRMSE UP = ', num2str(mean(FANRMSE(:,1))*100,2), '%'])
colorbar
subplot(1,2,2);
view3dp(1, FA_adjSUP, 'cscale', [0,10]);
colormap(jet);
title(['avg FA-NRMSE adjSUP = ', num2str(mean(FANRMSE(:,2))*100,2), '%']);
colorbar





