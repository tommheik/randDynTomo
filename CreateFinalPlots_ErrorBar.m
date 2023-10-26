% Plotting codes by Tatiana and Luca
% Slight changes by T H   2023

notLoaded = true;

if notLoaded
    close all
    clear all
    
    lab_data = 'Stempo'; % 'Cartoon', 'Gel'
    lab_noise = 'decreasing'; % 'fixed', 'decreasing'
    lab_power = 'p1'; % 'p1', 'p32'
    lab_Nsamp = 5;
    lab_reso = 280; % 256; % 128;
    lab_transform = 'Wavelet'; % 'Shearlet'/'Wavelet';
    dateString = '230814-190250'; % '230915-230929_combined';


   load(sprintf('results/%s/Results_%s_%s_%s_img%d_Nsamp%03d_%sConstrained_v1_Cluster_%s.mat',lower(lab_transform),lab_data,lab_noise,lab_power,lab_reso,lab_Nsamp,lab_transform,dateString));
end


%%
ColorMatrix = [0    0.4470    0.7410];


inds = 1:length(numAngles);

numAngles = numAngles(inds);
Nangsamp = numel(numAngles);

indicator = BregDist(inds,1:Nsamp); % 

indicator_sum = sum(indicator,2)/Nsamp;
indicator_sum = indicator_sum';

% Fitting curve
start = 1;
yy = log(indicator_sum(start:end)');
AA = [ones(Nangsamp-start+1,1) log(numAngles(start:end)')];
cc = AA\yy;
comparison = exp(AA*cc);


f = figure;
f.Color = [1 1 1]; % White background
f.Position = [100 100 1280 720];
clf
ph1 = loglog(numAngles(start:end),indicator_sum,'-s','Linewidth',2,'markersize',10);
set(ph1,'Color',ColorMatrix(1,:));
hold on;
loglog(numAngles(start:end),comparison,'k--','Linewidth',1.5);

stand_dev = zeros(1,Nangsamp);
for i = 1:Nangsamp
    tmp = indicator(i,:);
    stand_dev(i) = sqrt(var(tmp));
    ph = loglog([numAngles(i),numAngles(i)],[indicator_sum(i)-stand_dev(i), indicator_sum(i)+stand_dev(i)],'-.','markersize',12);
    set(ph,'Color',ColorMatrix(1,:));
end
 

up = indicator_sum - stand_dev;
down = indicator_sum + stand_dev;
fg = fill([numAngles'; flipud(numAngles')],[down';flipud(up')],[0    0.4470    0.7410],'linestyle','none');
set(fg,'facealpha',.15)

% Legend and axes
xlabel('N','interpreter','latex');
N2 = ['$N^{',num2str(cc(2)),'}$ (fit)'];
% N3 = ['$N^{',num2str(exp_power),'}$ (theo)'];
leg = legend({'Bregman distance',N2},'Location','SW'); %,N3
set(leg,'Interpreter','latex');

% title(['Min error: ',num2str(min(indicator_sum))]);

set(gca,'fontsize',32);
inset = 0.02; % Magic number spacing around the plot
set(gca,'LooseInset',max(inset,get(gca,'TightInset')))

savename = sprintf('%s_%s_%s_%s_Nsamp%03d',lab_data,lab_noise,lab_power,lab_transform,lab_Nsamp);
fprintf('Saving as: %s.png \n',savename)
print('-dpng',fullfile('plots',savename))

