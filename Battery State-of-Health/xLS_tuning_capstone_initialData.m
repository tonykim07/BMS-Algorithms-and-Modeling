% Plot the clean and noisy SOC used by this project 
addpath readonly
load readonly/shortDataset.mat
plot(0:199999,100*noisyData, 'b');
hold on; 
plot(0:199999,100*cleanData, 'c', 'linewidth',2);
xlabel('Time (seconds of operation while BMS operational)');
ylabel('State of charge (%)');
title('Noisy SOC estimates used for total-capacity updates (zoom for detail using "xlim")')
legend('Noise-added SOC data','Clean (noise-free) SOC data','location','southwest')
figure;

% Plot the true capacity versus time
load readonly/Qtrue.mat
plot(Qtrue);
xlabel('Time (seconds of operation while BMS operational)');
ylabel('True cell capacity (Ah)');
title('Capacity fading over time for the simulation');
figure;


