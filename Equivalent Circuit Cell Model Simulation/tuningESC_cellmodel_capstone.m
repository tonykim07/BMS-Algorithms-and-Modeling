addpath readonly 

load readonly/pulseData.mat
load readonly/pulseModel.mat

T = 25; % Test Temperature
deltaT = 1; % sampling period for data
Q = getParamESC('QParam', T, model); 

tk = pulseData.time;    % testing time
ik = pulseData.current; % testing current
vk = pulseData.voltage; % testing voltage
rcValues = tuneModel;

R0 = rcValues(1)/1000; % convert milliohms to ohms
R = rcValues(2:2:end)/1000; % convert these also
C = rcValues(3:2:end)*1000; % convert kF to F
RCfact = exp(-deltaT./(R.*C));

% Simulate the dynamic states of the model
iRk = zeros(length(RCfact),1); % initial resistor currents
vCk = 0*pulseData.voltage; % initialize capacitor voltages
if ~isempty(RCfact)
  for k = 2:length(vCk)
    iRk = diag(RCfact)*iRk + (1-RCfact)*ik(k-1); % update resistor current
    vCk(k) = R'*iRk; % compute capacitor voltage
  end
end

% Simulate SOC state
z0 = SOCfromOCVtemp(pulseData.voltage(1),25,model); 
zk = z0-cumsum([0;ik(1:end-1)])*deltaT/(Q*3600); 

% Compute voltage estimate
vest = OCVfromSOCtemp(zk,25,model) - vCk - ik.*R0;

% Compare against measured voltage, compute RMSE in mV
rmse = 1000*sqrt(mean((vk - vest).^2));

% Plot some results 
subplot(2,1,1); plot(tk/60,vk,tk/60,vest); 
title('Voltage estimation'); grid on
xlabel('Time (min)'); ylabel('Voltage (V))'); 
legend('True measured voltage','Model voltage','location','southeast');

subplot(2,1,2); plot(tk/60,1000*(vk-vest)); 
title('Voltage estimation errors');
xlabel('Time (min)'); ylabel('Voltage error (mV)');
grid on

% output message for results
fprintf('Your tuning values produced an RMS voltage-prediction error of %g mV\n',rmse);

% Function that specifies the resistor and capacitor values for Extended
% Thevenin Model 

function rcValues = tuneModel
% [R0; R1; C1; R2; C2]
    rcValues = [9.4; 11.7; 49; 5.9; 3.5];
end
