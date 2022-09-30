addpath readonly
load ./readonly/E2model; 

load ./readonly/packData
% Plot histogram of cell SOC in this pack
histogram(100*packData.storez); 
title('Histogram of battery-pack individual SOC values before balancing');
xlabel('State-of-charge (%)'); ylabel('Count (in 100-cell battery pack)');
xlim([85 95])

% Results of this code
rBalance = tuneBalancer;
[storez_output,pk,maxpk,tbal,saveV] = cellBalance(packData,model,rBalance);

t = 1:length(pk);
subplot(2,2,1);
plot(t/3600,pk); xlabel('Time (h)'); ylabel('Total pack balancing power'); 
title('Total pack power versus time');

subplot(2,2,2)
plot(t/3600,maxpk); xlabel('Time (h)'); ylabel('Max. cell balancing power'); 
title('Max. cell balancing power versus time');

subplot(2,2,3)
histogram(100*storez_output); % Plot histogram of post-balancing cell SOC in this pack (in percent)
title('Post-balancing histogram of SOC values'); xlabel('State-of-charge (%)');
ylabel('Count (in 100-cell battery pack)'); xlim([85 95]);

subplot(2,2,4)
plot(t/3600,saveV'); % Plot trace of all cell voltages over time
title('Cell voltages during balancing'); xlabel('Time (h)'); ylabel('Voltage (V)');

fprintf('During balancing, maximum pack balancing power was %fW (should be less than 10W)\n',max(pk));
fprintf('During balancing, maximum cell balancing power was %fW (should be less than 0.1W)\n',max(maxpk));
fprintf('Time to balance was %fh (should be less than 29.5h)\n',tbal/3600);

% Balancing simulation functions
function [storez_output, pk, maxpk, tbal, saveV] = cellBalance(packData, model, Rbal)

    maxtime = 36*3600;
    pk = zeros(maxtime,1);
    maxpk = zeros(maxtime,1);
    saveV = zeros(length(packData.storez),maxtime);

    % Initialize states for ESC cell model
    z = packData.storez;
    irc = packData.storeirc; 

    % Default initialization for cells within the pack 
    T = packData.T;       % Temperature for each cell, assumed to be constant
    Tsd = packData.Tsd;   % Self-discharge "temperature" for each cell
    leak = packData.leak; % Leakage current for each cell
    q = packData.q;       % Total capacity for each cell (all may be different)
    rc = packData.rc;     % R-C time constant for each cell
    r = packData.r;       % Diffusion-resistor values for each cell
    r0 = packData.r0;     % Series resistance values for each cell

    % Run actual simulations
    for k = 1: maxtime
        % Calculate cell voltages
        v = OCVfromSOCtemp(z,T,model); % get OCV for each cell
        v = v - r.*irc; % add in capacitor voltages (ignore hysteresis to simplify sim)
        saveV(:,k) = v(:); % save all cells' terminal voltage for later analysis
        % Cell Simulation
        ik = zeros(size(v)); % no string current
        % Simulate self discharge via variable resistor in parallel
        rsd = ((-20+0.4*Tsd).*z + (35-0.5*Tsd))*1e3; ik = ik + v./rsd;
        % Simulate leakage current
        ik = ik + leak;
        % Check to see which cells are 2mV or more above minimum cell voltage
        checkBalance = (v - min(v)) - 2e-3 >= 0;
        if sum(checkBalance) == 0, % balancing is complete, so return
          saveV = saveV(:,1:k-1);
          pk = pk(1:k-1);
          maxpk = maxpk(1:k-1);
          tbal = k-1;
          storez_output = z; % Output only final SOC for each cell after 4 hours of simulation time
          return
        end
        % cells 2mV or more from the minimum will have array value 1, otherwise 0
        % Add balancing resistors and calculate resulting cell current
        v_balance = v.*checkBalance;
        % Set non-balance cell voltage to 0 for calculation (to ensure balance current = 0 for no-balance cells)
        i_balance = (v_balance./Rbal); % Current calculated for balance cell, with parallel resistor 
        ik = ik + i_balance; % Add balance current to externally applied cell current
        
        % Compute power
        pk(k) = sum(i_balance.^2*Rbal);    % total power dissipated by all cells being balanced in pack 
        maxpk(k) = max(i_balance.^2*Rbal); % maximum single-cell power dissipated
        % Calculate new SOC for each cell
        z = z - (1/3600)*ik./q; % Update each cell SOC
        % Update diffusion-resistor currents
        irc = rc.*irc + (1-rc).*ik; % Update capacitor voltages
        if mod(k,3600) == 0
          fprintf('Completed %dh balancing (%d "unbalanced" cells remain).\n',round(k/3600),sum(checkBalance));
        end
    end
    % If we get to this point, then must have balanced for more than 200h (not good)
    tbal = k;
    storez_output = z; % Output only final SOC for each cell 
end
function rBalance = tuneBalancer
    rBalance = 170; % hand tuned value
end



