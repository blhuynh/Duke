clear; clc;

%% -------- PART 2: INTRINSICALLY BURSTING IZHIKEVICH NEURON --------
% Initialize Parameters
% Time (ms)
dt = 0.1;
tstop = 1000;
tvec = 0:dt:tstop;
delay = 5;

% Istim
Istim = zeros(1,length(tvec));
Istim(5/dt+1:end) = 8;
Istim((tstop-100)/dt+1:end)= 0;

% ODE parameters
a = 0.02; % time scale of recovery variable u
b = 0.2; % sensitivity of recovery variable u to subthreshold fluctuations of membrane potential v
c = -55; % after-spike reset value of membrane potential v caused by fast high-threshold K+ conductances
d = 4; % after-spike reset of recovery variable u caused by slow high-threshold Na+ and K+ conductances
v0 = -70;
u = b * v0;

% Iterate through time loop
v = [v0 zeros(1,length(tvec)-1)];
for k = 2:length(tvec)
    I = Istim(k);
    
    v(k) = v(k-1) + dt*(0.04*(v(k-1)^2)+5*v(k-1)+140-u+I);
    u = u + dt*a*(b*v(k-1) - u);
    
    if v(k) >= 30
        v(k) = c;
        u = u + d;
    end
end

% Plot Figure
figure(1); clf; hold on
subplot(2,1,1); hold on
plot(tvec,v)
xlabel('Time (ms)'); ylabel('v(t)')
title('Intrinsically Bursting Izhikevich Neuron')
setfont

subplot(2,1,2); hold on
plot(tvec,Istim)
xlabel('Time (ms)'); ylabel('I_{stim}(t)')
setfont

print -dpng bme515_hw1_part2.png
close(1);

%% ------------------------ PART 3: SYNPASE ------------------------
if true
    pmax = 1; % peak value
    taus = 5; % (ms)
    ps = @(t) pmax*t.*exp(1-t/taus)/taus;
    
    figure(2); clf; hold on
    plot(tvec,ps(tvec))
    xlabel('Time (ms)'); ylabel('P_s(t)')
    title('Open Postsynaptic Channel Probability')
    setfont
    
    print -dpng bme515_hw1_part3_1.png
    close(2)
    
    % series of spikes
    synapses = [ps(tvec-1); ps(tvec-5); ps(tvec-9); ps(tvec-15)];
    synapses(synapses<0) = 0;
    
    figure(3); clf; hold on
    subplot(2,1,1); hold on
    plot(tvec,synapses(1,:),'k-')
    plot(tvec,synapses(2,:),'k--')
    plot(tvec,synapses(3,:),'k-.')
    plot(tvec,synapses(4,:),'k:')
    axis([0 tstop 0 pmax])
    legend('t=1','t=5','t=9','t=15')
    xlabel('Time (ms)'); ylabel('P_s(t)')
    title('Probability of Open Postsynaptic Channel')
    setfont
    
    subplot(2,1,2); hold on
    plot(tvec,sum(synapses),'k-')
    xlabel('Time (ms)'); ylabel('P_s(t)')
    setfont
    print -dpng bme515_hw1_part3_2.png
    
    close(3)
end

%% ---------------- PART 4: TWO-NEURON OSCILLATING NETWORK ----------------
u = b * c * ones (1,2);
% Iterate through time loop
v = repmat([c zeros(1,length(tvec)-1)],[1 2]);
ps = @(t,pmax,taus) pmax*t.*exp(1-t/taus)/taus;

count = 1;
for var2 = 4
    for var3 = 3
        pmax = [var2 1];
        taus = [var3 5];
        Ps = [ps(tvec,pmax(1),taus(1)); ps(tvec,pmax(2),taus(2))];
        Ps(Ps<0) = 0;
        
        for var1=19
            clear v u Istim
            v = [repmat(v0,[1,2]); zeros(numel(tvec)-1,2)];
            u = ones(1,2) * (b * v0);
            
            Ps = [ps(tvec,pmax(1),taus(1)); ps(tvec,pmax(2),taus(2))];
            Ps(Ps<0) = 0;
            
            Istim = zeros(2,length(tvec));
            Istim(1,5/dt+1:end) = 8;
            Istim(2,(5+var1)/dt+1:end) = 8;
            Istim(:,(tstop-25)/dt+1:end) = 0;
            
            for k = 2:length(tvec)
                v(k,1) = v(k-1,1) + dt*(0.04*v(k-1,1)^2+5*v(k-1,1)+140-u(1)+Istim(1,k)-Ps(1,k));
                v(k,2) = v(k-1,2) + dt*(0.04*v(k-1,2)^2+5*v(k-1,2)+140-u(2)+Istim(2,k)-Ps(2,k));
                
                u(1) = u(1) + dt*a*(b*v(k-1,1) - u(1));
                u(2) = u(2) + dt*a*(b*v(k-1,2) - u(2));
                
                % AP spike on Neuron 1
                if v(k,1) >= 30
                    v(k,1) = c;
                    u(1) = u(1) + d;
                    newPS = ps(tvec-tvec(k),pmax(1),taus(1));
                    newPS(newPS<0) = 0;
                    Ps(1,k:end) = Ps(1,k:end) + newPS(1:length(Ps(1,k:end)));
                    clear newPS
                end
                
                % AP spike on Neuron 2
                if v(k,2) >= 30
                    v(k,2) = c;
                    u(2) = u(2) + d;
                    newPS = ps(tvec-tvec(k),pmax(2),taus(2));
                    newPS(newPS<0) = 0;
                    Ps(2,k:end) = Ps(2,k:end) + newPS(1:length(Ps(2,k:end)));
                    clear newPS
                end
                
                %     if mod(k,25)==0; fprintf([num2str(k),'\n']); end
            end
            
            % Plot Figure
            figure(4); clf; hold on
            subplot(3,1,1); hold on
            plot(tvec,v(:,1),'k')
            plot(tvec,v(:,2),'b')
            %         vline(30)
            %         vline(122)
            xlabel('Time (ms)'); ylabel('v(t)')
            title(['Two-Neuron Oscillating Network'])
            setfont; axis tight;
            
            subplot(3,1,2); hold on
            plot(tvec,Istim)
            xlabel('Time (ms)'); ylabel('I_{stim}(t)')
            setfont; axis([0 max(tvec) 0 1.1*max(max(Istim))]);
            
            subplot(3,1,3); hold on
            plot(tvec,Ps(1,:),'k')
            plot(tvec,Ps(2,:),'b')
            xlabel('Time (ms)'); ylabel('P_s(t)')
            legend('Neuron 1','Neuron 2')
            setfont; axis tight;
            print -dpng bme515_hw1_part4.png
            
            %             fprintf(sprintf('%g %g %g \n',var1, var2, var3))
            %             keyboard
            %% -------------- Part 5: QUANTIFICATION OF THE OSCILLATION ---------------
            timeBlock = [25/dt : (tstop-25)/dt];
            [pks1,locs1] = findpeaks(v(timeBlock,1),'minpeakheight',0);
            [pks2,locs2] = findpeaks(v(timeBlock,2),'minpeakheight',0);
            
            for k = 1:floor(length(locs1)/2)-1
                period1(k) = (timeBlock(locs1(2*(k+1)-1))-timeBlock(locs1(2*k-1))) * dt *1e-3; % ms -> s
                %     fprintf(sprintf('Neuron 1: period of %g ms\n',period1(k)))
            end
            
            for k = 1:floor(length(locs2)/2)-1
                period2(k) = (timeBlock(locs2(2*(k+1)-1))-timeBlock(locs2(2*k-1))) * dt * 1e-3; % ms -> s
                %     fprintf(sprintf('Neuron 2: period of %g ms\n',period2(k)))
            end
            
            % fprintf(sprintf('Neuron 1: mean period of %g ms\n',mean(period1(:))))
            % fprintf(sprintf('Neuron 2: mean period of %g ms\n',mean(period2(:))))
            
            
            % Phase Offset
            avgFreq = 1 ./ [mean(period1*1e-3) mean(period2*1e-3)];
            timeloc = [1:length(period2)];
            for k = 1:length(timeloc)
                timeDifference(k) = [locs1(timeloc(k)) - locs2(timeloc(k))] * dt * 1e-3;
            end
            
            for k = 1:length(timeloc)
                offset(k) = [timeDifference(k)/period2(k)] * 360;
            end
            
            fprintf(sprintf('%2.2f \n', sum(sqrt((offset-180).^2)) ))
            cost(count,:) = [var1 var2 var3 sum(sqrt((offset-180).^2)) ];
            count=count+1;
        end
    end
end



