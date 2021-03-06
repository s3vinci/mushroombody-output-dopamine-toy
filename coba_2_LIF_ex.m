function [w_stdp,time,V,gi,ge,SM,no_spikes,t_no,sp_trains,i_train,e_train] =  snn_coba_2_LIF_ex(no_runs,duration,w12_0,w21_0,A1,A2,poisson_rate,stdp_active,in_ext_syn_ratio)
tic;
% this network is tuned to produce firing rate at 20 hz when poisson firing
% rate is 20 hz
%% General simulation variables
%Learning Rate
eta=0.001;
% Scaling factor of the excitatory synaptic conductance in units of the
% leak (gleak = 10 nS) vogels
gBarEx=0.014;
% Scaling factor of the inhibitory synaptic conductance in units of the leak (gleak = 10 nS)
gBarIn=0.035;
%gui_simulation;
%handles; % for GUI purposes
% how many runs should the simulation run for

% number of LIF neurons simulated
%% spike monitor
%SM;
%SM(nCells).spikeTimes;

% show how much time elapsed
%% STDP variables
gMax = 2;
% maximum potentiation amplitude
A_plus = A1;
%maximum depression amplitude
A_min = A2;
%% LIF neuron variables
% membrane potential
tau_m = 20;
% resting potential for LIF neuron
V_r  =- 60;
% threshold at which neuron fires
V_th = -50;
% synaptic time constant
tau_syn = 5;
% excitatory synapse potential
tau_e = 5 ;
% inhibitory synapse potential
tau_i = 5;
% number of poisson inhibitory neurons
p_inh = 200;
%number of poisson excitatory neurons
p_ex =  800;
%rate of firing for Poisson neurons
p_rate = poisson_rate;
% time step
dt =0.1;
%membrane voltage of spikes
V_peak = -40 ;
% Excitatory reversal potential, in mV. EAMPA
Ee  = 0;
% Inhibitory reversal potential, in mV. EGABA
Ei = -80;
% excitatory conductance time constant

% excitatory synapse potential
% Scaling factor of the excitatory synaptic conductance in units of the leak (gleak = 10 nS)
g_Bar_Ex = 0.014 ;
% Scaling factor of the inh synaptic conductance in units of the leak (gleak = 10 nS)
g_Bar_In = 0.035;
% potentiation time constant ms
tau_p = 20;
% depression time constant
tau_min = 20;
% STDP weights
%w12_0 = 1;
%w21_0 = 1;
%% initialize vectors
% counts the number of spikes for each neuron
nCells = 2;
SM = struct();
SM(nCells,no_runs).spikeTimes = [];
livePlot = 0;
%% declaration of variables
% synaptic time constant
%total time of simulation ms
totalTime = duration;
% vector that holds each time step
time = 0:dt:totalTime;
% number of poisson excitatory neuronsxw
%%initialize vectors
% membrane potential vector
V = zeros(length(time),nCells);
% excitatory conductance vector
ge = zeros(length(time),nCells);
% inhibitory conductance vector
gi = zeros(length(time),nCells);
% vector that records the number of spikes for each vector
no_spikes = zeros(1,nCells);
%
x_mr = zeros(length(time),nCells);
%
y_mr = zeros(length(time),nCells);

w_stdp = zeros(length(time),nCells);

sp_trains = zeros(length(time),nCells);

% Inhibitory spike trains comprised of 50 Poisson inputs for each neuron
% i_train dims = no_ex_poissson_neurons*no_LIF_neurons x length(time)
i_train = zeros(nCells*p_inh,length(time));
e_train = zeros(nCells*p_ex,length(time));

% vector that stores that strengths for the synapses generated by poisson
% spike trains
synapse_poiss=zeros(nCells,1000);
% vector that stores the strengths of synapses that connect integrate and
% fire neurons
synapse_LIF = zeros(nCells,nCells);
% 1xnoCells column vector 1xnoCells that holds the firing rate of the LIF neurons

% Leak conductance (Everything else is normalized in respect to the leak.)
gLeak=1;
%% Synapse Tuning
for i=1:800
    synapse_poiss(:,i) = 0.4;
end
for i=801:1000
    synapse_poiss(:,i) = 0.6;
end
if nargin == 4
    fprintf('yohoo');
    w12_0;
    w_stdp(1,1) = w12_0;
    w_stdp(1,2) = w21_0;
end
t_no = 1;
no_spikes = zeros(1,nCells);

V = zeros(length(time),nCells);
x_mr = zeros(length(time),nCells);
y_mr = zeros(length(time),nCells);

%initialize spike train
e_train(:,:) = p_rate*dt/1000 >=rand(nCells*p_ex,length(time));
i_train(:,:) = p_rate*dt/1000 >=rand(nCells*p_inh,length(time));
% set the synaptic strengths to their initial values
w_stdp(1,1) = w12_0;
w_stdp(1,2) = w21_0;
sp_trains = zeros(length(time),nCells);
% set the initial voltage to the resting state
V(1,:) = V_r;
% time simulation has been running for
t_running = 0;
% Refractory period for the spike trains.
t_Ref=5;
% (== t(ime)o(f)l(ast)(o)utput (s)pike (Keeping track of the output cell's refractory period )
tolos = zeros(1,nCells);
for  t = 2:length(time)
    t_running = t_running + dt;
      % verify if neuron 1 has reached threshold
    if ((t_running - tolos(1)) < t_Ref)
        %  fprintf('it works - neuron 1');
        V(t,1) = V_r;
    else
        if V(t-1,1) > V_th % a spike has occured
            V(t-1,1) = V_peak;
            % reset the voltage to the resting memebrane potential
            V(t,1) = V_r;
            % update number of spikes for this neuron
            no_spikes(1) = no_spikes(1) + 1;
            % record the time of the spike
            SM(1,t_no).spikeTimes = [SM(1,t_no).spikeTimes t];
            sp_trains(t,1) = 1;
            % ... set the refractory counter to the current time step
            tolos(1)=t_running;
        else
            sp_trains(t,1) = 0;
            V(t,1) = V(t-1,1) +  dt/tau_m*( (V_r-V(t-1,1)+ge(t-1,1)*(Ee-V(t-1,1))) + gi(t-1,1)*(Ei-V(t-1,1 ))) ;
        end
    end
    % verify if neuron 2 has reached threshold
    if ((t_running - tolos(2)) < t_Ref)
        %fprintf('it works - neuron 1');
        V(t,2) = V_r;
    else
        if V(t-1,2) > V_th
            V(t-1,2) = V_peak;
            V(t,2) = V_r;
            no_spikes(2) = no_spikes(2) + 1;
            % record the time of the spike
            SM(2,t_no).spikeTimes = [SM(2,t_no).spikeTimes t];
            sp_trains(t,2) = 1;
            % refractory period set
            tolos(2)=t_running;
        else
            sp_trains(t,2) = 0;
            V(t,2) = V(t-1,2) +  dt/tau_m*( (V_r-V(t-1,2)+ge(t-1,2)*(Ee-V(t-1,2))) + gi(t-1,2)*(Ei-V(t-1,2 ))) ;
             % g_Ex = g_Ex*exp_GEx;
             %g_In = g_In*exp_GIn;
            %  gEx = gEx + gBarEx * sum_SynapseE(i);
            %   gIn = gIn + gBarIn * sum_SynapseI(i);
           
             % Meaning: if the cell is not refractory, ...
            %gTot = gLeak + gEx + gIn;
            % calculate the total membrane conductance,
            %tauEff=taumem/gTot;
            % and the effective time constant, as well as...
           % % the membrane potential that V strives towards.
          %  V(t) = ((gLeak=1*VRest + gEx * EAMPA+ gIn*EGABA)/gTot) + (V(t-1) - ((gLeak*0-VRest + gEx * EAMPA-0 = 0+ gIn*EGABA = -80-0)/gTot))*exp(-dt/(taumem/gLeak + gEx + gIn));
               %  V(t) = VInf + (V(t-1) - VInf)*exp(-dt/tauEff); % is this
               %  perhaps the solution ? ?????? ?? ?? ? ? ? ? ? 
        end
    end
    
     % zeros(nCells*p_inh,length(time));
    i_spikes_index_1 = find(i_train(1:200,t)) + 800; % shift the indices by the number of excitatory cells
   % find index for spikes in excitatory poisson train 
    e_spikes_index_1 = find(e_train(1:800,t));   
    % calculate sum of their contribution to ex / inh conductance
    inh_spikes_weight= sum(synapse_poiss(1,i_spikes_index_1)*g_Bar_In);
    ex_spikes_weight = sum(synapse_poiss(1,e_spikes_index_1)*g_Bar_Ex);
    if stdp_active
        %STDP is calculated
        d_x = sp_trains(t,1);
        d_y = sp_trains(t,2);
        
        x_mr(t,1) = x_mr(t-1,1) + dt*(-x_mr(t-1,1)/tau_p+d_x);
        y_mr(t,1) = y_mr(t-1,1) + dt*(-y_mr(t-1,1)/tau_min+d_y);
        w_stdp(t,1) = w_stdp(t-1,1)+dt*(-A_min*y_mr(t,1)*d_x+A_plus*x_mr(t,1)*d_y);
        
        w_stdp(t,1) = min( w_stdp(t,1),gMax);                        % now ensure no W is above maximum allowed
        w_stdp(t,1)= max(w_stdp(t,1),0);
        %STDP ACTIVE
       % ge(t,1) = ge(t-1,1) - dt/tau_e*ge(t-1,1) + g_Bar_Ex*e_sp_no*0.1+ g_Bar_Ex*w_stdp(t,1)*sp_trains(t,1)*4;
       % LIF excitatory EPSP is 5 times stronger than poisson input
         ge(t,1) = ge(t-1,1) - dt/tau_e*ge(t-1,1) +ex_spikes_weight+ g_Bar_Ex*w_stdp(t,1)*sp_trains(t,2)*in_ext_syn_ratio;
          g_Bar_Ex = g_Bar_Ex+(w_stdp(t,1) - w_stdp(t-1,1));
    else
           ge(t,1) = ge(t-1,1) - dt/tau_e*ge(t-1,1) + ex_spikes_weight; %0.0495
    end
    %gi(t,1) = gi(t-1,1) - dt/tau_i*gi(t-1,1)+ g_Bar_In*i_sp_no*0.35; %0.1 %j.nCells;  % for each input
     gi(t,1) = gi(t-1,1) - dt/tau_i*gi(t-1,1)+ inh_spikes_weight; %0.1 %j.nCells;  % for each input
   % gi(t,1) = 300;
    i_spikes_index_2 = find(i_train(201:end,t)) + 800; % shift the indices by the number of excitatory cells
   % find index for spikes in excitatory poisson train 
    e_spikes_index_2 = find(e_train(801:end,t));   
    % calculate sum of their contribution to ex / inh conductance
    %i_spikes_index_2
    %synapse_poiss(2,i_spikes_index_2)
    inh_spikes_weight2= sum(synapse_poiss(2,i_spikes_index_2)*g_Bar_In);
    ex_spikes_weight2 = sum(synapse_poiss(2,e_spikes_index_2)*g_Bar_Ex);
    if stdp_active
        %STDP is calculated
        d_x = sp_trains(t,2);
        d_y = sp_trains(t,1);
        
        x_mr(t,2) = x_mr(t-1,2) + dt*(-x_mr(t-1,2)/tau_p+d_x);
        y_mr(t,2) = y_mr(t-1,2) + dt*(-y_mr(t-1,2)/tau_min+d_y);
        w_stdp(t,2) = w_stdp(t-1,2)+dt*(-A_min*y_mr(t,2)*d_x+A_plus*x_mr(t,2)*d_y);
        
        w_stdp(t,2) = min( w_stdp(t,2),gMax);                        % now ensure no W is above maximum allowed
        w_stdp(t,2)= max(w_stdp(t,2),0);

        %ge(t,2) = ge(t-1,2) - dt/tau_e*ge(t-1,2) + g_Bar_Ex*e_sp_no*0.1 + g_Bar_Ex*w_stdp(t,2)*sp_trains(t,2)*4;
        ge(t,2) = ge(t-1,2) - dt/tau_e*ge(t-1,2) +ex_spikes_weight2 + g_Bar_Ex*w_stdp(t,2)*sp_trains(t,1)*in_ext_syn_ratio;
        g_Bar_Ex = g_Bar_Ex+(w_stdp(t,2) - w_stdp(t-1,2));
    else
        ge(t,2) = ge(t-1,2) - dt/tau_e*ge(t-1,2) + ex_spikes_weight2;
    end
   gi(t,2) = gi(t-1,2) - dt/tau_i*gi(t-1,2)+ inh_spikes_weight2; %j.nCells;  % for each input
    
    %gi(t,2) = 300;
        
end
timeElapsed = toc;
fprintf('1 second simulating 2 neurons with STDP and Poisson input took %f\n',timeElapsed);
%plotter.plotWeights(time,(w_stdp(:,1,t_no)),(w_stdp(:,2,t_no)))
t_no = t_no + 1;
end
%end