close force all; clc; clear all;
tic;

scale = 1;
p_ex = 10;
p_rate = 20;
dt=0.1;
duration= 10000;

time=0: dt : duration;
% LIF model params for neurons
Vrest=-70;
Vex=0;
Rm=10;
Vth=-54;
Vthpn = -60;
Vthkc = -40;
Vthon = 0;

Vthda = -30;


Vreset=-80;
cm=10;  % 10nF/mm^2
tauex=10;
Gleak=1.0; %1uS/mm^2
dGex=0.5; %0.5uS/mm^2
dGinh=0.1; %0.5uS/mm^2
exp_fac=exp(-dt/tauex);

w_exp_fac=exp(-dt/5);
%number of pn neurons
pnn = 150;
%number of orn neurons
ornn = 1000;
% number of KC neurons
kcn = 800;
%number of dopaminergic neurons
dan = 1;
% initializing LIF model params for orn neurons
pp.orn.V=zeros(ornn,1);
pp.orn.w=zeros(ornn,1);
pp.orn.Gex=zeros(ornn,1);
pp.orn.Ginh=zeros(ornn,1);
pp.orn.V(:)=Vrest;
% initializing LIF model params for PN neurons
pp.pn.V=zeros(pnn,1);
pp.pn.Vinf= zeros(pnn,1) ;
pp.pn.Gex=zeros(pnn,1);
pp.pn.Ginh=zeros(pnn,1);
pp.pn.V(:)=Vrest;
% initializing LIF model params for KC neurons
pp.kc.V= zeros(kcn,1);
pp.kc.Vinf= zeros(kcn,1);
pp.kc.Gex=zeros(kcn,1);
pp.kc.Ginh=zeros(kcn,1);
pp.kc.V(:)=Vrest;
pp.kc.gex(:)=0;
% initializing LIF model params for single ON neuron
pp.on.V= zeros(1,length(time));
pp.on.Vinf= 0;
pp.on.V(:)=Vrest;
pp.on.Gex=zeros(1,1);
pp.on.Ginh=zeros(1,1);
% initializing LIF model params for single DA neuron
pp.da.V= zeros(1,length(time));
pp.da.Vinf= 0;
pp.da.V(:)=Vrest;
pp.da.Gex=zeros(1,1);
pp.da.Ginh=zeros(1,1);

spikeTimes_orn = cell(ornn,1);
spikeTimes_pn = cell(pnn,1);
spikeTimes_kc = cell(kcn,1);

spikeTimes_on_da = cell(2,1);
%spikeTimes_da = cell(1,1);

ornperpn = randsample([6 5],150,true);
wornpn =             false(ornn, pnn);

Iex=zeros(1,length(time));
Iinh=zeros(1,length(time));

[o,p]= loadonpn();
orn_rates = randsample(o(20,:),1000,'true');
orn_rates = orn_rates';
%orn_rates(1:200,:) = 0;
spikes = o(10,:);

%num_kcs = length(pns_per_kc);
connectivity = cell(150, 1);

for i=1:150
    connectivity{i} = randsample(ornn,ornperpn(i),true);
end
for i=1:150
    % assign connections of random weight (most will be empty)
    
    % connectivity{i} is the array of PNs indices assigned to each KC
    wornpn(connectivity{i}, i) = 1;
end
wornpn = normc(wornpn);
%wornpn = wornpn * 5;
pnperkc = randsample([7 5 6],800,true);
wpnkc =             false(pnn, kcn);
connectivity = cell(800, 1);
for i=1:800
    connectivity{i} = randsample(pnn,pnperkc(i),true);
    
end
for i=1:800
    % assign connections of random weight (most will be empty)
    % connectivity{i} is the array of PNs indices assigned to each KC
    wpnkc(connectivity{i}, i) = 1;
end
wpnkc = normc(wpnkc);
wkcon =             ones(kcn,1);


% STDP weights

A_plus = 0.5;
A_minus = 0.5;
x_mr = zeros(length(time),kcn);
y_mr = zeros(length(time),kcn);
w_stdp =  zeros(length(time),kcn);
for i=1:length(time)-1
    if mod(i,2000) == 0
        i
        figure(1)
        plotSpikeRaster(spikeTimes_orn,'PlotType','vertline');
        ylabel('Orn No')
        title('ORN Raster Plot');
        set(gca,'XTick',[]);
        figure(2);
        plotSpikeRaster(spikeTimes_pn,'PlotType','vertline');
        ylabel('PN No')
        title('PN Raster Plot');
        set(gca,'XTick',[]);
        figure(3);
        plotSpikeRaster(spikeTimes_kc,'PlotType','vertline');
        ylabel('KC No')
        title('KC Raster Plot');
        set(gca,'XTick',[]);
        
        figure(4); plot(pp.on.V(1,:))
        
        figure(5); plot(pp.da.V(1,:),'r')
        title('DA Neuron V');
        
        figure(10)
        plotSpikeRaster(spikeTimes_on_da,'PlotType','vertline');
        disp('DA neuron spiked x times in 2 seconds');
        disp('ON neuron spiked x times in 2 seconds');
        drawnow;
    end
    
    ex_train = orn_rates*dt/1000 >=rand(ornn,1);
    inh_train = p_rate*dt/1000 >=rand(ornn,1);
    
    ex_train_kc = orn_rates*dt/1000 >=rand(ornn,1);
    inh_train_kc = p_rate*dt/1000 >=rand(ornn,1);
    
    ex_train_da = 40*dt/1000 >=rand;
    inh_train_da = p_rate*dt/1000 >=rand;
    
    % pp.kc.Gex(:,i+1)=ex_train(:,i+1)*dGex + pp.kc.Gex(:,i)*exp_fac;
    %   pp.kc.Ginh(:,i+1)=inh_train(:,i+1)*dGinh + pp.kc.Ginh(:,i)*exp_fac;
    %    myce;
    
    
    pp.orn.Gex(:)=ex_train*dGex + pp.orn.Gex(:)*exp_fac;
    pp.orn.Ginh(:)=inh_train*dGinh + pp.orn.Ginh(:)*exp_fac;
    pp.orn.Gtotal=Gleak+pp.orn.Gex(:) - pp.orn.Ginh(:);
    taueff=(cm) ./(pp.orn.Gtotal);
    pp.orn.Vinf=(Gleak*Vrest+pp.orn.Gex(:)*Vex ) ./ (pp.orn.Gtotal);
    pp.orn.V(:)=pp.orn.Vinf(:)+(pp.orn.V(:)-pp.orn.Vinf ).*exp(-dt ./ taueff);
    %   pp.kc.V(   find(  pp.kc.V(:)>Vth ) )= 50;
    orn_logical = pp.orn.V(:) >= Vth;
    sp =   (  find(  pp.orn.V(:)>Vth )) ;
    for indsp = 1:length(sp)
        spikeTimes_orn{sp(indsp)} = [ spikeTimes_orn{sp(indsp)} i];
    end
    pp.orn.w(find ( pp.orn.V(:)>Vth)) =    pp.orn.w(find ( pp.orn.V(:)>Vth)) + 5;
    % pp.orn.w = pp.orn.w * w_exp_fac;
    
    pp.orn.V( find (pp.orn.V(:)>Vth )    )=Vreset;
    pn_gex = orn_logical' * wornpn;
    pp.pn.Gex(:)=pn_gex'*dGex + pp.pn.Gex(:)*exp_fac;
    pp.pn.Gtotal=Gleak+pp.pn.Gex(:);
    taueff=(cm) ./(pp.pn.Gtotal);
    pp.pn.Vinf=(Gleak*Vrest+pp.pn.Gex(:)*Vex) ./ (pp.pn.Gtotal);
    pp.pn.V(:)=pp.pn.Vinf(:)+(pp.pn.V(:)-pp.pn.Vinf ).*exp(-dt ./ taueff);
    sp =   (  find(  pp.pn.V(:)>Vthpn )) ;
    for indsp = 1:length(sp)
        spikeTimes_pn{sp(indsp)} = [ spikeTimes_pn{sp(indsp)} i];
    end
    pn_logical =  pp.pn.V(:)>Vthpn;
    
    pp.pn.V( find (pp.pn.V(:)>Vthpn )    )=Vreset;
    
    % sum(find(pn_logical))
    %pn_logical=  randsample([0,0,0,1,0,0,0,00],150,'true');
    %pn_logical = zeros(150,1);
    %pn_logical(1:10,:) = 1;
    kc_gex = pn_logical' * wpnkc;
    pp.kc.Gex(:)=kc_gex'*dGex + pp.kc.Gex(:)*exp_fac;
    pp.kc.Gtotal=Gleak+pp.kc.Gex(:);
    taueff=(cm) ./(pp.kc.Gtotal);
    pp.kc.Vinf=(Gleak*Vrest+pp.kc.Gex(:)*Vex) ./ (pp.kc.Gtotal);
    pp.kc.V(:)=pp.kc.Vinf(:)+(pp.kc.V(:)-pp.kc.Vinf ).*exp(-dt ./ taueff);
    sp =   (  find(  pp.kc.V(:)>Vthkc )) ;
    for indsp = 1:length(sp)
        spikeTimes_kc{sp(indsp)} = [ spikeTimes_kc{sp(indsp)} i];
    end
    kc_logical =  pp.kc.V(:)>Vthkc;
    pp.kc.V( find (pp.kc.V(:)>Vthkc )    )=Vreset;
    
    on_gex = kc_logical' * wkcon;
    pp.on.Gex(:)=on_gex'*dGex + pp.on.Gex(:)*exp_fac;
    pp.on.Gtotal=Gleak+pp.on.Gex(:);
    taueff=(cm) ./(pp.on.Gtotal);
    pp.on.Vinf=(Gleak*Vrest+pp.on.Gex(:)*Vex) ./ (pp.on.Gtotal);
    pp.on.V(:,i+1)=pp.on.Vinf(:)+(pp.on.V(:,i)-pp.on.Vinf ).*exp(-dt ./ taueff);
    if pp.on.V(:,i)>Vthon
        pp.on.V(  1, i       )= 50;
        pp.on.V( 1, i+1     )=Vreset;
        spikeTimes_on_da{1} = [ spikeTimes_on_da{1} i];
    end
    %pp.on.V(   find(  pp.on.V(:,i)>Vthon ) ,i       )= 50;
    %pp.on.V( find (pp.on.V(:,i)>Vthon ), i+1     )=Vreset;
    %calculating membrane potential for DA neuron
    %keyboard;
    on_logical =  pp.on.V(  1, i       )== 50;
    
    on_gex =on_logical * 100;
      
    pp.da.Gex(:)=ex_train_da*dGex + on_gex + pp.da.Gex(:)*exp_fac;
    
    pp.da.Gtotal=Gleak+pp.da.Gex(:);
    taueff=(cm) ./(pp.da.Gtotal);
    pp.da.Vinf=(Gleak*Vrest+pp.da.Gex(:)*Vex) ./ (pp.da.Gtotal);
    pp.da.V(:,i+1)=pp.da.Vinf(:)+(pp.da.V(:,i)-pp.da.Vinf ).*exp(-dt ./ taueff);
    %pp.da.V(i)
    if pp.da.V(:,i)>Vthda
        pp.da.V(  1, i       )= 50;
        pp.da.V( 1, i+1     )=Vreset;
      %  keyboard;
        spikeTimes_on_da{2} = [ spikeTimes_on_da{2} i];
    end
    %pp.da.V(   find(  pp.da.V(:,i)>Vthon ) ,i)= 50;
    %pp.da.V( find (pp.da.V(:,i)>Vthon ), i+1)=Vreset;  
    x_mr(t,:) = x_mr(t-1,:) + dt*(-x_mr(t-1,:)/tau_p+d_x);
    y_mr(t,:) = y_mr(t-1,:) + dt*(-y_mr(t-1,:)/tau_min+d_y);
    w_stdp(t,:) = w_stdp(t-1,:)+dt*(-A_min*y_mr(t,:)*d_x+A_plus*x_mr(t,:)*d_y);
end
%mean ( sum(ex_train,2) )
%figure; plot(pp.kc.V(1,:))
t1 = toc
tic;
t2=     toc
