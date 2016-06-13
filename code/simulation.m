% Part 1: initializing 
clear;

% initializing time steps for simulation
numSteps = 500;

% initializing frequency for square function
deltaT = 0.3;
%setting limit
minLimx = 0.1;
maxLimY = 1.5;
% parameters
beta = 0.5;
alpha = 0.5;
gamma = 0.5;
zeta = 0.5;
thetaFs = 0.5;
thetaFl = 0.5;
eta = 0.5;
psi = 0.5;
omegaA1 = 0.3;
omegaA2 = 0.3;
omegaA3 = 0.3;
epsilon = 0.5;

%Part 2: declaring vectors for all nodes and set initial value to external factors
Op = zeros(1,numSteps); %Observing penalty
Sp = zeros(1,numSteps); %Social support
Pe = zeros(1,numSteps); %Uncomfortable physical environment
Sv = zeros(1,numSteps); %Stressfull event
Np = zeros(1,numSteps); %Negative personality
Ac = zeros(1,numSteps); %Aggressive cues

%declaring intial values for temporal properties
Fl(1) = 0.3;
Rs(1) = 0.1;
As(1) = 0.3;

%setting values to external properties
for t=1:numSteps
    Op(t) = 0.9; %Observing penalty
    Sp(t) = 0.9; %Social support
    Pe(t) = 0.9; %Uncomfortable physical environment
    Sv(t) = 0.9; %Stressfull event
    Np(t) = 0.9; %Negative personality
    Ac(t) = 0.9; %Aggressive cues
end

%initializing instantaneous functions at time t=1
Aw(1) = beta * Op(1); %awareness
Ds(1) = (alpha * Sv(1) + (1-alpha) * Pe(1))*(1-Sp(1)*Hd(1)); %discomfort 
Hd(1) = (zeta*Aw(1)+(1-zeta)*(1-Np(1)))*(1-Sc(1)); %hardiness 
Sc(1) = (gamma*Ds(1)+(1-gamma)*Np(1))*(1-Hd(1)); %susceptibility
Ir(1) = Ds(1)*(1-Rs(1)); %irritation
Fs(1) = thetaFs * Ds(1) + (1-thetaFs)*Ir(1); %short term frustration
Ss(1) = eta*Ir(1)+(1-eta)*Sc(1); %stress
Ag(1) = (omegaA1 * Fl(1) + omegaA2 * Sc(1) + omegaA3 * Ss(1))*(1-Rs(1)); %anger

%Part 3: equations
for t = 2:numSteps
    Aw(t) = beta * Op(t); %awareness
    Ds(t) = (alpha * Sv(t) + (1-alpha) * Pe(t))*(1-Sp(t)*Hd(t)); %discomfort 
    Hd(t) = (zeta*Aw(t)+(1-zeta)*(1-Np(t)))*(1-Sc(t)); %hardiness 
    Sc(t) = (gamma*Ds(t)+(1-gamma)*Np(t))*(1-Hd(t)); %susceptibility
    
    
    
    Rs(t) = Rs(t-1) + psi *(Hd(t) - Rs(t-1)) * (1-Rs(t-1))*Rs(t-1)*deltaT; %resilience
    Ir(t) = Ds(t)*(1-Rs(t)); %irritation
    Fs(t) = thetaFs * Ds(t) + (1-thetaFs)*Ir(t); %short term frustration
    Ss(t) = eta*Ir(t)+(1-eta)*Sc(t); %stress
    Fl(t) = Fl(t-1) + thetaFl * (Fs(t) - Fl(t-1))*(1-Fl(t-1))*Fl(t-1)*deltaT; %long-term frustration
    Ag(t) = (omegaA1 * Fl(t) + omegaA2 * Sc(t) + omegaA3 * Ss(t))*(1-Rs(t)); %anger
    
    
    As(t) = As(t-1) + epsilon * (Ag(t) - As(t-1)) * Ac(t) * (1-As(t-1)) * As(t-1) * deltaT; %agressiveness
end
    %Part 4: plotting
    hold on
    t=1:numSteps;
    
    subplot(2,1,1);
    y = plot (t, As, 'r:');
    
    xlabel('time steps'); ylabel('Aggresiveness');
    xlim([0 numSteps]);ylim([minLimx maxLimY]);
    hold off;
    %legend(y,'Aggresiveness');
    











