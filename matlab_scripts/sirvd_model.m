%{
This script solves the SIRD model using Euler's method.

The system of DEs is:
dS/dt=-B*I*S/N
dI/dt=B*I*S/N-GI-MI
dR/dt=GI
dV/dt=US
dD/dt=MI

S is susceptible
I is infected
R is recovered (or removed)
V is vaccinated
D is dead
N is total population

At any time, S+I+R+V+D=N
and dS/dt+dI/dt+dR/dt+dV/dt+dD/dt=0

%}

clear;
N=1000; % Population size
h=.01; % Step size

I(1)=3; % initial infected
S(1)=N-I(1); % initial susceptible
D(1)=0; % intial dead
V(1)=0; % initial vaccinated
R(1)=0; % intial recovered

B=.4; % this is beta, 
G=.035; % this is gamma
M=.005; % this is mu
U=.010; % this is nu

n=10000; % number of t-values
t_max=h*n; % adjust n and h to adjust t_max

t=0:h:t_max;

for i=1:n
    S(i+1)=S(i)+h*(-1*B*I(i)*S(i))/N;
    I(i+1)=I(i)+h*((B*I(i)*S(i))/N-(G+M)*I(i));
    R(i+1)=R(i)+h*G*I(i);
    V(i+1)=V(i)+h*U*S(i);
    D(i+1)=D(i)+h*M*I(i);
end

plot(t,S,t,I,t,R,t,V,t,D)
legend('S','I','R','V','D');
xlabel('time');
ylabel('Individuals');
title('SIRVD Model with \beta =.4, \gamma=.035, \mu=.005, \nu=.010');
