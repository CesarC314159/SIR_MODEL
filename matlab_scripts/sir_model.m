%{
This script solves outputs a graph of the SIR model using Euler's method.

The system of DEs is:
dS/dt=-B*I*S/N
dI/dt=B*I*S/N-GI-MI
dR/dt=GI

S is susceptible
I is infected
R is recovered
N is total population

At any time, S+I+R=N
and dS/dt+dI/dt+dR/dt=0

%}

clear;
N=1000; % Population size
h=.01; % Step size

I(1)=3; % initial infected
S(1)=N-I(1); % initial susceptible
R(1)=0; % intial recovered

B=.4; % this is beta
G=.035; % this is gamma

n=10000; % number of t-values
t_max=h*n; % adjust n and h to adjust t_max

t=0:h:t_max;

for i=1:n
    S(i+1)=S(i)+h*(-1*B*I(i)*S(i))/N;
    I(i+1)=I(i)+h*((B*I(i)*S(i))/N-(G)*I(i));
    R(i+1)=R(i)+h*G*I(i);
end

plot(t,S,t,I,t,R)
legend('S','I','R');
xlabel('time');
ylabel('Individuals');
title('SIR Model with \beta =.4, \gamma=.035');
