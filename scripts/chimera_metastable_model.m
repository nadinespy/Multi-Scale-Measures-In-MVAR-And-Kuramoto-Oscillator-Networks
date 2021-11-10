function [X, sigma_chi, synchrony] = chimera_metastable_model(coupling_matrix, npoints, beta, intra_comm_size, n_communities)
% Produces metastable chimera-like states in a modular network of
% oscillators. The model is inspired by Abrams, et al., PRL 2008, and is
% parameterised by b and A, according to their paper.
% For more details of this system, see:
%
% Shanahan, M. (2010). Metastable Chimera States in Community-Structured
% Oscillator Networks. Chaos 20, 013108.
%
% Written by Murray Shanahan, August 2009 / December 2014

b = beta;				% was 0.1
a = pi/2-b;				% phase lag (alpha): 1.4708

n0 = intra_comm_size;				     % intra-community size
n1 = n_communities;					% number of communities
d0 = intra_comm_size; 
d1 = intra_comm_size;				      % numbers of connections at different community levels
		
N = n_communities*intra_comm_size;   % total number of oscillators: 256
M = n_communities;				       % number of lowest level communities (what's that?): 8
	
w = ones(1,N);	% identical natural frequencies
h = 0.05;			% Runge-Kutta method step size
T = npoints;		% number of time steps
ws = 2;			% window size for downsampling synchrony data


X(1,:) = rand(1,N)*2*pi-pi;	% random initial phases for all 256 initial oscillators (X will grow in rows each of which is a time-step; columns 1-8 are oscillators of 1st community, columns 9-16 are oscillators of 2nd community and so on)
Dmean = d0+d1;					% average connections per oscillator
synchrony = zeros(T,M);					% synchrony data

tot_synchrony = zeros(1,M);				% running total for synchrony over window
Sync = zeros(T/ws,M);				% downsampled synchrony data
hby2 = 0.5*h;

sigma_chi = zeros(1, T);

for t=2:T
   %disp(t)
   % Simulate the Kuramoto model
   temp1 = X(t-1,:);			% take phases of oscillators from former time-point
   
   for j=1:1/h						% for j = 1:20
      temp2 = temp1;
      
      for i=1:N
		
         % Numerical simulation using 4th-order Runge-Kutta method
         rk1 = w(i);
         for k=1:N
            if coupling_matrix(i,k) ~= 0
               rk1 = rk1 + coupling_matrix(i,k)*sin(temp2(k)-temp1(i)-a)/Dmean;
            end
         end
         rk2 = w(i);
         for k=1:N
            if coupling_matrix(i,k) ~= 0
               rk2 = rk2 + coupling_matrix(i,k)*sin(temp2(k)-(temp1(i)+hby2*rk1)-a)/Dmean;
            end
         end
         rk3 = w(i);
         for k=1:N
            if coupling_matrix(i,k) ~= 0
               rk3 = rk3 + coupling_matrix(i,k).*sin(temp2(k)-(temp1(i)+hby2*rk2)-a)/Dmean;
            end
         end
         rk4 = w(i);
         for k=1:N
            if coupling_matrix(i,k) ~= 0
               rk4 = rk4 + coupling_matrix(i,k).*sin(temp2(k)-(temp1(i)+h*rk3)-a)/Dmean;
            end
	   end
	   
         temp1(i) = temp1(i) + h*(rk1+2*rk2+2*rk3+rk4)/6;
      end
            
   end
   X(t,:) = temp1;		% add sukcessively rows (time-points) of the 256 oscillators

   % Compute synchrony within communities
   for i = 1:M
      for j = 1:n0
         x1 = X(t,(i-1)*n0+j);
         synchrony(t,i) = synchrony(t,i)+exp(x1*sqrt(-1));		% adding all synchrony values of each oscillator belonging to the same community
      end
   end
   synchrony(t,:) = abs(synchrony(t,:)/n0);					% taking the average
   
   synchrony_temp = mean(synchrony(t,:)');
  
   sigma_chi(1,t) = var(synchrony(t,:));
   
   X(t,:) = mod(X(t,:)+pi,2*pi)-pi;	   % normalise phases

end

X = X';
synchrony = synchrony';

global_chi = mean(sigma_chi);
lambda = mean(synchrony, 2);
global_lambda = mean(lambda);
end
