
function [phase, chi, synchrony] = sim_km_oscillators(npoints, beta, intra_comm_size, n_communities, kuramoto_coupling_matrix)
% Produces metastable chimera-like states in a modular network of
% oscillators. The model is inspired by Abrams, et al., PRL 2008, and is
% parameterised by b and A, according to their paper.
% For more details of this system, see:
%
% Shanahan, M. (2010). Metastable Chimera States in Community-Structured
% Oscillator Networks. Chaos 20, 013108.
%
% Written by Murray Shanahan, August 2009 / December 2014

b = beta;				                        % was 0.1
a = pi/2-b;				                        % phase lag (alpha): 1.4708

n0 = intra_comm_size;						% intra-community size
n1 = n_communities;						% number of communities
d0 = intra_comm_size; 
d1 = intra_comm_size;						% numbers of connections at different community levels
		
N = n_communities*intra_comm_size;				% total number of oscillators: 256
M = n_communities;						% number of lowest level communities (what's that?): 8
	
w = ones(1,N);	                                    % identical natural frequencies
h = 0.05;			                              % Runge-Kutta minterethod step size
T = npoints;		                              % number of time steps
ws = 2;			                              % window size for downsampling synchrony data

% rng(1);
phase(1,:) = rand(1,N)*2*pi-pi;				% random initial phases for all 256 initial oscillators (phase will grow in rows each of which is a time-step; columns 1-8 are oscillators of 1st community, columns 9-16 are oscillators of 2nd community and so on)
Dmean = d0+d1;							% average connections per oscillator
synchrony = zeros(T,M);						% synchrony data
hby2 = 0.5*h;

chi = zeros(1, T);

for t = 2:T
   % Simulate the Kuramoto model
   temp1 = phase(t-1,:);					% take phases of oscillators from former time-point
   
   % do the Runge-Kutta method 1/h (i.e., 20) times, each time using the previously estimated phases as a starting point
   for j = 1:1/h						      
      temp2 = temp1;						% take previously estimated phases as a starting point
      
	% calculate the phase for each oscillator at t
      for i = 1:N
		
		% numerical simulation using 4th-order Runge-Kutta method for oscillator i
		rk1 = w(i);
		for k=1:N
			if kuramoto_coupling_matrix(i,k) ~= 0
				rk1 = rk1 + kuramoto_coupling_matrix(i,k)*sin(temp2(k)-temp1(i)-a)/Dmean;
			end
		end
		rk2 = w(i);
		for k=1:N
			if kuramoto_coupling_matrix(i,k) ~= 0
				rk2 = rk2 + kuramoto_coupling_matrix(i,k)*sin(temp2(k)-(temp1(i)+hby2*rk1)-a)/Dmean;
			end
		end
		rk3 = w(i);
		for k=1:N
			if kuramoto_coupling_matrix(i,k) ~= 0
				rk3 = rk3 + kuramoto_coupling_matrix(i,k).*sin(temp2(k)-(temp1(i)+hby2*rk2)-a)/Dmean;
			end
		end
		rk4 = w(i);
		for k=1:N
			if kuramoto_coupling_matrix(i,k) ~= 0
				rk4 = rk4 + kuramoto_coupling_matrix(i,k).*sin(temp2(k)-(temp1(i)+h*rk3)-a)/Dmean;
			end
		end
	   
		temp1(i) = temp1(i) + h*(rk1+2*rk2+2*rk3+rk4)/6;	% temp1 contains the phase of each oscillator at t
      end
            
   end
   
   phase(t,:) = temp1;							% add sukcessively rows (time-points) of the N oscillators

   % compute synchrony within communities
   for i = 1:M
      for j = 1:n0
         x1 = phase(t,(i-1)*n0+j);
         synchrony(t,i) = synchrony(t,i)+exp(x1*sqrt(-1));	% add all synchrony values of each oscillator belonging to the same community
      end
   end
   synchrony(t,:) = abs(synchrony(t,:)/n0);			% take the average
   
   %synchrony_temp = mean(synchrony(t,:)');
   
%    % compute pairwise synchrony between communities
%    psi = zeros(M,M);
%    for i = 1:M
%       for j = i:M
%          psi(i,j) = mean(abs((synchrony(1:T,i)+synchrony(1:T,j))/2));
%          psi(j,i) = psi(i,j);
%       end
%    end
  
   chi(1,t) = var(synchrony(t,:)');
   
   phase(t,:) = mod(phase(t,:)+pi,2*pi)-pi;			% normalise phases

end

phase = phase';

% cov_err  = eye(N, N);
% mu = zeros(N,1);

% % adding error
% E = mvnrnd(mu, cov_err, npoints)';	
% phase = phase + E;

synchrony = synchrony';

end
