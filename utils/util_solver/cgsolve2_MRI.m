% cgsolve.m
%
% Solve a symmetric positive definite system Ax = b via conjugate gradients.
%
% Usage: [x, res, iter] = cgsolve(x0, A, b, tol, maxiter, verbose)
%
% x0 - Initial solution
%
% A - Either an NxN matrix, or a function handle.
%
% b - N vector
%
% tol - Desired precision.  Algorithm terminates when 
%    norm(Ax-b)/norm(b) < tol .
%
% maxiter - Maximum number of iterations.
%
% x_orig - original signal (for selecting the best solution)
% 
% verbose - If 0, do not print out progress messages.
%    If and integer greater than 0, print out progress every 'verbose' iters.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function [x, bestx, res, iter, residuals, rms_table] = cgsolve2_MRI(x0, A, b, tol, maxiter, x_orig, verbose)
global ROW COL T_frames STATIC_MASK

static_index = find(STATIC_MASK);
dynamic_index = ~STATIC_MASK;
N_S = length(static_index);
N_D = ROW*COL-N_S;

residuals = []; 
rms_table = [];

if (nargin < 6), verbose = 1; end

implicit = isa(A,'function_handle');

x = x0;
if (implicit), 
    r = b - A(x0);  
else
    r = b - A*x0;  
end
d = r;
delta = r'*r;
delta0 = b'*b;
numiter = 0;
bestx = x;
bestres = sqrt(delta/delta0); 
bestrms = inf;
while ((numiter < maxiter) && (delta > tol^2*delta0))

  % q = A*d
  if (implicit), 
      q = A(d);  
  else
      q = A*d;  
  end
 
  alpha = delta/(d'*q);
  x = x + alpha*d;
   
  if (mod(numiter+1,50) == 0)
    % r = b - Aux*x
    if (implicit), r = b - A(x);  else  r = b - A*x;  end
  else
    r = r - alpha*q;
  end
  
  deltaold = delta;
  delta = r'*r;
  beta = delta/deltaold;
  d = r + beta*d;
  numiter = numiter + 1;
  res = sqrt(delta/delta0);
  if (res < bestres)
    % bestx = x;
    bestres = res;
  end    
  
  rms = sqrt(norm(x-x_orig)^2/length(x));
  % In the recorded experiments, I used wrong definition for rms
  % as sqrt(norm(x-x_orig)/length(x));
  
  residuals = [residuals; res];
  rms_table = [rms_table; rms];
  
  if rms < bestrms && numiter > 5
      bestx = x;
      bestrms = rms;        
  end 
  
  % figure(151);
  % subplot(121); plot(residuals);
  % subplot(122); plot(rms_table);
  
  if ((verbose) && (mod(numiter,verbose)==0))
      disp(sprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e', ...
          numiter, sqrt(sum(bestres.^2)), sqrt(sum((delta/delta0)))));
      x_s = x(1:N_S,:);
      x_d = x(N_S+1:end,:);
      for frame = 1:1%T_frames;
          % f_frame = [fvec_s(N_S_down+1:N_S,:); fvec_d((frame-1)*N_D+1:frame*N_D,:); fvec_s(1:N_S_down,:)];
          f_frame = zeros(ROW, COL);
          f_frame(static_index) = x_s;
          f_frame(dynamic_index) = x_d((frame-1)*N_D+1:frame*N_D);
          figure(1001);
          % subplot(121); 
          imagesc(abs(f_frame)); axis image; drawnow;
          %subplot(122); imagesc(abs(fftshift(fft2(f_frame)))); axis image;
      end
  end
  
end

if (verbose)
  disp(sprintf('cg: Iterations = %d, best residual = %14.8e', numiter, bestres));
end
% x = bestx;
res = bestres;
iter = numiter;

