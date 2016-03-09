function [y,Wt,J] = convbss(x,T,Q,K,N,verbose,iter,printevery,eta,ds);
%CONVBSS demixes convolutively mixed signals.
%
% This implements an algorithm for blind source separation of
% convolutive mixtures that has been published in 
%
%    Lucas Parra, Clay Spence, 
%    "Convolutive blind source separation of non-stationary sources", 
%    IEEE Trans. on Speech and Audio Processing pp.320-327, May 2000.
%
% and it has been patented (US patent US6167417).
%
% The input arguments:
%   x       - the mixed signals (along the columns)
%   T       - the length of the FFT (default==512)
%   Q       - the length of the filter in time domain (default==128)
%   K       - the number of matrices to diagonalize (default==5)
%   N       - the number of intervals used to estimate each
%             cross-power-matrix (default==7)
%   verbose - if >0 the progress is printed out (default==1)
%   iter    - the maximum number of iterations (default==1000)
%   printevery - how often should be something printed (default==10)
%   eta     - the learning rate (default==1.0)
%   ds      - the number of sources to be estimated
%             (default==number_of_columns_of_x)
%
%Copyright notice:
%
%  author: Stefan Harmeling
%          GMD FIRST Berlin
%              Germany
%
%THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE of GMD FIRST Berlin.
%  The purpose of this software is the dissemination of
%  scientific work for scientific use. The commercial
%  distribution or use of this source code is prohibited, 
%  as is any modification of its content.
%The copyright notice does not evidence any actual or intended 
%publication of this code.
%
%Term and Termination:
%
%This License shall continue for as long as you use the Software. 
%However, it will terminate if you fail to comply with any of its 
%terms and conditions. You agree, upon termination, to discontinue
%using, and to destroy all copies of, the Software. 
%The Limitations of Warranties and Liability set out below shall 
%continue in force even after any termination.
%
%Limitation of Warranties and Liability:
%
%THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR
%IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%HEREBY DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
%DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
%OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
%HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
%STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
%ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%POSSIBILITY OF SUCH DAMAGE.
%
%----------- 2001 GMD FIRST Berlin - Stefan Harmeling -----------
%



% the size is needed very early
[rx,cx] = size(x);

% assign default values to the unassigned
if ~exist('T')|isempty(T),                 T = 512; end
  if mod(T,2), error('T must be multiple of 2'), end
  hT = T/2+1;        % the relevant part of the FFT
if ~exist('Q')|isempty(Q),                 Q = 128; end
if ~exist('K')|isempty(K),                 K = 5; end
if ~exist('N')|isempty(N),                 N = floor(rx/T/K); end
if ~exist('verbose')|isempty(verbose),     verbose = 1; end
if ~exist('iter')|isempty(iter),           iter = 1000; end
if ~exist('printevery')|isempty(printevery), printevery = 10; end
if ~exist('eta')|isempty(eta),             eta = 1.0; end
if ~exist('ds')|isempty(ds),               ds = cx; end

%  we want some of the large variables to be global
global Rx Rxk m dsdshT dsdxhT

% the mixed signals must have zero mean and variance one
x=(x-repmat(mean(x),length(x),1))./std(x(:)); 

%x = (x-repmat(mean(x),[rx 1]))./repmat(std(x),[rx 1]);
dx = cx;  % the number of senors is the number of columns
if ds>dx, error('ds is not allowed to be greater than dx'), end

% calculate the cross-power-spectrum Rx
Rx = cps(x,T,hT,K,N,dx,verbose);  % size(Rx) == [dx dx hT K]

% calculate some straightforward power normalization
Rxk = sum(Rx,4);   % sum over the k matrices
m = 1./repmat(sum(sum(Rxk.*conj(Rxk),2),1),[ds dx 1]);

% initialize the demixing filter to the identity filter
Wt = zeros(ds,dx,T);    % Wt = W in time domain
Wt(:,:,1) = eye(ds,dx);
Wf = shortfft(Wt);      % Wf = W in frequency domain

% define the indices of the diagonal of some matrices
dsdshT = find(repmat(eye(ds,ds),[1 1 hT]));
dsdxhT = find(repmat(eye(ds,dx),[1 1 hT]));

oldJ = inf;

% do the gradient descent
while iter > 0
 if verbose, fprintf('iter=%d         \r',iter), end

 if verbose & (mod(iter,printevery)==0 | iter==1)
   % calculate the error
   J = calcJ(Wf,hT,K);

   % check whether we still made progress
   if oldJ<J
     % we are no longer making progress, let's break this loop 
     break
   else
     oldJ = J;
   end
   fprintf('iter=%d J=%f \n',iter,J)
 end

 % calculate the gradient of the error
 dWf = calcdWf(Wf,ds,dx,hT,K);

 % update the unmixing filter Wf
 Wf = Wf - eta*dWf;

 % project the filter matrix such that Wt(:,:,Q:end) == 0
 Wt = shortifft(Wf);
 Wt(:,:,Q+1:end) = 0;
 Wf = shortfft(Wt);

 iter = iter - 1;
end

% calculate the error
J = calcJ(Wf,hT,K);

% finally calculate the demixed signals
if verbose, fprintf('calculate the demixed signals (this might take a while)...'), end
Wt = shortifft(Wf);
y = demix(Wt,x);
if verbose, fprintf('done\n'), end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = demix(Wt,x)
% this is very slow, better would be something 
% like overlap-save, see Haykin or Oppenheim/Schaefer
% it is assumed that size(Wt) = [ds dx T]
%                    size(x)  = [?? ds]

[ds dx T] = size(Wt);
[rx cx]   = size(x);

% check for consistency
if ds ~= cx
  error('Wt and x must fit together')
end

% allocate memory for y
y = zeros(rx,ds);

% do the convolution using matlab's built-in filter
for ii = 1:ds
 for jj = 1:dx
   y(:,ii) = y(:,ii) + filter(Wt(ii,jj,:),1,x(:,jj));
 end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rx = cps(x,T,hT,K,N,dx,verbose);
% calculate the cross-power-spectrum Rx

% check whether these parameters can be used with signals of length rx,
% this check is done by checking whether the last index that will be
% accessed during the calculation of Rx exists
if K*T*N > size(x,1)
 error('the signals are too short for the chosen parameters')
end

% do the calculation
if verbose, disp('estimate the cross-power spectrum of x'), end
Rx = zeros(dx,dx,hT,K);
XX = zeros(dx,dx,hT);
for k=1:K
 for n=1:N
   if verbose, fprintf('%d/%d %d/%d  \r',k,K,n,N); end
   indx = (k-1)*T*N + (n-1)*T + (1:T);
   X=fft(x(indx,:));
   % since we are only interested in filter in frequency domain
   % that correspond to real valued filters in time domain,
   % we can ignore the negative frequencies
   X=X(1:hT,:);
   for row = 1:dx, for col = 1:dx
     XX(row,col,:) = reshape(conj(X(:,col)).*X(:,row),[1 1 hT]);
   end, end
   Rx(:,:,:,k) = Rx(:,:,:,k) + XX;
 end;
end;
Rx = Rx/N/T;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Wf = shortfft(Wt)
Wf = fft(Wt,[],3);
Wf = Wf(:,:,1:end/2+1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Wt = shortifft(Wf)
Wt = real(ifft(cat(3,Wf,conj(Wf(:,:,[end-1:-1:2]))),[],3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = fast_ip(A,B)
% for high-(more than three)-dimensional arrays with a lot
% of elements in those higher dimensions it might be very 
% expensive to calculate the inner product.  this function 
% implements it in a fast way
[rA cA sA] = size(A);
[rB cB sB] = size(B);
if cA~=rB, error('the inner dimensions of A and B do not agree'), end
if sA~=sB, error('the higher dimensions of A and B do not agree'), end
C = zeros(rA,cB,sB);
for row=1:rA, for col=1:cB
 C(row,col,:) = sum(reshape(A(row,:,:),[cA 1 sA]) .* B(:,col,:),1);
end, end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ah = hermite3(A)
% for three dimensional arrays the function ' can not be used
% to conjugate-transpose (hermite) a matrix
Ah = permute(conj(A),[2 1 3]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = calcJ(Wf,hT,K);
global Rx dsdshT
% calculate the error term for diagonalization
J = 0;
for k = 1:K
 Eomega = fast_ip(fast_ip(Wf,Rx(:,:,:,k)),hermite3(Wf));
 Eomega(dsdshT) = 0;  % this uses the optimal lambdaS
 % note that trace(A*A') == A(:)' * A(:)
 J = J + Eomega(:)' * Eomega(:);
end
J = real(J)/K/hT;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dWf = calcdWf(Wf,ds,dx,hT,K);
global Rx m dsdshT dsdxhT

% calculate the gradient of the objective function
dWf = zeros(ds,dx,hT);
hWf = hermite3(Wf);
for k = 1:K
 WfRx = fast_ip(Wf,Rx(:,:,:,k));
 Eomega = fast_ip(WfRx,hWf);
 Eomega(dsdshT) = 0;  % this uses the optimal lambdaS
 dWf = dWf + fast_ip(Eomega,WfRx);
end
dWf = 2*dWf;

%% apply the straightforward power normalization
dWf = m.*dWf;

% enforce unit gains of the diagonal filters
dWf(dsdxhT) = 0;

