function [wh_x_nu,wh_x_re,x_ce,alphah]=KLIEP22(x_de,x_nu,x_re,sigma_chosen,b,fold)

% x_re=x_nu;

% Kullback-Leiblar importance estimation procedure (with cross validation)
%
% Estimating ratio of probability densities
%   \frac{ p_{nu}(x) }{ p_{de}(x) }
% from samples
%    { xde_i | xde_i\in R^{d} }_{i=1}^{n_{de}} 
% drawn independently from p_{de}(x) and samples
%    { xnu_i | xnu_i\in R^{d} }_{i=1}^{n_{nu}} 
% drawn independently from p_{nu}(x).
%
% Usage:
%       [wh_x_de,wh_x_re]=KLIEP(x_de,x_nu,x_re,sigma_chosen,b)
%
% Input:
%    x_de:         d by n_de sample matrix corresponding to `denominator' (iid from density p_de)
%    x_nu:         d by n_nu sample matrix corresponding to `numerator'   (iid from density p_nu)
%    x_re:         (OPTIONAL) d by n_re reference input matrix
%    sigma_chosen: (OPTIONAL) positive scalar representing Gaussian kernel width;
%                  if omitted (or zero), it is automatically chosen by cross validation
%    b:            (OPTINLAL) positive integer representing the number of kernels (default: 100);
%
% Output:
%    wh_x_de:      estimates of density ratio w=p_nu/p_de at x_de
%    wh_x_re:      estimates of density ratio w=p_nu/p_de at x_re (if x_re is provided)
%
% (c) Masashi Sugiyama, Department of Compter Science, Tokyo Institute of Technology, Japan.
%     sugi@cs.titech.ac.jp,     http://sugiyama-www.cs.titech.ac.jp/~sugi/software/KLIEP/
  
  if nargin<2
    error('number of input arguments is not enough!!!')
  end

  [d,n_de]=size(x_de);
  [d_nu,n_nu]=size(x_nu);
  if d~=d_nu
    error('x_de and x_nu must have the same dimension!!!')
  end
  
  if nargin<4 || isempty(sigma_chosen)
    sigma_chosen=0;
  elseif sigma_chosen<0
    error('Gaussian width must be positive')
  end

  if nargin<5 || isempty(b)
    b = 100;
  end  

  if nargin<6 || isempty(fold)
    fold=5;
  end  

  disp('Run KLIEP')
  %%%%%%%%%%%%%%%% Choosing Gaussian kernel center `x_ce'
  rand_index=randperm(n_nu);
  b=min(b,n_nu);
  x_ce=x_nu(:,rand_index(1:b)); 
  
%   R = mvnrnd(mean(x_nu),cov(x_nu),b)

%   if sigma_chosen==0
%     %%%%%%%%%%%%%%%% Searching Gaussian kernel width `sigma_chosen'
%     sigma=1000;
%     
%     epsilon_list = log10(sigma)-1:-1:0;
%     score = zeros(length(epsilon_list),9,fold,1);
%     wd_opt = zeros(length(epsilon_list),9);
%     for j = 1:length(epsilon_list)
%         epsilon = epsilon_list(j);
%       for k = 1:9 %iteration
%         sigma = sigma-10^epsilon;
%         wd_opt(j,k) = sigma;
%         
%         cv_index=randperm(n_nu);
%         cv_split=floor([0:n_nu-1]*fold./n_nu)+1;
%         
%         
%         X_de=kernel_Gaussian(x_de,x_ce,sigma);
%         X_nu=kernel_Gaussian(x_nu,x_ce,sigma);
%         mean_X_de=mean(X_de,1)';
%         for i=1:fold
%           alpha_cv=KLIEP_learning(mean_X_de,X_nu(cv_index(cv_split~=i),:));
%           wh_cv=X_nu(cv_index(cv_split==i),:)*alpha_cv;
%           score(j,k,i)=mean(log(wh_cv));
%         end
%         mean(score(j,k),3)
%       end %iteration
%     end %epsilon
%     score_mean = mean(score,3);
%     sigma_chosen = wd_opt(score_mean == max(score_mean(:)));
%     fprintf('sigma = %g\n',sigma_chosen)
%   end    
  
  %%%%%%%%%%%%%%%% Computing the final solution `wh_x_de'
  X_de=kernel_Gaussian(x_de,x_ce,sigma_chosen);
  X_nu=kernel_Gaussian(x_nu,x_ce,sigma_chosen);
  mean_X_nu=mean(X_nu,1)';
  alphah=KLIEP_learning(mean_X_nu,X_de);
  wh_x_nu=(X_nu*alphah)';

  if nargin<3 || isempty(x_re)
    wh_x_re=nan;
  else
    X_Re=kernel_Gaussian(x_re,x_ce,sigma_chosen);
    wh_x_re=(X_Re*alphah)';
  end
