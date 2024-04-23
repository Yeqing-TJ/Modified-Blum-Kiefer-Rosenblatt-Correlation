function UTILITY = modifiedBKRcorr(X,Y)

% INPUT
% X: an n x p matrix
% Y: an n x 1 vector

% OUTPUT: FOUR RESULTS WILL BE REPORTED
% UTILITY.HOEF = int {cov(I(X <= x), I(Y <= y))}^2 dF(x,y)
% UTILITY.BKRC = int {cov(I(X <= x), I(Y <= y))}^2 dF(x)dF(y)
% UTILITY.DIAG = int {corr(I(X <= x), I(Y <= y))}^2 dF(x,y)
% UTILITY.FULL = int {corr(I(X <= x), I(Y <= y))}^2 dF(x)dF(y)

% REFERENCE
% [1] Yeqing Zhou and Liping Zhu (2018) Model-Free Feature Screening 
%         for Ultrahigh Dimensional Data through a Modified
%         Blum-Kiefer-Rosenblatt Correlation. Statistica Sinica, 28(3): 1351-1370.

% [2] Hoeffding, W. (1948). A non-parametric test of independence. 
%         The Annals of Mathematical Statistics, 19(4):546?557.

% [3] Blum, J. R., Kiefer, J., and Rosenblatt, M. (1961). Distribution 
%         free tests of independence based on the sample distribution 
%         function. The Annals of Mathematical Statistics, 32(2):485?498.



% EXAMPLE
% n = 200;
% p = 2000;
% ActiveSet = [1,2,20];
% rho = 0.9;
% Sig = rho.^abs((1:p)'*ones(1,p) - ones(p,1)*(1:p));
% x = mvnrnd(zeros(p,1),Sig,n);
% beta = zeros(p,1);
% beta(ActiveSet) = 1;
% z = x;
% z(:,ActiveSet(end)) = x(:,ActiveSet(end)).^2;
% y = z * beta + normrnd(0,1,n,1);
% 
% tic,UTILITY = modifiedBKRcorr(x,y);toc
% 
% ModelSize = (1:3)*n/log(n)
% [MMS, SelSize] = ScreeningSelection(UTILITY.FULL,ActiveSet,ModelSize)


% COPYRIGHT
% LIPING ZHU @ OCT 25 2015

[n, p] = size(X);                                                          % size of X
[m, d] = size(Y);                                                          % size of Y
if n ~= m || d ~= 1
    error('Y must be univariate and has the same length as X')
end

Aone = ones(n,1); oneA = ones(1,n); 
UTILITY.FULL = zeros(p,1);                                                 % initialization
UTILITY.DIAG = zeros(p,1);                                                 % initialization
UTILITY.HOEF = zeros(p,1);                                                 % initialization
UTILITY.BKRC = zeros(p,1);                                                 % initialization

indY = Y(:)*oneA <= Aone*Y(:)';                                            % I( Y(ii) <= Y(jj) ) 

FY = mean(indY);
ctrY = indY - Aone*FY;


for kk = 1:p
    indX = X(:,kk)*oneA <= Aone * X(:,kk)';                                % I(X(ii,kk) <= X(jj,kk))

    FX = mean(indX);
    ctrX = indX - Aone*FX;

    UTILITY.HOEF(kk) = mean((mean(ctrY.*ctrX)).^2);                        % HOEFFDING correlation using integral of dF(X,Y)
    UTILITY.BKRC(kk) = mean(mean((ctrY'*ctrX/n).^2));                      % BKR correlation using integral of dF(X)dF(Y)  

    Rcor = corr(indX,indY);
    Rdiag = diag(Rcor);
    Rcor(isnan(Rcor) == 1) = [];
    Rdiag(isnan(Rdiag) == 1) = [];
    UTILITY.FULL(kk) = mean(mean(Rcor(:).^2));                             % modified BKR correlation using integral of dF(X)dF(Y)
    UTILITY.DIAG(kk) = mean(mean(Rdiag(:).^2));                            % modified BKR correlation using integral of dF(X,Y)
end





    
% EOF    