
function [B R2 res]=ols(Y,X,constant)
%[B R2 res]=ols(Y,X,constant)
%
% Estimate Linear Regression Model through OLS 
%
% INPUT:
% Y is a column vector N x 1 (Response variable)
% X is matrix N x P where each column is a predictor.
% constant is a boolean value equal to: 1 to estimate the model with the intercept and 0 otherwise 
%
% OUTPUT:
% B is the vector of beta coefficients
% R2 is the goodness of fit
% res is the vector of residuals

n=size(X,1);
if constant==1
    X=[ones(n,1) X];
end

XX=inv(X'*X)*X';
H=X*XX;

B=XX*Y;
Yhat=H*Y;
res=Y-Yhat;
DevRes=sum(res.^2);
DevTot=var(Y,1)*n;
R2=1-(DevRes/DevTot);