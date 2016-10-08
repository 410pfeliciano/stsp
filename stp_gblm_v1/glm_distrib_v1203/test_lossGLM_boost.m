for i=1:100

% Fisher Scoring Test...

n = 1000;
p = 100;
b = randn(p,1);
b(50:end)=0;
X = randn(n,p)/10;
% loss_fun = @lossGLM_gauss_iden;
% y = X*b+randn(size(X,1),1);

% loss_fun = 'lossGLM_poiss_exp';
% y = poissrnd(exp(X*b));

loss_fun = @lossGLM_poiss_log;
y = poissrnd(log(1+exp(X*b)));

opts=optimset('Gradobj','on','Hessian','on');

b0 = randn(p,1);
bmle = fminunc(loss_fun,b0,opts,X,y,0);

llhd=[];
bfs = b0;
maxiter=100;
for iter=1:maxiter
    [llhd(iter),dx,H] = feval(loss_fun,bfs,X,y,0);
    bfs = bfs - pinv(H)*dx;
end

diffN(i)=norm(bfs-bmle);
end

%%

for i=1:size(path,2)
    llhdt(i) = feval(loss_fun,path(:,i),X,yt);
end