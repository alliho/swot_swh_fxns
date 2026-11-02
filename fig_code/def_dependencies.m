%% define stats

getnn = @(x,y) find(~isnan(x) & ~isnan(y));
getnnn = @(x,y,z) find(~isnan(x) & ~isnan(y) & ~isnan(z));
biasfxn = @(x,y) nanmean((y-x));
rmsefxn = @(x,y) sqrt(nanmean((x-y).^2));
crmsefxn = @(x,y) sqrt(rmsefxn(x,y).^2 - biasfxn(x,y).^2);
normrmsefxn = @(x,y) sqrt(nanmean((x-y).^2))./nanmean(abs(x));
absrmsefxn = @(x,y) sqrt(nanmean((abs(x)-abs(y)).^2));
nfxn = @(x,y) sum(~isnan(x) & ~isnan(y));
ccfxn = @(x,y) mode(diag(flipud(corrcoef(x(getnn(x,y)),y(getnn(x,y))))));
sifxn = @(x,y) 1./nanmean(x(getnn(x,y))) .* sqrt( 1./length(x(getnn(x,y))) .* nansum( ( (y(getnn(x,y))-nanmean(y(getnn(x,y)))) - (x(getnn(x,y))-nanmean(x(getnn(x,y))))   ).^2 ) );;
linfxn = @(x,y) polyfit(x(getnn(x,y)), y(getnn(x,y)), 1);
getind = @(var,i) var(i); 
intcptfxn = @(x,y) getind(linfxn(x,y),2); 
slpfxn = @(x,y) getind(linfxn(x,y),1); 
slp0fxn = @(x,y) y(getnn(x,y))/x(getnn(x,y));
