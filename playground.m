x = [0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1];
y = [0,0.7071,1,0.7071,0,-0.7071,-1,-0.7071,0];
k = 5; % quintic spline as you desire
sp = spapi( optknt(x,k), x, y );
xx = 0:.005:1; % desired sampling
yy = fnval(xx,sp);

plot (x,y,'*',xx,yy)
figure
%,diff(yy,3),diff(yy,4)) %
fnplt(fnder(fnder(fnder(sp))))
plot(diff(yy,3)); % third difference
figure
plot(diff(yy,4)); % fourth difference