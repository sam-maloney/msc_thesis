dt = (15+16)/16/3;
dt = round(dt,3);

err = (3*dt/tacc - 0.46)^2;
err = err + (9*dt/tacc - 1.39)^2;
err = err + (16*dt/tacc - 2.47)^2;
disp(dt);
disp(err);

% dt = (16+17)/17/3;
% dt = round(dt,3);
% 
% err = (3*dt/tacc - 0.46)^2;
% err = err + (9*dt/tacc - 1.39)^2;
% err = err + (16*dt/tacc - 2.47)^2;
% disp(dt);
% disp(err);