min_err = 1000;
opt_dt = 1;

for dt = 4.365:0.001:4.375
    err = (1*dt/tacc - 1.04)^2;
    err = (2*dt/tacc - 2.09)^2 + err;
    err = (3*dt/tacc - 3.13)^2 + err;
    if (err < min_err)
        min_err = err;
        opt_dt = dt;
    end
end

disp(opt_dt)