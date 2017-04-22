%% fdtd_tm
mstep=100;
n = 1;
file_id = fopen('data_tm/e_0001.raw');
mydata{1} = fread(file_id, [60 60], 'double');
F1(n) = getframe(gcf);
fclose(file_id);

for k = 5:5:1800
    n = n+1;
    filename = sprintf('data_tm/e_%4.4d.raw', k);
    file_id = fopen(filename);
    mydata{n} = fread(file_id, [60 60], 'double');
    fclose(file_id);
    mesh(mydata{n});
    zlim([-1 1])
    drawnow
    F1(n) = getframe(gcf);
end

fig = figure;
set(fig,'doublebuffer','on')
movie(fig,F1,1)