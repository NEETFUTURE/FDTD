%% fdtd1d
mstep=1800;
n = 1;
fig = figure;

for k = 10:5:1800
    myfilename = sprintf('e_%4.4d.txt', k);

    mydata{n} = importdata(myfilename);
    plot(mydata{n});
    ylim([-1 1]);
    drawnow
    F1(n) = getframe(gcf);
    n = n+1;
end
%%
fig = figure;
set(fig,'doublebuffer','on')
movie(fig,F1,1)


%% fdtd_tm
mstep=100;
n = 1;
mydata{n} = reshape(importdata('data_tm/e_0001.txt'),[],60);
F1(n) = getframe(gcf);

for k = 5:5:1800
    n = n+1;
    myfilename = sprintf('data_tm/e_%4.4d.txt', k);
    mydata{n} = reshape(importdata(myfilename),[],60);
    mesh(mydata{n});
    zlim([-1 1])
    drawnow
    F1(n) = getframe(gcf);
end

%% ƒtƒŒ[ƒ€‚ğÄ¶
fig = figure;
set(fig,'doublebuffer','on')
movie(fig,F1,1)
%% “®‰æ‚ğ•Û‘¶
a = VideoWriter('movie.avi');
open(a)
writeVideo(a,F1)
close(a)

