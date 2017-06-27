
fig = figure;
files = dir('data_3d_raw');
num = size(files);
n = 0;

for i = 1:num(1)
    if(files(i).isdir == 1)
        continue
    end
    n = n+1;
    fname = strcat(files(i).folder,"\",files(i).name);
    id = fopen(fname, 'r');
    mydata{i} = reshape(fread(id,'double'),[71,71]);
    surf(mydata{i});
    zlim([-1 1]);
    drawnow
    F1(n) = getframe(gcf);
    fclose(id);
end

%%
fig = figure;
set(fig,'doublebuffer','on')
movie(fig,F1,1)
