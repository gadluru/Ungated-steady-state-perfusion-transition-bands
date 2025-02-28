function show_yt(img)
figure
img = squeeze(img);
nof = size(img,3);
i = 1;
count = 0;
img = brighten(img,0.4);
while i
    j = mod(i,nof);
    if j == 0
        j = nof;
        count = count + 1;
        if count == 10
            break
        end
    end 
    imagesc(abs(img(:,:,j)))
    colormap gray
    axis image
    title(num2str(j))
    % brighten(0.4)
    drawnow
    i = i+1;
    pause(0.1)
end