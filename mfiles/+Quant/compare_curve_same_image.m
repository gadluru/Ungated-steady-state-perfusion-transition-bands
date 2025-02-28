function value = compare_curve_same_image(Image,n)

if nargin == 1
    n = 1;
end

Image = (squeeze(Image));

temp = figure;
set(temp,'Position',[10,10,1000,1000]);
imagesc(abs(sum(Image,3)));
colormap gray;
brighten(0.6);
axis image;
axis off;
impixelinfo;

for i=1:n
    mask = roipoly();
    nop = sum(sum(mask));
    value(:,i) = squeeze(sum(sum(mask.*abs(Image),1),2))/nop;
end


figure
hold on
for i=1:n
    plot(value(:,i));
    string(i,:) = ['curve',num2str(i)];
end
legend(string);
