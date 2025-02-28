function Image = sr_backward(Image,siz)

[sx,sy,nof,~,~,~] = size(Image);
Image = reshape(Image,[sx*sy*nof,siz.nSMS*siz.nset]);
Image = Image(:,siz.order);
Image = reshape(Image,[sx,sy,nof,siz.nSMS*siz.nset]);

end