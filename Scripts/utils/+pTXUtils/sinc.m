function y = sinc(x)

% sinc function

y = sin(x)./x;
y(x==0) = 1;