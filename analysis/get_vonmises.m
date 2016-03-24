function y = get_vonmises(x,xdata)
% function y = get_vonmises(x,xdata)
% y = x(2) + x(3)*exp(x(4).*(cos((xdata-x(1))/180*pi)-1))+ x(5)*exp(x(4).*(cos((xdata-x(1))/180*pi+pi)-1));

    y = x(2) + x(3)*exp(x(4).*(cos((xdata-x(1))/180*pi)-1))+ x(5)*exp(x(4).*(cos((xdata-x(1))/180*pi+pi)-1));
end
