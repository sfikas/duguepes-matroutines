function y = rgb2lab(input)
% Convert RGB to LAB values
% RGB is _uint8_,  ranging 0..255 on all channels
% LAB is _double_, ranging 0..100 on L*, -127..128 on a* & b*.
%
% G.Sfikas 15/5/2006
%
y = lab2double(applycform(uint8(input), makecform('srgb2lab')));
return;