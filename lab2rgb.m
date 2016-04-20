function y = lab2rgb(input);
% Convert LAB to RGB values
% RGB is _uint8_,  ranging 0..255 on all channels
% LAB is _double_, ranging 0..100 on L*, -127..128 on a* & b*.
%
% G.Sfikas 15/5/2006
%
y = uint8(applycform(lab2uint8(input), makecform('lab2srgb')));
return;