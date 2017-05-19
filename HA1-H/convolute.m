function [val] = convolute(segment, filter);

val = 0;
for k = -1:1
    for l = -1:1
        val = val + segment(2 - l, 2 - k)*filter(2 - l, 2 - k);
        res = segment(2 - l, 2 - k)*filter(2 - l, 2 - k);
        disp([num2str(segment(2 - l, 2 - k)), ' * ', num2str(filter(2 - l, 2 - k)), ' Res: ', num2str(res), ' Val: ', num2str(val)]);
    end
end