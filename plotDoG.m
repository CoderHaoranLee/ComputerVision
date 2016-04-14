function plotDoG( DoG )

octave_num = length(DoG);
[~,~,s] = size(DoG{1});
figure
hold on
for i = 1:octave_num*s
    subplot(octave_num,s,i)
    s_id = mod(i,s);
    if(s_id == 0)
        s_id = s;
    end
    o=floor((i-1)/s)+1;
    imshow(uint8(255 * mat2gray(DoG{o}(:,:,s_id))))
end

end

