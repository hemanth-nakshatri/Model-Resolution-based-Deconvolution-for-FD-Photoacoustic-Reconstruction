function [norm] = normalize_img(img)
    max1 = max(max(img));
    min1 = min(min(img));
    norm = (img - min1) / (max1 - min1);
end

