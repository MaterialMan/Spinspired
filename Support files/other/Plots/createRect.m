function coord = createRect(width,height)

coord  = zeros(4,2);
coord(2,2) = width;
coord(3,1) = height;
coord(4,:) = [height width ];
end