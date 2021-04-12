function coord = createRect(width,height)

coord  = zeros(4,2);
coord(2,2) = width;
coord(3,:) = [height width];
coord(4,1) = height;

end