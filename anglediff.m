function angle = anglediff(vec1r,vec1c,vec2r,vec2c)

dp = vec1r*vec2r + vec1c*vec2c;
mag1 = (vec1r^2 + vec1c^2)^(0.5);
mag2 = (vec2r^2 + vec2c^2)^(0.5);

angle = acos(dp/(mag1*mag2));