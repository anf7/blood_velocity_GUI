function smask = findstep(dilateimg, iimg, step)

%#codegen

smask = dilateimg & (iimg == step);