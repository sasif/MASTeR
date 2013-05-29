function [u v] = handle_OF(im1, im2, method, opts) 

switch method
    case 'Sun'
        method = opts.method;
        para = opts.para;
        uv = estimate_flow_interface(im1, im2, method,para);
        u = uv(:,:,1);
        v = uv(:,:,2);
    case 'Liu'
        para = opts.para;
        [u v]=Coarse2FineTwoFrames(im1,im2,para);
end