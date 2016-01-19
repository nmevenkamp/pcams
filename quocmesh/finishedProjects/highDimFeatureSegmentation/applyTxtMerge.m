function applyTxtMerge ( imgDir, segDir )
    segDirTxtMerge=[segDir '/TxtMerge']
    if ~exist(segDirTxtMerge)
        mkdir(segDirTxtMerge)
    end

    ws=60; % window size

    f1=fspecial('log',[5,5],.8);
    f2=fspecial('log',[7,7],1.2);
    f3=fspecial('log',[9,9],1.8);
    f4=gabor_fn(3.5,pi/2);
    f5=gabor_fn(3.5,0);
    f6=gabor_fn(3.5,pi/4);
    f7=gabor_fn(3.5,-pi/4);
    f8=gabor_fn(2.5,pi/2);
    f9=gabor_fn(2.5,0);
    f10=gabor_fn(2.5,pi/4);
    f11=gabor_fn(2.5,-pi/4);

    N=80;
    for n=1:N
        n
        i1 = floor((n - 1) / 4) + 1;
        i2 = floor(mod(n - 1, 4) / 2) + 1;
        i3 = mod(mod(n - 1, 4), 2) + 1;
        if exist([segDir sprintf('/seg%d_%d_%d.png', i1, i2, i3)],'file')
            img=imread([imgDir sprintf('/tm%d_%d_%d.png', i1, i2, i3)]);

            cf=makecform('srgb2lab');
            Ilab=applycform(img,cf);
            Ig1=subImg(Ilab(:,:,1),f4,f5,f6,f7,f8,f9,f10,f11);
            Ig=cat(3,single(Ilab),Ig1);   

            res0=double(imread([segDir sprintf('/seg%d_%d_%d.png', i1, i2, i3)]));
            res0(:)=res0(:)+1;

            Ig=single(Ig);
            [N1,N2,bn]=size(Ig);
            sh_mx=SHcomp(bn,floor(ws/2),Ig);
            tmp=SHedge_1s(floor(ws/2),1,sh_mx);
            EdgeMap=tmp./max(tmp(:));

            [res_cc]=TxtMerge(res0, EdgeMap);

            imwrite(uint8(res_cc),[segDirTxtMerge sprintf('/seg%d_%d_%d.png', i1, i2, i3)])
        end
    end
end