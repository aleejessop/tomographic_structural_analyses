function threshold=NodalSurfaceThresholdForGivenVolumeFraction(volFrac,surfaceName)
    steps=[2  0.5 0.25 0.125 0.0625 0.03 0.015 0.0075 0.003 0.0015 0.001];
    nvoxs=[30 40  50   60    75     90   100   120    150   170    200  ];
    lowerbound=-10;
    upperbound=10;
    step=steps(1);
    nvox=nvoxs(1);
    for i=1:length(steps)
        step=steps(i);
        nvox=nvoxs(i);
        t0=lowerbound;
        binary=createNodalSurface(surfaceName,[nvox,nvox,nvox], 1/nvox, 1, [1 0 0], [0 1 0], [0 0 0], t0);
        phi0=length(find(binary==0))/nvox^3;
        disp([lowerbound,upperbound,phi0,nvox]);
        for t1=lowerbound+step:step:upperbound
            binary=createNodalSurface(surfaceName,[nvox,nvox,nvox], 1/nvox, 1, [1 0 0], [0 1 0], [0 0 0], t1);
            phi1=length(find(binary==0))/nvox^3;
            disp([t1,lowerbound,upperbound,phi1,nvox]);
            if phi1 > volFrac
                lowerbound=t0;
                upperbound=t1;
                break;
            end
            t0=t1;
        end
    end
    threshold=(lowerbound+upperbound)/2;
    threshold;
end

