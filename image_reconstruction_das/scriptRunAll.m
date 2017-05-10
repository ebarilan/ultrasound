% methodsOpt = [1,4,7];
methodsOpt = 7;
for m = methodsOpt
    for i = 2
        for j = 1:2
            script_reconstruct_images_from_das(i,j,m);
            close all;
        end
    end
end