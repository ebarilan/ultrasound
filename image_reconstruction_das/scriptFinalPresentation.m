ImageOpt = [1,2];
ProccesOpt = [4, 6]; %[1, 4, 6];
numAngels = 75;
IsEmbedded = [0,1];
acquisition_type = 1;
for ii = IsEmbedded
    for i = ImageOpt
        for j = ProccesOpt
            phantom_type = i;
            data_type = j;
            tic;
            script_reconstruct_images_from_das(acquisition_type , phantom_type , data_type, ii, numAngels);
            fprintf('data = %d. Embedded = %d. Time = %.2f [sec]\n',data_type, ii, toc/60)
        end
    end
end