ImageOpt = [1,2];
ProccesOpt = 7;%[4, 6]; %[1, 4, 6];
numAngels = 5;
IsEmbedded = 1;
acquisition_type = 1;
for i = ImageOpt
    for ii = IsEmbedded
        for inumAngels = numAngels
            for j = ProccesOpt
                phantom_type = i;
                data_type = j;
                tic;
                script_reconstruct_images_from_das(acquisition_type , phantom_type , data_type, ii, inumAngels);
                fprintf('data = %d. Embedded = %d. Time = %.2f [sec]\n',data_type, ii, toc/60)
            end
        end
    end
end