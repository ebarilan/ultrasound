


ImageOpt = 1; %[1,2]
ProccesOpt = 6; %[1, 4, 6];
numAngels = 11;
IsEmbedded = 1;
acquisition_type = 1; %[1,2] - simulation, experiment


phantom_type = ImageOpt;
data_type = ProccesOpt;
script_reconstruct_images_from_das(acquisition_type , phantom_type , data_type, IsEmbedded, numAngels);
% fprintf('data = %d. Embedded = %d. Time = %.2f [sec]\n',data_type, ii, toc/60)

