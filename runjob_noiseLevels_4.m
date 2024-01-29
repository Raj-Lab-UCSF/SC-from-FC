
data_MICA = '/wynton/home/rajlab/fabdelnour/data/fMRI_DK';

%data_MICA = '../../Data/fMRI_MICA'; % Data location

numSubj = 50;
parpool(numSubj);

results_struct_noise = struct;

sparcity = [0.01 0.2 0.5 1];

jj=4;
parfor subj = 1:numSubj
	%for jj = 1:length(sparcity)
    fprintf('Simulation %d/%d \n',subj,numSubj);
    fprintf('Noise %d/%d \n',jj,length(sparcity));
    [ outw, empiricalSC, noisySC, x0_out, err, err2] = SC_fMRI_MICA_noise( subj , data_MICA , sparcity(jj) );

    results_struct_noise(subj).outw = outw; 
    results_struct_noise(subj).x0_out = x0_out; 
    results_struct_noise(subj).err = err; 
    results_struct_noise(subj).err2 = err2; 
    results_struct_noise(subj).empiricalSC = empiricalSC;
    results_struct_noise(subj).noisySC = noisySC;
	%end
	 
end

delete( gcp('nocreate') );

save('results_struct_noise_4.mat' , 'results_struct_noise' );