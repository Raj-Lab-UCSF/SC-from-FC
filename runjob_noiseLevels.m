
data_MICA = '/wynton/home/rajlab/fabdelnour/data/fMRI_DK';

%data_MICA = '../../Data/fMRI_MICA'; % Data location

parpool(20);

results_struct_noise = struct;

numSubj = 5;
sparcity = [0.01 0.2 0.5 1];

parfor subj = 1:numSubj
	parfor jj = 1:length(sparcity)
    fprintf('Simulation %d/%d \n',subj,numSubj);
    fprintf('Noise %d/%d \n',jj,length(sparcity));
    [ outw, empiricalSC, x0_out, err, err2] = SC_fMRI_MICA_noise( subj , data_MICA , sparcity(jj) );

    results_struct_noise(subj,jj).outw = outw; 
    results_struct_noise(subj,jj).x0_out = x0_out; 
    results_struct_noise(subj,jj).err = err; 
    results_struct_noise(subj,jj).err2 = err2; 
    
	end
	results_struct_noise(subj).empiricalSC = empiricalSC; 
end

delete( gcp('nocreate') );

save('results_struct_noise.mat' , 'results_struct_noise' );