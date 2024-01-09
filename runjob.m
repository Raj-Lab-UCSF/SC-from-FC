
data_MICA = '/wynton/home/rajlab/fabdelnour/data/fMRI_DK';

%data_MICA = '../../Data/fMRI_MICA'; % Data location

parpool(25);

results_struct = struct;

numSubj = 50;

parfor subj = 1:numSubj
    fprintf('Simulation %d/%d \n',subj,numSubj);
    [ outw, empiricalSC, x0_out, err, err2] = SC_fMRI_MICA( subj , data_MICA );

    results_struct(subj).outw = outw; 
    results_struct(subj).empiricalSC = empiricalSC; 
    results_struct(subj).x0_out = x0_out; 
    results_struct(subj).err = err; 
    results_struct(subj).err2 = err2; 
end


delete( gcp('nocreate') );

%save('results_struct.mat' , 'results_struct');

save('results_struct.mat' , 'results_struct' );