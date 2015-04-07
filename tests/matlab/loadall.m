function loadall(dirname, savename)

nj50=dataset('File', fullfile(dirname, 'nj50_polytomysolver.csv'));
nj75=dataset('File', fullfile(dirname,'nj75_polytomysolver.csv'));
nj90=dataset('File', fullfile(dirname,'nj90_polytomysolver.csv'));
nj95=dataset('File', fullfile(dirname,'nj95_polytomysolver.csv'));
nj100=dataset('File',fullfile(dirname, 'nj100_polytomysolver.csv'));
nj101=dataset('File', fullfile(dirname,'nj101_polytomysolver.csv'));
njtree=dataset('File', fullfile(dirname,'nj_treefix.csv'));
njlkl=dataset('File', fullfile(dirname,'nj_likelihood.csv'));

nj_alt50=dataset('File', fullfile(dirname, 'nj_alt50_polytomysolver.csv'));
nj_alt75=dataset('File', fullfile(dirname,'nj_alt75_polytomysolver.csv'));
nj_alt90=dataset('File', fullfile(dirname,'nj_alt90_polytomysolver.csv'));
nj_alt95=dataset('File', fullfile(dirname,'nj_alt95_polytomysolver.csv'));
nj_alt100=dataset('File',fullfile(dirname, 'nj_alt100_polytomysolver.csv'));
nj_alttree=dataset('File', fullfile(dirname,'nj_alt_treefix.csv'));
nj_altlkl=dataset('File', fullfile(dirname,'nj_alt_likelihood.csv'));

save(savename);
clear all;

end
