%folder="/n/scratch3/users/j/jaw34/projects/orion/orion_1_74/rawcsv"
folder="/tmp/mat"

% load Orion 1-40
load("/home/jaw34/projects/orion/mat/matlab-Orion_CRC_allbacthes-20220602.mat")
for i = 1:40
  tablename = sprintf('dataC%02d', i)
  if exist(tablename, 'var')
      writetable(eval(tablename), fullfile(folder,sprintf('%s.csv',tablename)));
  else
      fprintf('%s does not exist\n', tablename);
  end
enqd
clear all
  
% load Orion 41-74
  folder="/tmp/mat";
load("/home/jaw34/projects/orion/mat/matlab-Orion_CRC_Cohort2and3-20230202.mat")
vars = who('dataLSP*');
for i = 1:length(vars)
    tablename = vars{i}
    writetable(eval(tablename), fullfile(folder, sprintf('%s.csv',tablename)));
end
clear all
    
% load CyCIF Immune
folder="/tmp/mat";
load ("/home/jaw34/projects/orion/mat/matlab-CRC_Immune_MC-20210826.mat");
vars = who('dataTNPCRC_*');
for i = 1:length(vars)
	  tablename = vars{i}
	  writetable(eval(tablename), fullfile(folder, sprintf('%s.immune.csv', tablename)));
end
clear all

# load CyCIF tumor
folder = "/tmp/mat";
load("/home/jaw34/projects/orion/mat/matlab-CRCWSI_Tumor_MC-20210817.mat");
vars = who('dataTNPCRC_*');
for i = 1:length(vars)
	  tablename = vars{i}
	  writetable(eval(tablename), fullfile(folder, sprintf('%s.tumor.csv', tablename)));
end
clear all
