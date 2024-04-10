function filepath = findfile_PA(varargin) 
directory = varargin{1};
pattern   = varargin{2};

filesAndFolders = dir(directory);     % Returns all the files and folders in the directory

filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory                    

numOfFiles = length(filesInDir);

i=1;


while(i<=numOfFiles)

      filename = filesInDir(i).name;                              % Store the name of the file
       
      if (contains(filename,pattern)&&contains(filename,'PA'))
          
          filepath=[directory filename]
          break;
      end
     
      i = i+1;

end
