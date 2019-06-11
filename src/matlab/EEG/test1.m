% Must rin `eeglab` before running this script

lis=[];
count=1;
EEG1=pop_biosig('p3.gdf');
EEG2=pop_biosig('p4.gdf');
c={'FP1','FP2','C3','C4','P7','P8','O1','O2','F7','F8','F3','F4','T7','T8','P3','P4','FP1','FP2','C3','C4','P7','P8','O1','O2','F7','F8','F3','F4','T7','T8','P3','P4'};
i = 1;
file= fopen('timepoints.txt', 'w');
fprintf(file,'All Epochs With Synchony');
while i<1005
    [y,record]=plotall(i,EEG1,EEG2);
    if y~=0
%       lis(count)=y;
       fprintf(file,'\n %d th epoch', y);
       fprintf(file,'\nthere are %d connections in this epoch ', length(record));
       for j= 1:length(record)
        lis=record{j};
        fprintf(file,'\n %s to %s ',c{lis(1)},c{lis(2)});
        
       end
       fprintf(file,'\n---------------------------');
       count=count+1;
    end
   
  pause(0.1);
    close all
    i=i+1;
end
fclose(file);