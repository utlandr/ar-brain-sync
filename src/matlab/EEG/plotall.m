function [val,points]= plotall(x,EEG1,EEG2)

val=0;
points={};
count=1;
a=(EEG1.data(:, x));
sec1 = transpose(a);


b=(EEG2.data(:, x));
sec2 = transpose(b);
inp=horzcat(sec1, sec2);
%disp(inp);


i=1;
j=1;
while i <33
    val1=inp(i);
    j=1;
    while j<33
        val2=inp(j);
        new(i,j)=abs(val2-val1);
       
       if ((i>16 && j<17) || (j>16 && i<17)) && (new(i,j)<=0.01) 
           val=x;
          points{count}=[i,j];
          count=count+1;
       end
       j=j+1;
    end
    i= i+1;
end
% disp(new);
% pop_eegplot(EEG1);
% pop_eegplot(EEG2);
% figure;
% word=sprintf('Pilot Study at %s ms',num2str(x));
%dualheadplot('p3.spl','values2plot',inp,'maplimits',[min(inp),max(inp)],'electrodes','on','cbar',0,'labels',2,'title',word);
%To plot the variables (eg. the power) on the scalp, with the conexion links ()eg synchrony, coherence, etc whose values are between 0.99 and 1
 dualheadplot('p3.spl','adjacencyMatrix',new,'thresholdLink',[0 0.01],'values2plot',inp,'maplimits',[min(inp),max(inp)],'electrodes','on','cbar',0,'labels',0);



end