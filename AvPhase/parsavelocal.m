function []=parsavelocal(outFile,dat)
%Support file save in parfor loops
%PARSAVELOCAL(F,D)
%Stores variable D as 'dat' in the file named F. Can be called from within PARFOR
%loops which do not allow direct use of SAVE.

%Soumya D. Mohanty, Mar 2013

save(outFile,'-struct','dat');