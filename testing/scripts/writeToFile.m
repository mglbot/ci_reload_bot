
% This file is used to generate golden input data for unit testing.
% DO NOT MODIFY

function writeToFile(path, toWrite)
%writes object to file in double prec, binary format readable by matrixIO.hpp
fd = fopen(path,'w');
fwrite(fd,toWrite,'double');
fclose(fd);
end

