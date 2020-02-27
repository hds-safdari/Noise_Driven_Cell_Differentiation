function writelog (type, message)

file = fopen(['logfile.txt'],type);

date=fix(clock);
fprintf(file,'%04d-%02d-%02d  %02d:%02d:%02d\t%s\n',date(1),date(2),date(3),date(4),date(5),date(6), message);

fclose(file);