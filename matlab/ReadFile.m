function [vCount, eigenValues, average, mid] = ReadFile(Name)
s = strcat('C:\Users\Billy\Documents\GitHub\Computational-Geometry-Project\resources\eigen\', Name, '.txt');
f1 = fopen(s, 'r' );
vCount = fscanf(f1, '%d', 1);
eigenValues = zeros(1, vCount);
for i = 1:vCount
    eigenValues(1, i) = fscanf( f1, '%f', 1 );
end
fclose(f1);
average = sum(eigenValues) / vCount;
mid = median(eigenValues);