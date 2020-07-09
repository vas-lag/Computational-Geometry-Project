clear all;
clc;
% files = ["armadillo_low_low", "b66_L2", "bone", "bunny_low", "cube", "dolphin", "dragon_low_low",...
%     "flashlight", "flashlightNoCentered", "hand2", "icosahedron", "phone_v02", "polyhedron",...
%     "suzanne", "teapotMultiMesh", "unicorn_low", "unicorn_low_low", "vvrlab"];
files = ["armadillo_low_low", "b66_L2", "bone", "bunny_low", "dolphin", "dragon_low_low",...
    "hand2", "phone_v02", "suzanne", "unicorn_low", "unicorn_low_low", "vvrlab"];

count = size(files, 2);
vCount = zeros(1, count);
eigenValues = cell(1, count);
average = zeros(1, count);
mid = zeros(1, count);
spacing = 0.4;
cmap = jet(12);
for i = 1:count
   [vCount(i), eigenValues{1, i}, average(i), mid(i)] = ReadFile(files(i));
   %figure
   divisions = zeros(1, ceil(max(eigenValues{i}) / spacing));
   for e = 1:vCount(i)
       index = ceil(eigenValues{i}(e) / spacing);
       if index == 0
           index = 1;
       end
       divisions(1, index) = divisions(1, index) + 1;
   end
   maxim = 0;
%    for e = 1:ceil(max(eigenValues{i}) / spacing)
%        if (divisions(e) > maxim)
%            maxim = divisions(e);
%        end
%    end
   for e = 1:ceil(max(eigenValues{i}) / spacing)
       maxim = maxim + divisions(e);
   end 
   for e = 1:ceil(max(eigenValues{i}) / spacing)
       divisions(e) = divisions(e) / maxim;
   end
   maxnum = spacing / 2 + (ceil(max(eigenValues{i}) / spacing) - 1) * spacing;
   x = (spacing / 2):spacing:maxnum;
   %subplot(4, 3, i);
   plot(x, divisions, 'Color', cmap(i,:));
   title(files(i));
   hold on;
end 