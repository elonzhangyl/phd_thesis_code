%%
% *This is an example of how to create a 3D bar chart in MATLAB&#174;* .
% 
% You can open this example in the <https://www.mathworks.com/products/matlab/live-editor.html 
% Live Editor> with MATLAB version 2016a or higher.
%
% Read about the <http://www.mathworks.com/help/matlab/ref/bar3.html |bar3|> function in the MATLAB documentation.
% For more examples, go to <http://www.mathworks.com/discovery/gallery.html MATLAB Plot Gallery>
%
% Copyright 2012-2018 The MathWorks, Inc.

% Load monthly temperature data
load MonthlyTemps temperatures months years

% Create the 3D bar chart
figure
bar3(temperatures)
axis([0 13 0 12 0 80])

% Add title and axis labels
title('Boston Monthly Temperatures 1900-2000')
xlabel('Month')
ylabel('Year')
zlabel('Temperature')

% Change the x and y axis tick labels
set(gca, 'XTickLabel', months)
set(gca, 'YTickLabel', years)
