%  Sample Bullseye Creation

% Create the AHA 17-segment bullseye
figure;
c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
set(c,'Color','k','LineWidth',2)

% Example 1 of filling the bullseye, vector by vector
fillBullseye(rand(4,1),0.5,1,-45,315);
fillBullseye([0.25 0.35 0.25 1 0.25 0.5],1,1.5,0,360);
fillBullseye( ( 1:10:360 )/360,1.5,2,0,360);
uistack(c,'top');

%Example 2 of filling the bullseye by creating a matrix.
figure;
[X,Y] = meshgrid([0:0.05:1],[0:0.05:1]);
mat = X.*Y;

c = createBullseye([0 0.5 1 0; 0.5 1 4 45; 1 1.5 6 0; 1.5 2 6 0]);
fillBullseye(mat);
uistack(c,'top');


