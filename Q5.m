clc
clear all

%%% input file which provides the values read by the sensor
table_ep = importdata('pa1-unknown-h-empivot.txt');
data_ep = table_ep.data;
textdata_ep = table_ep.textdata;
Ng_ep = str2double(textdata_ep(1));
N_frame = str2double(textdata_ep(2));

G1 = data_ep(1:Ng_ep,:);
avgG1 = mean(G1,1);
g1 = G1 - avgG1;

en = 0;
Fmatrix = cell(N_frame,1);
for i = 1:N_frame
    start = en + 1;
    en = Ng_ep * i;
    G = data_ep(start:en, :);
    F = T2TR(g1, G);
    Fmatrix{i} = F;
end

[~, Ppivot] = pivotCalibration(Fmatrix)

