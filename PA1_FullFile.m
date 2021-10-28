clc
clear all

k = 'g'; % define which number of file to use, 'a' to 'g' are for debugging, 'h' to 'k' are for testing.
x = 'debug'; % define which type of file to use, 'unknown' or 'debug'.
%% Q4
%%% input file which describes the calibration object
filename_1 = ['pa1-',x,'-',k,'-calbody.txt'];
table_cb = importdata(filename_1);
data_cb = table_cb.data;
textdata_cb = table_cb.textdata{1};
Nd_cb = str2double(textdata_cb(1));
Na_cb = str2double(textdata_cb(4));
Nc_cb = str2double(textdata_cb(7:8));

data_d = data_cb(1:Nd_cb,:);
data_a = data_cb(Nd_cb+1:Nd_cb+Na_cb,:);
data_c = data_cb(Nd_cb+Na_cb+1:size(data_cb,1),:);

%%% input file which provides the values read by the sensor
filename_2 = ['pa1-',x,'-',k,'-calreadings.txt'];
table_cr = importdata(filename_2);
data_cr = table_cr.data;
textdata_cr = table_cr.textdata{1};
ND_cr = str2double(textdata_cr(1));
NA_cr = str2double(textdata_cr(4));
NC_cr = str2double(textdata_cr(7:8));
num_frame = str2double(textdata_cr(11));

%%% Calculation of C^expected
en = 0;
expectC = [];
for i = 1:num_frame
    start = en + 1;
    en = i * (ND_cr + NA_cr + NC_cr);
    data_D = data_cr(start:start+ND_cr-1,:);
    data_A = data_cr(start+ND_cr:start+ND_cr+NA_cr-1,:);
    data_C = data_cr(start+ND_cr+NA_cr:en,:);
    FD = T2TR(data_d, data_D);
    FA = T2TR(data_a, data_A);
    [~, ~, infD] = InvF(FD);
    [~, ~, Ftrans] = FrameTrans(infD, FA);
    for j = 1:length(data_c)
        [~, transc] = TriDxform(Ftrans, data_c(j, :)');
        expectC = [expectC; transc'];
    end
end

%%% Validation of C^expected
output_cmp = strcmp('debug',x);
if output_cmp == 1
    filename_output = ['pa1-debug-',k,'-output1.txt'];
    output = importdata(filename_output);
    outData = output.data(3:size(output.data,1), :);
    disp("Compare real expect C to example data: ");
    disp(abs((outData - expectC) - 0 ) < 0.02*outData)
else
end

%% Q5
%%% input file which provides the values read by the sensor
filename_3 = ['pa1-',x,'-',k,'-empivot.txt'];
table_ep = importdata(filename_3);
data_ep = table_ep.data;
textdata_ep = table_ep.textdata;
Ng_ep = str2double(textdata_ep(1));
N_frame = str2double(textdata_ep(2));

%%% Using the first "frame" of pivot calibration data to define a local "probe" coordinate system and use this to comput g_j
G1 = data_ep(1:Ng_ep,:);
avgG1 = mean(G1,1);
g1 = G1 - avgG1;

%%% Getting all F_G[k] such that G_j = F_G[k] * g_j
en = 0;
Fmatrix = cell(N_frame,1);
for i = 1:N_frame
    start = en + 1;
    en = Ng_ep * i;
    G = data_ep(start:en, :);
    F = T2TR(g1, G);
    Fmatrix{i} = F;
end

%%% Calculate the P_pivot (P_dimple) using calibration method
[~, Ppivot] = pivotCalibration(Fmatrix)

%% Q6
format longG
%%% input file which provides the values read by the sensor
filename_4 = ['pa1-',x,'-',k,'-optpivot.txt'];
table_op = importdata(filename_4);
data_op = table_op.data;
textdata_op = table_op.textdata{1};

Nd_op = str2double(textdata_op(1));
Nh_op = str2double(textdata_op(4));
N_frame = str2double(textdata_op(7:8));


%%% calculate the H after transformed into EM tracker coordinates.
new_H = cell(N_frame,1);
for i = 1:N_frame
    st = (Nd_op + Nh_op) * (i-1) + 1;
    en = st + Nd_op - 1;
    D = data_op(st:en, :);
    H = data_op(en+1:(Nd_op + Nh_op) * i, :);
    [~,~,invFd] = InvF(FD);
    Ht = [];
    for j = 1:Nh_op
        H_ele = H(j, :);
        [~, Ht_ele] = TriDxform(invFd, H_ele');
        Ht = [Ht; Ht_ele'];
    end
    new_H{i} = Ht;
end

%%% using the same method as Q5 to get h_1
H1 = new_H{1};
avgH1 = mean(H1,1);
h1 = H1 - avgH1;

%%% Getting all F_H[k] such that H_j = F_H[k] * h_j
Fmatrix = cell(N_frame,1);
for i = 1:N_frame
    HH = new_H{i};
    Fh = T2TR(h1, HH);
    Fmatrix{i} = Fh;
end

%%% Calculate the P_pivot (P_dimple) using calibration method in EM tracker
%%% coordinates
[~, Ppivot_EM] = pivotCalibration(Fmatrix)

%%% Export data
FULLDATA = [Ppivot Ppivot_EM expectC'];
file_name = ['pa1-real-',x,'-',k,'-output.txt'];
fileID = fopen(file_name,'wt+');
fprintf(fileID,strcat('27, 8, ', file_name, '\n'));
fprintf(fileID,'%.2f %.2f %.2f\n',FULLDATA);
fclose(fileID);