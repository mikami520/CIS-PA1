clc
clear all

format longG

%%% input file which describes the calibration object
table_cb = importdata('pa1-debug-g-calbody.txt');
data_cb = table_cb.data;
textdata_cb = table_cb.textdata{1};
Nd_cb = str2double(textdata_cb(1));
Na_cb = str2double(textdata_cb(4));
Nc_cb = str2double(textdata_cb(7:8));

data_d = data_cb(1:Nd_cb,:);
data_a = data_cb(Nd_cb+1:Nd_cb+Na_cb,:);
data_c = data_cb(Nd_cb+Na_cb+1:size(data_cb,1),:);

%%% input file which provides the values read by the sensor

table_cr = importdata('pa1-debug-g-calreadings.txt');
data_cr = table_cr.data;
textdata_cr = table_cr.textdata{1};
ND_cr = str2double(textdata_cr(1));
NA_cr = str2double(textdata_cr(4));
NC_cr = str2double(textdata_cr(7:8));
num_frame = str2double(textdata_cr(11));
en = 0;
expectC = [];
for i = 1:num_frame
    start = en + 1;
    en = i * (ND_cr + NA_cr + NC_cr);
    data_D = data_cr(start:start+ND_cr-1,:);
    data_A = data_cr(start+ND_cr:start+ND_cr+NA_cr-1,:);
    data_C = data_cr(start+ND_cr+NA_cr:en,:);
    FD = T2TR(data_d, data_D)
    FA = T2TR(data_a, data_A);
    [~, ~, infD] = InvF(FD);
    [~, ~, Ftrans] = FrameTrans(infD, FA);
    for j = 1:length(data_c)
        [~, transc] = TriDxform(Ftrans, data_c(j, :)');
        expectC = [expectC; transc'];
    end
end

output = importdata('pa1-debug-g-output1.txt');
outData = output.data(3:size(output.data,1), :);
disp(abs((outData - expectC) - 0 ) < outData * 0.02)
%%% input file which provides the values read by the sensor
table_op = importdata('pa1-debug-g-optpivot.txt');
data_op = table_op.data;
textdata_op = table_op.textdata{1};

Nd_op = str2double(textdata_op(1));
Nh_op = str2double(textdata_op(4));
N_frame = str2double(textdata_op(7:8));

% Fd
new_H = cell(N_frame,1);
for i = 1:N_frame
    st = (Nd_op + Nh_op) * (i-1) + 1;
    en = st + Nd_op - 1;
    D = data_op(st:en, :);
    H = data_op(en+1:(Nd_op + Nh_op) * i, :);
    Ht = [];
    [~,~,invFd] = InvF(FD);
    for j = 1:Nh_op
        H_ele = H(j, :);
        [~, Ht_ele] = TriDxform(invFd, H_ele');
        Ht = [Ht; Ht_ele'];
    end
    new_H{i} = Ht;
end

H1 = new_H{1};
avgH1 = mean(H1,1);
h1 = H1 - avgH1;


Fmatrix = cell(N_frame,1);
for i = 1:N_frame
    HH = new_H{i};
    Fh = T2TR(h1, HH);
    Fmatrix{i} = Fh;
end

[~, Ppivot] = pivotCalibration(Fmatrix)


