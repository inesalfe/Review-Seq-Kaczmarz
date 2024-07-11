%% Setup

clc;
clear;

addpath('../../../data/ct/AIRToolsII/')
addpath('../../../data/ct/')
AIRToolsII_setup('temporary')

%% Data

N_pixels = 225;
M = 100703;
N = 50625;
seed = 1;

%% Import Kaczmarz Solution

max_it_CK = 1891001;
max_it_CK_box_proj = 511001;
max_it_RK = 745001;
max_it_RK_box_proj = 591001;
max_it_SRKWOR = 579001;
max_it_SRKWOR_box_proj = 680001;

filename = "../../outputs/ct_real/CK_sol_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_CK) + "_" + int2str(seed) + ".txt";
x_sol_CK = load(filename);

filename = "../../outputs/ct_real/CK_box_proj_sol_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_CK_box_proj) + "_" + int2str(seed) + ".txt";
x_sol_CK_box_proj = load(filename);

filename = "../../outputs/ct_real/RK_sol_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_RK) + "_" + int2str(seed) + ".txt";
x_sol_RK = load(filename);

filename = "../../outputs/ct_real/RK_box_proj_sol_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_RK_box_proj) + "_" + int2str(seed) + ".txt";
x_sol_RK_box_proj = load(filename);

filename = "../../outputs/ct_real/SRKWOR_sol_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_SRKWOR) + "_" + int2str(seed) + ".txt";
x_sol_SRKWOR = load(filename);

filename = "../../outputs/ct_real/SRKWOR_box_proj_sol_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_SRKWOR_box_proj) + "_" + int2str(seed) + ".txt";
x_sol_SRKWOR_box_proj = load(filename);

% Figures

figure1 = figure(1);
imagesc(reshape(x_sol_CK,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_val max_val])
filename_fig = "png/CK_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_CK) + int2str(seed) + ".png";
saveas(figure1,filename_fig);

figure2 = figure(2);
imagesc(reshape(x_sol_CK_box_proj,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_val max_val])
filename_fig = "png/CK_box_proj_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_CK_box_proj) + int2str(seed) + ".png";
saveas(figure2,filename_fig);

figure3 = figure(3);
imagesc(reshape(x_sol_RK,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_val max_val])
filename_fig = "png/RK_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_RK) + int2str(seed) + ".png";
saveas(figure3,filename_fig);

figure4 = figure(4);
imagesc(reshape(x_sol_RK_box_proj,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_val max_val])
filename_fig = "png/RK_box_proj_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_RK_box_proj) + int2str(seed) + ".png";
saveas(figure4,filename_fig);

figure5 = figure(5);
imagesc(reshape(x_sol_SRKWOR,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_val max_val])
filename_fig = "png/SRKWOR_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_SRKWOR) + int2str(seed) + ".png";
saveas(figure5,filename_fig);

figure6 = figure(6);
imagesc(reshape(x_sol_SRKWOR_box_proj,N_pixels,N_pixels))
colorbar
% caxis manual
% caxis([min_val max_val])
filename_fig = "png/SRKWOR_box_proj_" + int2str(M) + "_" + int2str(N) + "_" + int2str(max_it_SRKWOR_box_proj) + int2str(seed) + ".png";
saveas(figure6,filename_fig);