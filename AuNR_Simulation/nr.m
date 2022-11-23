close all; clc;clearvars;




%%%% Reading exp. value of RI from Johnson & Christy (1972) given at the end of program
experimental_refractive_index = data_john_christy();
wavelength_john = experimental_refractive_index (:,1);
n_exp = experimental_refractive_index (:,2);
k_exp = experimental_refractive_index (:,3);
%%%% reading experimental UV-Vis values from external file
filename = 'data'; %%% UV_VIS spectrum file name (without header):
filename = [filename '.txt']; %%%% UV-Vis file format (wavelength, Absorption (nm))
UV_Vis_data_exp = load (filename);
% finding row number corresponding to wavelength 400 nm, 617 nm & 1000
...wavelength_index is a subroutine, whose body is defined in the end
index_400 = wavelength_index(UV_Vis_data_exp (:,1),400);
index_618 = wavelength_index(UV_Vis_data_exp (:,1),618);
index_900 = wavelength_index(UV_Vis_data_exp (:,1),900);
%extracting data from UV-VIS spectum file
%LR stands for Longitudinal Resonance
absorb_exp_LR = UV_Vis_data_exp (index_900:index_618,2);
absorb_exp1_normalized = absorb_exp_LR./max(absorb_exp_LR);
absorb_exp_complete_spectrum_normalized = UV_Vis_data_exp(index_900:index_400,2)./max(UV_Vis_data_exp (index_900:index_400,2));
lambda_complete = UV_Vis_data_exp (index_900:index_400,1); %% lambda = 400 nm to 1000 nm
lambda_LR = UV_Vis_data_exp (index_900:index_618,1);
len_lamba = length(lambda_LR);
%interpolation of experimental refractive index data over 300 to 900 at step of 1

%calculating relative permittivity from interpolated refractive index data


%%%%%%%%%%%%%%%% Applying Gans theory %%%%%%%%%%%%%%


N = 1; % Initial contribution
V = 1; %volume

AR = 90.645/47.723;%aspect ratio
epsilon_medium = 2;%dielectric constant of surrouding medium 
 
 absorb_exp1_norm_new = absorb_exp1_normalized;

%%%% for plotting complete curve i.e lambda = 400 nm to 900 nm
%interpolation of experimental epsilon data over 300 to 900 at step of 1
ref_real_interpolated = interp1(wavelength_john, n_exp,lambda_complete,'spline');
ref_imag_interpolated = interp1(wavelength_john, k_exp,lambda_complete,'spline');
%calculating relative permittivity from ref. index data
epsilon_exp_complete = (ref_real_interpolated+1i.*ref_imag_interpolated).^2;
epsilon_real_interpolated = real(epsilon_exp_complete);
epsilon_imag_interpolated = imag(epsilon_exp_complete);
% intialization corresponding to complete lambda
%Absob_coff_transpose = zeros(length(AR),length(lambda_complete));
%Absob_coff_transpose_normalized = zeros(length(AR),length(lambda_complete));
% Calculating complete UV-Vis form fitted contribution

 elong = sqrt (1-(1/(AR^2)));
 PA = ((1-elong^2)/elong^2)*(1./(2*elong)* log ((1+elong)/(1-elong))-1);
 PB = (1-PA)/2;
 PC = PB;
 coff_A = ((1/PA^2).*epsilon_imag_interpolated)./((epsilon_real_interpolated + ((1-PA)/PA)*epsilon_medium).^2 + epsilon_imag_interpolated.^2);
 coff_B = ((1/PB^2).*epsilon_imag_interpolated)./((epsilon_real_interpolated + ((1-PB)/PB)*epsilon_medium).^2 + epsilon_imag_interpolated.^2);
 coff_C = coff_B;
 Absob_coff_transpose  = (2*pi*N*V*(epsilon_medium)^(3/2))./(3.*lambda_complete) .*(coff_A + coff_B +coff_C);
 Absob_coff_transpose_normalized  = Absob_coff_transpose  ./max(Absob_coff_transpose) ;

Absob_coff_norm = transpose(Absob_coff_transpose_normalized);
fitted_curve_1= Absob_coff_norm;
% calculating mean & standard deviation

%%%%% plotting figures %%%%%%%%%%%%%%
% plotting uv-vis spectrum
figure ();
p = plot (lambda_complete,absorb_exp_complete_spectrum_normalized,lambda_complete,fitted_curve_1,'Linewidth',6);
p(1).LineStyle = '-';
p(1).LineWidth = 1;
p(1).Color ='[0,0,1]';
p(2).LineStyle = '--';
p(2).LineWidth = 1;
p(2).Color ='[1,0,0]';
set ( gca , 'FontSize' , 30 , 'fontweight' , 'b' , 'FontName' , 'Times New Roman') ;
xlabel ( 'Wavelength (nm)' , 'FontSize' , 30 , 'fontweight' , 'b' ,'FontName' , 'Times New Roman') ;
ylabel ( 'Absorption (a.u.)', 'FontSize' , 30 , 'fontweight' , 'b', 'FontName' , 'Times New Roman' ) ;
legend ( 'Experimental' ,'Gans theory approximation') ;
xlim([400,1000]);
ylim([0,1]);
dimension = [0.2 0.5 0.3 0.3];

%%%% GOld refractive index were taken from following reference
% P.B. Johnson and R.W. Christy, “Optical constants of the noble metals,” Phys. Rev. B, vol. 6, pp. 4370, 1972.
% Data format [wavelength (nm), real (refractive index), img (refractive_index)]


function y = data_john_christy()
y = [983.968254 0.22 6.350
 891.942446 0.17 5.663
 821.0596026 0.16 5.083
 755.9756098 0.14 4.542
 704.4318182 0.13 4.103
 659.4680851 0.14 3.697
 616.8159204 0.21 3.272
 582.0657277 0.29 2.863
 548.5840708 0.43 2.455
 520.9243697 0.62 2.081
 495.9200000 1.04 1.833
 471.4068441 1.31 1.849
 450.8363636 1.38 1.914
 430.4861111 1.45 1.948
 413.2666667 1.46 1.958
 397.3717949 1.47 1.952
 381.4769231 1.46 1.933
 367.8931751 1.48 1.895
 354.2285714 1.5 1.866
 342.4861878 1.48 1.871
 331.4973262 1.48 1.883
 320.3617571 1.54 1.898
 310.726817 1.53 1.893
 300.9223301 1.53 1.889
 292.4056604 1.49 1.878];
end
% wavelength subroutine Body
function y = wavelength_index(wavelen, wave)
[~, index_400] = min (abs(wavelen(:,1)-wave));
y = index_400;
end

%reference
%https://www.nature.com/articles/s41598-019-53621-4#Equ1