addpath('..\..\Lattice Boltzmann\Practice_MiniLBM\prototype_D2Q9_rectangular\test_cases')
couette_inc_data = readtable('cfi_rk4_test.csv');
couette_comp_data = readtable('cfc_rk4_test.csv');
%poiseuille_inc_data = readtable('poiseuille_incompressible_L2_error_sorted.csv');
%poiseuille_comp_data = readtable('poiseuille_compressible_L2_error_sorted.csv');

%%
%load('couette_compressible_L2_error.mat');
%load('couette_incompressible_L2_error.mat');
%load('poiseuille_compressible_L2_error.mat');
%load('poiseuille_incompressible_L2_error.mat');
%figure(1);
%  loglog(couette_inc_data.Re,couette_inc_data.L2,'.','Color',[0 0.4470 0.7410]);
%  hold on;
%  loglog(couette_comp_data.Re,couette_comp_data.L2,'.','Color',[0.8500 0.3250 0.0980]);
%  hold on;
%  loglog(poiseuille_inc_data.Re,poiseuille_inc_data.L2,'.','Color',[0.9290 0.6940 0.1250]);
%  hold on;
%  loglog(poiseuille_comp_data.Re,poiseuille_comp_data.L2,'.','Color',[0.4940 0.1840 0.5560]);
%  hold on;
%  grid on;
%  xline(1.5e15,'--r',{'Couette incompressible crash (1.5e15)'},'LabelHorizontalAlignment','left');
%  xline(8.34e7,'--r',{'Couette compressible crash (8.34e7)'},'LabelHorizontalAlignment','left');
%  xline(5.166e19,'--r',{'Poiseuille incompressible crash (5.166e19)'},'LabelHorizontalAlignment','left');
%  xline(918.11916,'--r',{'Poiseuille compressible crash (918.11916)'},'LabelHorizontalAlignment','left');
%  yline(10e-2,'--',{'10e-2 error'},'LabelHorizontalAlignment','left');
%  yline(10e-5,'--',{'10e-5 error'},'LabelHorizontalAlignment','left');
%  xlabel('Reynolds number');
%  ylabel('L_2 error');
%  legend('Couette, incompressible','Couette, compressible','Poiseuille, incompressible','Poiseuille, compressible','Location','southwest');

%figure(2);
%  loglog(couette_inc_data.Re,couette_inc_data.L2,'.','Color',[0 0.4470 0.7410]);
%  hold on;
%  loglog(couette_comp_data.Re,couette_comp_data.L2,'.','Color',[0.8500 0.3250 0.0980]);
%  hold on;
%  loglog(poiseuille_inc_data.Re,poiseuille_inc_data.L2,'.','Color',[0.9290 0.6940 0.1250]);
%  hold on;
%  loglog(poiseuille_comp_data.Re,poiseuille_comp_data.L2,'.','Color',[0.4940 0.1840 0.5560]);
%  hold on;
%  grid on;
%  yline(10e-2,'--',{'10e-2 error'},'LabelHorizontalAlignment','left');
%  yline(10e-5,'--',{'10e-5 error'},'LabelHorizontalAlignment','left');
%  xlabel('Reynolds number');
%  ylabel('L_2 error');
%  legend('Couette, incompressible','Couette, compressible','Poiseuille, incompressible','Poiseuille, compressible','Location','southwest')%

%figure(3);
%  tiledlayout(2,2);
  
%  nexttile;
%  loglog(couette_inc_data.Re,couette_inc_data.L2,'.','Color',[0 0.4470 0.7410]);
%  title('Couette, Incompressible');
%  xline(1.5e15,'-r',{'Couette incompressible', 'crash (1.5e15)'},'LabelHorizontalAlignment','left','LabelOrientation','horizontal');
%  yline(10e-2,'--',{'10e-2 error'},'LabelHorizontalAlignment','left');
%  yline(10e-5,'--',{'10e-5 error'},'LabelHorizontalAlignment','left');
%  xlabel('Reynolds number');
%  ylabel('L_2 error');
%  grid on;

%  nexttile;
%  loglog(couette_comp_data.Re,couette_comp_data.L2,'.','Color',[0.8500 0.3250 0.0980]);
%  title('Couette, Compressible');
%  xline(8.34e7,'-r',{'Couette compressible', 'crash (8.34e7)'},'LabelHorizontalAlignment','left','LabelOrientation','horizontal');
%  yline(10e-2,'--',{'10e-2 error'},'LabelHorizontalAlignment','left');
%  yline(10e-5,'--',{'10e-5 error'},'LabelHorizontalAlignment','left');
%  xlabel('Reynolds number');
%  ylabel('L_2 error');
%  grid on;

%  nexttile;
%  loglog(poiseuille_inc_data.Re,poiseuille_inc_data.L2,'.','Color',[0.9290 0.6940 0.1250]);
%  title('Poiseuille, Incompressible');
%  xline(5.166e19,'-r',{'Poiseuille incompressible', 'crash (5.166e19)'},'LabelHorizontalAlignment','left','LabelOrientation','horizontal');
%  yline(10e-2,'--',{'10e-2 error'},'LabelHorizontalAlignment','left');
%  yline(10e-5,'--',{'10e-5 error'},'LabelHorizontalAlignment','left');
%  xlabel('Reynolds number');
%  ylabel('L_2 error');
%  grid on;

%  nexttile;
%  loglog(poiseuille_comp_data.Re,poiseuille_comp_data.L2,'.','Color',[0.4940 0.1840 0.5560]);
%  title('Poiseuille, Compressible');
%  xline(918.11916,'-r',{'Poiseuille compressible', 'crash (918.11916)'},'LabelHorizontalAlignment','left','LabelOrientation','horizontal');
%  yline(10e-2,'--',{'10e-2 error'},'LabelHorizontalAlignment','left');
%  yline(10e-5,'--',{'10e-5 error'},'LabelHorizontalAlignment','left');
%  xlabel('Reynolds number');
%  ylabel('L_2 error');
%  grid on;
