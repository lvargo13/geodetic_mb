% rgiid,period,area,dhdt,err_dhdt,dvoldt,err_dvoldt
% read in data
t1 = readtable('EOSS_RGI.csv'); 
t2 = readtable('dh_18_rgi60_pergla_rates.csv'); 
rgi_all = t2.rgiid; 
years = t2.period; 
dhdt_all = t2.dhdt; %  t2.dhdt or t2.dvoldt
dhdt_all_err = t2.err_dhdt; % t2.err_dhdt or t2.err_dvoldt


% pull out data
dhdt = zeros(height(t1),20);
dhdt_err = zeros(height(t1),20);
y1 = 1999;
y2 = 2000;

for i =1:height(t1)
    rgi_inds = find(strcmp(t1.RGIId{i},rgi_all)); 
    period = years(rgi_inds(1):rgi_inds(end)); 
    dh = dhdt_all(rgi_inds(1):rgi_inds(end));
    dhe = dhdt_all_err(rgi_inds(1):rgi_inds(end));
    
    for ii = 1:20
        st_want = [num2str(y1+ii) '-01-01_' num2str(y2+ii) '-01-01']; 
        tmp = find(strcmp(st_want,period)); 
        dhdt(i,ii) = dh(tmp);
        dhdt_err(i,ii) = dhe(tmp); 
    end
end

dhdt_eoss = sum(dhdt); 
dhdt_err_eoss = sum(dhdt_err); 

% plot
figure; errorbar(2001:2020, dhdt_eoss, dhdt_err_eoss,'.','LineWidth',2.5) ; hold on
plot(2001:2020, dhdt_eoss, 'ko--'); 
xlabel('years'); ylabel('dVol (m3)'); title('EOSS all annual volume change')

figure; imagesc(dhdt)
set(gca, 'YTick',1:40, 'YTickLabel',t1.name)   
set(gca, 'XTick',1:20, 'XTickLabel',2001:2020)
colorbar; title('EOSS all volume change (m)')

figure; imagesc(dhdt_err)
set(gca, 'YTick',1:40, 'YTickLabel',t1.name)   
set(gca, 'XTick',1:20, 'XTickLabel',2001:2020)
colorbar; title('EOSS all volume change (m) uncertainty')

% dV Brewster
dhdt_brew = dhdt(26,:);
dhdt_err_brew = dhdt_err(26,:);

figure; errorbar(2001:2020, dhdt_brew, dhdt_err_brew,'.','LineWidth',2.5) ; hold on
plot(2001:2020, dhdt_brew, 'ko--'); 
xlabel('years'); ylabel('dVol (m3)'); title('Brewster annual volume change (m3)')



% witout Tasman
dhdt(16,:) = [];
dhdt_eoss_subTas = mean(dhdt); 
dhdt_err(16,:) = [];
dhdt_err_eoss_subTas = mean(dhdt_err); 
nm = t1.name([1:15,17:end]); 

figure; errorbar(2001:2020, dhdt_eoss_subTas, dhdt_err_eoss_subTas,'.','LineWidth',2.5) ; hold on
plot(2001:2020, dhdt_eoss_subTas, 'ko--'); 
xlabel('years'); ylabel('dVol (m3)'); title('EOSS (minus Tasman) annual volume change (m3)')

figure; imagesc(dhdt)
set(gca, 'YTick',1:39, 'YTickLabel',nm)   
set(gca, 'XTick',1:20, 'XTickLabel',2001:2020)
colorbar; title('EOSS (minus Tasman) volume change (m3)')

figure; imagesc(dhdt_err)
set(gca, 'YTick',1:39, 'YTickLabel',nm)   
set(gca, 'XTick',1:20, 'XTickLabel',2001:2020)
colorbar; title('EOSS (minus Tasman) volume change (m3) uncertainty')
