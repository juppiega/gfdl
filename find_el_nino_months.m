function [el_nino, la_nina] = find_el_nino_months(timestamps)

D = load('nino34.txt');
yr_data = D(:,1); yr_data = repmat(yr_data, 1,12)'; yr_data = yr_data(:);
mo_data = repmat((1:12)',size(D,1),1);
t_mat = zeros(length(mo_data),6); 
t_mat(:,1) = yr_data; t_mat(:,2) = mo_data; t_mat(:,3) = 1;
t_data = datenum(t_mat);
nino34 = D(:,2:end)'; 
nino34 = smooth(nino34(:),3);
nino34_norm = nino34(:) / std(nino34(end-12*40:end));
N = length(timestamps);

el_nino = false(N,1);
la_nina = false(N,1);
timestamps = timestamps + 0.5;

for i = 1:N
    mo = find(t_data < timestamps(i),1,'last');
    if nino34(mo) > 1
        el_nino(i) = true;
    elseif nino34(mo) < -1
        la_nina(i) = true;
    end
end

end