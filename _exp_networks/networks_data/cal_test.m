net_id=10;

Net=f_read_Net(net_id);
for ii=1:100
    try
        load(strcat('../', num2str(ii)));
        f_wt_method(Net, ii);
%         clear mex;
    catch
        0;
    end
end







