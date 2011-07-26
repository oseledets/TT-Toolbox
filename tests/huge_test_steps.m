
% huge_norms_Au = cell(5,6);
huge_results = cell(6,6);
huge_psi = cell(6,6);
huge_eta = cell(6,6);
for d0t=4:9
    for d0x=4:9
%         test_steps;
	test_full_KN;
        
%        huge_norms_Au{d0t-5,d0x-3} = results;
        huge_eta{d0t-3,d0x-3} = eta;
        huge_psi{d0t-3,d0x-3} = psi;
        huge_results{d0t-3,d0x-3} = global_results;
    end
end;
