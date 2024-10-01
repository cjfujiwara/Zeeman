function B_out = K_rf_freq2field(f_in)

[out,hF_S,hF_P1,hF_P3]=K_zeeman;

close(hF_S);
close(hF_P1);
close(hF_P3);

Es12=out.Es12;
dE=Es12(2,:)-Es12(1,:);
B=out.B;


Brange=[179 205];

i1=find(B>=Brange(1),1);
i2=find(B>=Brange(2),1);

B=B(i1:i2);
dE=dE(i1:i2);

% hF=figure;
% plot(B,dE)

B_out=interp1(dE,B,f_in);


end

