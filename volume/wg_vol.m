[tsg,t,z,r,b,tit,param,planes,traces,samples,fs]=pp_reader('vol_scat.asc');
tim=[0.0:1/fs:0.0+(samples-1)/fs];
beam = [-90:90];
[tsg_beam]=beamform(tsg,z,fs,beam,1500); 
[a,b]=wavei(dba(tsg_beam)',tim,beam,-60,0);
title('Volume Scattering -  Waveguide');
xlabel('Time (s)');
ylabel('Beam (deg)');
