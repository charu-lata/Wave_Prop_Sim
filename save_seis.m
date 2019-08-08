%%%saves seisx seisz hdr at the end of total time of propagation


save_dir = '../v3/OUT/seis.mat';

%%%%%variables
seisx = Ut';
seisz = Wt';

for i=1:trace.num_traces
%    hdr(i).channel = trace.channel(i);
%    hdr(i).vred = red_vel;
   hdr(i).x = trace.x(i);
   hdr(i).z = trace.z(i);
   hdr(i).range = trace.range(i);
   hdr(i).t_window = trace.tlen;
   hdr(i).nsamp = trace.nt;
   hdr(i).shot = i;
   hdr(i).static = 0;
   hdr(i).sampint = dt;
end
eval(['save ',save_dir,' hdr seisx seisz velp vels rhog xaxis zaxis'])
