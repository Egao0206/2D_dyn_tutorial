
function out = Process_results(model_name,Het_para,FltX,stress,dc,nor_str,...
    mus,mud,x,y,v,indx,vsat,Fault_length, OUTxseis,OUTyseis,time,OUTv,FltV,FltD,it)

% Source properties
STF = zeros(size(FltV,2),1);
for i=1:size(FltV,2)
    for p=2:length(FltX)
        STF(i) = STF(i) + FltV(p,i)*(FltX(p)-FltX(p-1));
    end
end

cuttent_tip = 0;
for p=2:length(FltX)
  if(FltD(p-1,it)>0.01 && FltD(p,it)<=0.01)
      current_tip = p;
      break;  
  end
end


clf
ax1 = subplot(3, 2, 1);
isel = find(abs(FltX)>=Fault_length);
stress(isel) = nan;
dc(isel) = nan;
if(strcmp(Het_para,'stress'))
    plot(FltX/1e3,stress/1e6);
    title('Initial shear stress')
    xlabel('Fault X (km)')
    ylabel('Shear stress (MPa)')
    hold on
    plot([0,Fault_length/1e3],[nor_str*mus/1e6, nor_str*mus/1e6],'k--')
    plot([0,Fault_length/1e3],[nor_str*mud/1e6, nor_str*mud/1e6],'k--')
    plot(FltX(current_tip)/1e3,stress(current_tip)/1e6,'r*')
    hold off
else
    plot(FltX/1e3,dc);
    title('Dc')
    xlabel('X')
    ylabel('Dc (m)')
    hold on
    plot(FltX(current_tip)/1e3,dc(current_tip),'r.');
    hold off
end
xlim([0,Fault_length/1e3]);


ax2 = subplot(3, 2, 3);
set(gca,'DefaultPatchEdgeColor','none');
h=patch(ax2,'Vertices',[x(:)/1e3 y(:)/1e3],'Faces',indx,'FaceVertexCData',v(:),'FaceColor','interp');
axis equal 
axis tight
title('SEM2D snapshot')
xlabel('X (km)')
ylabel('Y (km)')
if exist('vsat','var') && ~isempty(vsat), set(gca,'CLim',vsat); end
colorbar('vert')

hold(ax2,'on')
plot(OUTxseis/1e3,OUTyseis/1e3,'r^')
%plot(Fault_length/1e3,0,'r+');
rectangle('Position',[0 -0.5 Fault_length/1e3 0.5],'FaceColor','r')
hold(ax2,'off')


ax3 = subplot(3, 2, 5);
NT = length(time);
ampli = 0.5*max(abs( OUTv(:) ));
%disp(sprintf('Amplitude (trace to trace) = %f',ampli))
if ampli>0
  offset = max(abs(diff(OUTxseis)));
  ampli = offset/ampli;
end
plot(time, OUTv*ampli +repmat(OUTxseis,1,NT) );
title('SEM Seismograms')
xlabel('Time (s)')
ylabel('Distance (m)')


ax4 = subplot(3, 2, 2);
Sliprate = FltV(:,it);
for i=1:length(Sliprate)
    if(Sliprate(i)<1e-8)
       Sliprate(i)=0.0;
    end
end
plot(FltX/1e3,Sliprate);
title('Slip rate')
xlabel('Fault X (km)')
ylabel('Slip rate (m/s)')
xlim([0,Fault_length/1e3]);

ax5 = subplot(3, 2, 4);
plot(FltX/1e3,FltD(:,it));
title('Slip')
xlabel('Fault X (km)')
ylabel('Slip (m)')
xlim([0,Fault_length/1e3]);

ax6 = subplot(3, 2, 6);
plot(time,STF(1:size(time)));
title('Source time function')
xlabel('Time (s)')
ylabel('STF')


% Spec = abs(fft(STF)/length(STF));
% f = 1/(time(2)-time(1))*(0:(length(STF))-1)/length(STF);
% 
% ax6 = subplot(3, 2, 6);
% loglog(f(1:int64(length(f)/2)),Spec(1:int64(length(f)/2)));
% title('Spectrum of STF')
% xlabel('Frequency (Hz)')
% ylabel('Spectrum')
% 

