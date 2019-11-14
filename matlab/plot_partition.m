clear all; close all

data = importdata('../partition.txt');
np = size(data,1)/6;
N = max(data(:,3));


figure(1)
hold on
panel_num = 1;
for id = (panel_num-1)*np+1:panel_num*np
   is = data(id,2); 
   ie = data(id,3);
   js = data(id,4);
   je = data(id,5);
%    fill([is ie ie is], [js js je je],'r', 'LineStyle','none')
   line([is ie], [js js],'Color','k')
   line([is is], [js je],'Color','k')
   text(is+N/100, js+N/100, num2str(id-1))
   
%    line([is ie], [je je],'Color','k')
%    line([ie,ie], [js,je],'Color','k')
end
line([1 N], [N N],'Color','k')
line([N N], [1 N],'Color','k')

panel_num = 2;
for id = (panel_num-1)*np+1:panel_num*np
   is = data(id,2)+N+0.1*N; 
   ie = data(id,3)+N+0.1*N;
   js = data(id,4);
   je = data(id,5);
   line([is ie], [js js],'Color','k')
   line([is is], [js je],'Color','k')
   text(is+N/100, js+N/100, num2str(id-1))
%    line([is ie], [je je],'Color','k')
%    line([ie,ie], [js,je],'Color','k')
end
line([1 N] +N+0.1*N, [N N],'Color','k')
line([N N] +N+0.1*N, [1 N],'Color','k')

panel_num = 3;
for id = (panel_num-1)*np+1:panel_num*np
   is = data(id,2)+2*(N+0.1*N); 
   ie = data(id,3)+2*(N+0.1*N);
   js = data(id,4);
   je = data(id,5);
   line([is ie], [js js],'Color','k')
   line([is is], [js je],'Color','k')
   text(is+N/100, js+N/100, num2str(id-1))
%    line([is ie], [je je],'Color','k')
%    line([ie,ie], [js,je],'Color','k')
end
line([1 N] +2*(N+0.1*N), [N N],'Color','k')
line([N N] +2*(N+0.1*N), [1 N],'Color','k')

panel_num = 4;
for id = (panel_num-1)*np+1:panel_num*np
   is = data(id,2)-1*(N+0.1*N); 
   ie = data(id,3)-1*(N+0.1*N);
   js = data(id,4);
   je = data(id,5);
   line([is ie], [js js],'Color','k')
   line([is is], [js je],'Color','k')
   text(is+N/100, js+N/100, num2str(id-1))   
%    line([is ie], [je je],'Color','k')
%    line([ie,ie], [js,je],'Color','k')
end
line([1 N] -1*(N+0.1*N), [N N],'Color','k')
line([N N] -1*(N+0.1*N), [1 N],'Color','k')

panel_num = 5;
for id = (panel_num-1)*np+1:panel_num*np
   is = data(id,2);%+N+0.1*N; 
   ie = data(id,3);%+N+0.1*N;
   js = data(id,4)-N-0.1*N;
   je = data(id,5)-N-0.1*N;
   line([is ie], [js js],'Color','k')
   line([is is], [js je],'Color','k')
   text(is+N/100, js+N/100, num2str(id-1))   
%    line([is ie], [je je],'Color','k')
%    line([ie,ie], [js,je],'Color','k')
end
line([1 N]             , [N N]-N-0.1*N,'Color','k')
line([N N]             , [1 N]-N-0.1*N,'Color','k')

panel_num = 6;
for id = (panel_num-1)*np+1:panel_num*np
   is = data(id,2);%+N+0.1*N; 
   ie = data(id,3);%+N+0.1*N;
   js = data(id,4)+N+0.1*N;
   je = data(id,5)+N+0.1*N;
   line([is ie], [js js],'Color','k')
   line([is is], [js je],'Color','k')
   text(is+N/100, js+N/100, num2str(id-1))   
%    line([is ie], [je je],'Color','k')
%    line([ie,ie], [js,je],'Color','k')
end
line([1 N] , [N N]+N+0.1*N,'Color','k')
line([N N] , [1 N]+N+0.1*N,'Color','k')

axis([-1.2*N 3*N+3*0.1*N -N-2*0.1*N 2*N+2*0.1*N])

% figure(2)
% hold on
% 
% panel_num = 1;
% for id = (panel_num-1)*np+1:panel_num*np
%    is = data(id,2); 
%    ie = data(id,3);
%    js = data(id,4);
%    je = data(id,5);
%    line([is ie], [js js], [1 1],'Color','k')
%    line([is is], [js je], [1 1],'Color','k')
%    line([is ie], [je je], [1 1],'Color','k')
%    line([ie,ie], [js,je], [1 1],'Color','k')
% end
% panel_num = 2;
% for id = (panel_num-1)*np+1:panel_num*np
%    is = data(id,2); 
%    ie = data(id,3);
%    js = data(id,4);
%    je = data(id,5);
%    line([N N], [js js],[is ie], 'Color','k')
%    line([N N], [js je],[is is], 'Color','k')
%    line([N N], [je je],[is ie], 'Color','k')
%    line([N N], [js,je],[ie,ie], 'Color','k')
% end
% 
% panel_num = 3;
% for id = (panel_num-1)*np+1:panel_num*np
%    is = data(id,2); 
%    ie = data(id,3);
%    js = data(id,4);
%    je = data(id,5);
%    line([is ie], [js js], [N N],'Color','k')
%    line([is is], [js je], [N N],'Color','k')
%    line([is ie], [je je], [N N],'Color','k')
%    line([ie,ie], [js,je], [N N],'Color','k')
% end
% 
% panel_num = 4;
% for id = (panel_num-1)*np+1:panel_num*np
%    is = data(id,2); 
%    ie = data(id,3);
%    js = data(id,4);
%    je = data(id,5);
%    line([1 1], [js js],[is ie], 'Color','k')
%    line([1 1], [js je],[is is], 'Color','k')
%    line([1 1], [je je],[is ie], 'Color','k')
%    line([1 1], [js,je],[ie,ie], 'Color','k')
% end
% 
% axis([-N 2*N -N 2*N -N 2*N])


