function animate_frame(A,t,options)
arguments
    A (3,3,:) double 
    t (1,:) double 
    options.SpeedUp (1,1) = 1 %default to animate every 10 timesteps 
    options.ShowTime (1,1) logical = true %default to not showing time x speed
end
% given A attitude marix solution (3,3,length(t)), animate rotating frame
% wrt inertial fixed frame (A mst trabsform FROM INERTIAL TO BODY!)


% inertial frame (fixed)
figure('Name','Animation of rotating frame '); 
hold on 
grid on
axis equal

quiver3(0,0,0,1,0,0,'DisplayName','i','Color','blue','LineWidth',1.5); 
quiver3(0,0,0,0,1,0,'DisplayName','j','Color','blue','LineWidth',1.5); 
quiver3(0,0,0,0,0,1,'DisplayName','k','Color','blue','LineWidth',1.5); 

text(1,0,0,'i','Color', 'blue'); 
text(0,1,0,'j','Color', 'blue'); 
text(0 ,0,1,'k','Color', 'blue'); 

%rotating frame initialization 
bx0=A(1,:,1); 
by0=A(2,:,1); 
bz0=A(3,:,1); 
bx=quiver3(0,0,0,bx0(1), bx0(2), bx0(3),"DisplayName",'b_x','Color','red','LineWidth',1.5); 
by=quiver3(0,0,0,by0(1), by0(2), by0(3),"DisplayName",'b_y','Color', 'red','LineWidth',1.5); 
bz=quiver3(0,0,0,bz0(1), bz0(2), bz0(3),"DisplayName",'b_z','Color','red','LineWidth',1.5); 
txt_bx=text(bx0(1), bx0(2), bx0(3), 'b_x', 'FontSize', 12, 'Color', 'red');
txt_by=text(by0(1), by0(2), by0(3), 'b_y', 'FontSize', 12, 'Color', 'red');
txt_bz=text(bz0(1), bz0(2), bz0(3), 'b_z', 'FontSize', 12, 'Color', 'red');

%set legend and limit for axis
legend show; 
xlim([-1.2 1.2]);
ylim([-1.2 1.2]); 
zlim([-1.2 1.2]);
view(3); %default to 3dview otherwise bad shit happens god knows why 

%animation (kind of)
step=options.SpeedUp; 
len=length(t); 

for i = 2:step:length(t)

   
    R=A(:,:,i); %instantaneus rotattion mat
    %new directions 
    x_b=R(1,:); 
    y_b=R(2,:); 
    z_b=R(3,:); 

    %update vectors and text 
    set(bx, 'UData',x_b(1), 'VData', x_b(2) , 'WData', x_b(3)); 
    set(by, 'UData',y_b(1), 'VData', y_b(2) , 'WData', y_b(3)); 
    set(bz, 'UData',z_b(1), 'VData', z_b(2) , 'WData', z_b(3)); 
    set(txt_bx, 'Position' , x_b); 
    set(txt_by, 'Position', y_b); 
    set(txt_bz, 'Position', z_b); 
    if options.ShowTime
    title(sprintf('TIME : %.2f' , t(i))); 
    end 
   
    drawnow limitrate; %this works better because chat told me so... no clue why  
    pause(0.05); 
    
    if (step+i) > len %to suppres index out of bounds...
        break; 
    end
end
    

end 