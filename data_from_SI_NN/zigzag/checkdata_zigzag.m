clear 
close all
clc

%%%%%%%%%%%%%%%%%%%%
% quick check program for free-runnning experiment data, process and filter measured data, and 
% save to csv file 
%  input:    exp data spreadsheet (.xlsx). analyze all files under './alldata/'
% output: figures check for measured result and data processing
%         csv file contain: location, ship velocity, heading, angilar velo,
%         rudder angle, prop rev, wind measured on ship
%         rev, 
%%%%%%%%%%%%%%%%%%%%%

%%%% chose best condition GNSS for postporcess
%%% TODO: this is foolish implementation, should be fix to choose best one
%%% automatic
chosenGNSS = 1; % enter GPS to output speed and loc of midship
                    % 1: Fore, 2: Aft, 3:3rd
useheading = 2; %  enter the combination of GPS or gyro for compute heading of ship
                         % 0:Gyro, 1:Fore to aft, 2: Fore to 3rd, 3:aft to 3rd (not reccomended because of ship ref. length)
useinitheading =2; %enter the combination of GPS for compute initial heading of ship
                         %1:Fore to aft, 2: Fore to 3rd, 3:aft to 3rd (not reccomended because of ship ref. length) 
%%% params of berth
theta_berth = 1.143942290114601e+02; %longitudinal angle to north !!this should be acculate, necessary to corelate GPS and gyro!!
lon_origin    = 135+31.484779127149817/60; % origin of longitude of berth (measured on 2020/6/17 1502 400-12000sec, sigma=0.0126m)
lat_origin    = 34+ 49.365327040533202/60; % origin of latetude of berth (measured on 2020/6/17 1502 400-12000sec, sigma=0.0272m)

%%% params of GPS
L_GPS  = [1.432 -0.735 -1.160] ; %longi distance from  midship to GNSS antenna, forward +
no_GPS = 3; % number of GNSS sensors used

%%% params of Anemometer
L_anemo = [0.859 0.141]; %longitudinal distans from midship to anemometer, foward +, [Anemo_fore Anemo_mid]
no_anemo = 2; % number of Anemometers used


delta_t = 1; %%% params of numerical differential (dt:second)
init_t  = 2; % priod to initialize heading by 2 GPS (sec)

%%% set default for plot
set(0,'defaultAxesFontSize',11);
set(0,'defaultAxesFontName','times new roman');
set(0,'defaultAxesXcolor','k');
set(0,'defaultAxesYcolor','k');
set(0,'defaultAxesZcolor','k');
set(0,'defaultTextFontSize',11);
set(0,'defaultTextFontName','times new roman');
set(0,'defaultTextColor','k');
set(0,'defaultLineLineWidth',1.2);
set(0,'defaultLegendFontSize',11);
set(0,'defaultLegendFontName','times new roman');
set(0,'DefaultAxesXGrid','on');
set(0,'DefaultAxesYGrid','on');

%%% path of directory 
files = dir("./rawdata/*.xlsx");
for i = 1:length(files)
%  for i = 8
%     read experiment rawdata
    xfileplace = ["rawdata/", files(i).name];
    fileplace = xfileplace(1)+xfileplace(2);
    data = readtable(fileplace);
%     
%     %initialize values
    ave_spe = zeros(height(data),1);
    dist_f2a_m  = zeros(height(data),1);
    dist_f23_m  = zeros(height(data),1);
    dist_a23_m  = zeros(height(data),1);
    
    course_f2a  = zeros(height(data),1);
    course_f23  = zeros(height(data),1);
    course_a23  = zeros(height(data),1);
    
    trav_fore_m = zeros(height(data),1);
    trav_aft_m = zeros(height(data),1);
    trav_3_m = zeros(height(data),1);
    
    course_fore = zeros(height(data),1);
    course_aft = zeros(height(data),1);
    course_3 = zeros(height(data),1);
     
    psi     = zeros(height(data),1);
    psi_hat = zeros(height(data),1);
    psi_hat_gyro = zeros(height(data),1);
    psi_hat_GNSS = zeros(height(data),no_GPS);
    trav_GPS = zeros(height(data),no_GPS);
    U_GPS = zeros(height(data),no_GPS);
    U_GPS_kf = zeros(height(data),no_GPS);
    U_GPS_lowpass = zeros(height(data),no_GPS);
    U_m = zeros(height(data),no_GPS);
    U_m_kf = zeros(height(data),no_GPS);
    U_m_lowpass = zeros(height(data),no_GPS);
   
    x_north_GPS = zeros(height(data),no_GPS);
    x_hat_GPS = zeros(height(data),no_GPS);
    x_hat = zeros(height(data),no_GPS);
    x_hat_diff = zeros(height(data),no_GPS);
    
    u_hat = zeros(height(data),no_GPS);
    u_hat_GPS = zeros(height(data),no_GPS);
    u_hat_kf = zeros(height(data),no_GPS);
    u_hat_GPS_kf = zeros(height(data),no_GPS);
    u_hat_lowpass = zeros(height(data),no_GPS);
    u_hat_GPS_lowpass = zeros(height(data),no_GPS);
    u_shipfix_kf = zeros(height(data),no_GPS);
    u_shipfix_lp = zeros(height(data),no_GPS);
    u_shipfix_raw = zeros(height(data),no_GPS);
    
    y_north_GPS = zeros(height(data),no_GPS);
    y_hat_GPS = zeros(height(data),no_GPS);
    y_hat = zeros(height(data),no_GPS);
    y_hat_diff = zeros(height(data),no_GPS);
    
    vm_hat = zeros(height(data),no_GPS);
    vm_hat_GPS = zeros(height(data),no_GPS);
    vm_hat_kf = zeros(height(data),no_GPS);
    vm_hat_GPS_kf = zeros(height(data),no_GPS);
    vm_hat_lowpass = zeros(height(data),no_GPS);
    vm_hat_GPS_lowpass = zeros(height(data),no_GPS);
    vm_shipfix_kf = zeros(height(data),no_GPS);
    vm_shipfix_lp = zeros(height(data),no_GPS);
    vm_shipfix_raw = zeros(height(data),no_GPS);
    
    beta_kf  = zeros(height(data),no_GPS);
    beta_lp  = zeros(height(data),no_GPS);
    
    dir_of_travel_GPS = zeros(height(data),no_GPS);
    dir_of_travel_GPS_diff =zeros(height(data),no_GPS);
    
    gamma_anemo= zeros(height(data),no_anemo);
    UA_anemo    = zeros(height(data),no_anemo);
    gamma_anemo_check= zeros(height(data),no_anemo);
    UA_anemo_check    = zeros(height(data),no_anemo);
    u_shipfix_anemo =zeros(height(data),no_anemo);
    v_shipfix_anemo = zeros(height(data),no_anemo);
    U_anemo =zeros(height(data),no_anemo);
    beta_anemo =zeros(height(data),no_anemo);
    UT_anemo = zeros(height(data),no_anemo);
    zeta_wind_anemo = zeros(height(data),no_anemo);
    UA_midship =zeros(height(data),1);
    gamma_a_midship = zeros(height(data),1);
    
    %%%%
    % comvert dd.mm to dd.dddd for lat & long
    lat_fore = data.LatitudeDeg_Fore(:)+data.LatitudeMin_Fore(:)/60;
    lon_fore = data.LongitudeDeg_Fore(:)+data.LongitudeMin_Fore(:)/60;

    lat_aft = data.LatitudeDeg_Aft(:)+data.LatitudeMin_Aft(:)/60;
    lon_aft = data.LongitudeDeg_Aft(:)+data.LongitudeMin_Aft(:)/60;
        
    lat_3 = data.LatitudeDeg_3(:)+data.LatitudeMin_3(:)/60;
    lon_3 = data.LongitudeDeg_3(:)+data.LongitudeMin_3(:)/60;
    
   %%
    %measurement period.-
    dt  = round(data.time(2)-data.time(1),2);
        
    %%% parameter of Kalman filter
    if dt == 0.10
        A = 1.0; b = 1.0; c = 1; % System
        Qa = 0.05; Ra = 5; % Assumed Noise
    
    elseif dt ==0.01
        A = 1.0; b = 1.0; c = 1; % System
        Qa = 0.01; Ra = 10;    
    end
    
    % numerical differential period (sec)
    delay = round(delta_t/dt);
    
    %%% compute x-y loc and speed for each timestep from GPS lat and long
    for j =1:height(data)
      
       % compute difference of location
        gps_long_f2a =[lon_fore(j), lon_aft(j)];
        gps_lat_f2a  =[lat_fore(j), lat_aft(j)];
        
        gps_long_f23 =[lon_fore(j), lon_3(j)];
        gps_lat_f23  =[lat_fore(j), lat_3(j)];
        
        gps_long_a23 =[lon_aft(j), lon_3(j)];
        gps_lat_a23  =[lat_aft(j), lat_3(j)];

        lon_traveled_fore =[lon_origin, lon_fore(j)]; 
        lon_traveled_aft  =[lon_origin, lon_aft(j)];  
        lon_traveled_3  =[lon_origin, lon_3(j)];   
        lat_traveled_fore =[lat_origin, lat_fore(j)];
        lat_traveled_aft  =[lat_origin, lat_aft(j)];  
        lat_traveled_3  =[lat_origin, lat_3(j)];  
        
        % distance between GPS sensor
        [ dist_f2a, co ] = Cal_dist_co( gps_long_f2a, gps_lat_f2a );
        dist_f2a_m(j) = dist_f2a*1852;
        course_f2a(j)  = deg2rad(co);
                
        [ dist_f23, co ] = Cal_dist_co( gps_long_f23, gps_lat_f23 );
        dist_f23_m(j) = dist_f23*1852;
        course_f23(j)  = deg2rad(co);
        
        [ dist_a23, co ] = Cal_dist_co( gps_long_a23, gps_lat_a23 );
        dist_a23_m(j) = dist_a23*1852;
        course_a23(j)  = deg2rad(co);
        
       % distance traveled from t=0 for FORE GPS sensor
        [ trav_fore, co_fore ] = Cal_dist_co( lon_traveled_fore, lat_traveled_fore );
        trav_fore_m(j) = trav_fore*1852;
        course_fore(j) = deg2rad(co_fore);
         
        % distance traveled from t=0 for AFT GPS sensor
        [ trav_aft, co_aft ] = Cal_dist_co( lon_traveled_aft, lat_traveled_aft );
        trav_aft_m(j) = trav_aft*1852;
        course_aft(j) = deg2rad(co_aft);
        
        % distance traveled from t=0 for 3rd GPS sensor
        [ trav_3, co_3 ] = Cal_dist_co( lon_traveled_3, lat_traveled_3 );
        trav_3_m(j) = trav_3*1852;
        course_3(j) =deg2rad(co_3);
    end
    
         % corelate GPS heading and gyro heading
        init_heading_GPS(1) = mean(course_f2a(2:init_t/dt))-pi;
        init_heading_GPS(2) = mean(course_f23(2:init_t/dt))-pi; 
        init_heading_GPS(3) = mean(course_a23(2:init_t/dt))-pi; 
        init_heading_to_north = rad2deg(init_heading_GPS)
        
         % heading by GNSS1
        psi_hat_GNSS(:,1) =  course_f2a-pi-deg2rad(theta_berth);
        psi_hat_GNSS(:,2) =  course_f23-pi-deg2rad(theta_berth);
        psi_hat_GNSS(:,3) =  course_a23-pi-deg2rad(theta_berth); 
     
    for j =1:height(data)
        %%% compute U, x-y position, drift 
        psi(j) = deg2rad(data.YawAng(j)); %radian
        psi_hat_gyro(j) = psi(j) -mean(psi(1:10))+ init_heading_GPS(useinitheading) - deg2rad(theta_berth); 
        if useheading == 0
            psi_hat(j) = psi_hat_gyro(j);
        elseif useheading ==1
            psi_hat(j)=  course_f2a(j)-pi-deg2rad(theta_berth);
        elseif useheading ==2
            psi_hat(j)=  course_f23(j)-pi-deg2rad(theta_berth);
        elseif useheading ==3
            psi_hat(j)=  course_a23(j)-pi-deg2rad(theta_berth);
        end
        
        % convert to [-pi, pi]
        if psi(j) > pi
            psi(j) = psi(j)- 2*pi;
        elseif psi(j) < -pi
            psi(j) = psi(j) + 2*pi;
        end
        if psi_hat(j) > pi
            psi_hat(j) = psi_hat(j)-2*pi;
        elseif psi_hat(j) < -pi
            psi_hat(j) = psi_hat(j) + 2*pi;
        end
        if psi_hat_gyro(j) > pi
            psi_hat_gyro(j) = psi_hat_gyro(j)-2*pi;
        elseif psi_hat_gyro(j) < -pi
            psi_hat_gyro(j) = psi_hat_gyro(j) + 2*pi;
        end
        
    end
    
    for k =1:no_GPS           
       for j =1:height(data)
            % convert to [-pi, pi]
            if psi_hat_GNSS(j,k) > pi
                psi_hat_GNSS(j,k) = psi_hat_GNSS(j,k)- 2* pi;
            elseif psi_hat_GNSS(j,k) < -pi
                psi_hat_GNSS(j,k) = psi_hat_GNSS(j,k) + 2*pi;
            end
            
            % GPS position on x -> north, y -> east coordinate
            if k==1
                x_north_GPS(j,k) = trav_fore_m(j) * cos(course_fore(j));
                y_north_GPS(j,k) = trav_fore_m(j) * sin(course_fore(j));
            elseif k==2
                x_north_GPS(j,k) = trav_aft_m(j) * cos(course_aft(j));
                y_north_GPS(j,k) = trav_aft_m(j) * sin(course_aft(j));
            elseif k==3
                x_north_GPS(j,k) = trav_3_m(j) * cos(course_3(j));
                y_north_GPS(j,k) = trav_3_m(j) * sin(course_3(j));
            end
            
           % convert to berth coordinate (x-hat,y-hat, berth longitudinal direction is x)
            x_hat_GPS(j,k) = x_north_GPS(j,k) * cos(deg2rad(theta_berth))...
                +y_north_GPS(j,k) * sin(deg2rad(theta_berth)); 
            y_hat_GPS(j,k) = -x_north_GPS(j,k) * sin(deg2rad(theta_berth))...
                +y_north_GPS(j,k) * cos(deg2rad(theta_berth));    
            % convert location GPS 2 midship 
            x_hat(j,k) = x_hat_GPS(j,k) - L_GPS(k)*cos(psi_hat(j));
            y_hat(j,k) = y_hat_GPS(j,k) - L_GPS(k)*sin(psi_hat(j));  
            

       end
    end
    for k = 1:no_GPS
        for j =1:height(data)
                % difference of MS location between antenna ref. to fore GPS(k=1):TODO(selet best gps automatically)
            if k==3
                x_hat_diff(j,k) = abs(x_hat(j,k)-x_hat(j,1));
                y_hat_diff(j,k) = abs(y_hat(j,k)-y_hat(j,1));
            elseif k==2 || k ==1
                x_hat_diff(j,k) = abs(x_hat(j,k)-x_hat(j,k+1));
                y_hat_diff(j,k) = abs(y_hat(j,k)-y_hat(j,k+1));    
            end
        end
    end
    for k =1:no_GPS 
        for j =1:height(data)
             % numerical deffentitation for ship velocity
            if delay < j && j < height(data)-delay % 2nd order central
                u_hat(j,k) = (x_hat(j+delay,k)-x_hat(j-delay,k))/(2*delta_t);
                vm_hat(j,k) = (y_hat(j+delay,k)-y_hat(j-delay,k))/(2*delta_t);
                U_m(j,k)  = sqrt(u_hat(j,k)^2 + vm_hat(j,k)^2);
                u_hat_GPS(j,k) = (x_hat_GPS(j+delay,k)-x_hat_GPS(j-delay,k))/(2*delta_t);
                vm_hat_GPS(j,k) = (y_hat_GPS(j+delay,k)-y_hat_GPS(j-delay,k))/(2*delta_t);
                U_GPS(j,k)  = sqrt(u_hat_GPS(j,k)^2 + vm_hat_GPS(j,k)^2);
            
            elseif j <= delay % 1st order upwind
                u_hat(j,k) = (x_hat(j+delay,k)-x_hat(j,k))/delta_t;
                vm_hat(j,k) = (y_hat(j+delay,k)-y_hat(j,k))/delta_t;
                U_m(j,k)  = sqrt(u_hat(j,k)^2 + vm_hat(j,k)^2);
                u_hat_GPS(j,k) = (x_hat_GPS(j+delay,k)-x_hat_GPS(j,k))/delta_t;
                vm_hat_GPS(j,k) = (y_hat_GPS(j+delay,k)-x_hat_GPS(j,k))/delta_t;
                U_GPS(j,k)  = sqrt(u_hat_GPS(j,k)^2 + vm_hat_GPS(j,k)^2);
            
            elseif j >= height(data)-delay % 1st order downwind
                u_hat(j,k) = (x_hat(j,k)-x_hat(j-delay,k))/delta_t;
                vm_hat(j,k) = (y_hat(j,k)-y_hat(j-delay,k))/delta_t;
                U_m(j,k)  = sqrt(u_hat(j,k)^2 + vm_hat(j,k)^2);
                u_hat_GPS(j,k) = (x_hat_GPS(j,k)-x_hat_GPS(j-delay,k))/delta_t;
                vm_hat_GPS(j,k) = (y_hat_GPS(j,k)-y_hat_GPS(j-delay,k))/delta_t;
                U_GPS(j,k)  = sqrt(u_hat_GPS(j,k)^2 + vm_hat_GPS(j,k)^2);     
            end
        end
    end
   
    %%% kalman filtering of speed
    xhat = zeros(height(data), no_GPS);
    y  = zeros(height(data), no_GPS);
   for k = 1: no_GPS
    for m = 1 : 6
        if m == 1
            y = u_hat(:,k);
        elseif m==2
            y = vm_hat(:,k);
        elseif m == 3
            y = U_m(:,k);
        elseif m == 4
            y = u_hat_GPS(:,k);
        elseif m == 5
            y = vm_hat_GPS(:,k);
        else           
            y = U_GPS(:,k);
        end 
        
    % initialize
        P=0;
        xhat = 0.0;
        xhat(1) = y(1);
    % update emtimations for each times step
        for l = 2 : height(data)
            [ xhat(l), P, G ] = kf( A, b, 0, c, Qa, Ra, 0, y( l ), xhat( l - 1 ), P);
        end
        
        if m == 1
            u_hat_kf(:,k) = xhat;
        elseif m ==2
            vm_hat_kf(:,k) =xhat;
        elseif m==3
            U_m_kf(:,k) =xhat;
        elseif m ==4
            u_hat_GPS_kf(:,k) =xhat;
        elseif m ==5
            vm_hat_GPS_kf(:,k) =xhat;
        else
            U_GPS_kf(:,k) =xhat;
        end
    end
   end
    % low pass filtering the speed 
    Fs = 1/dt;
    fil = 0.2; % filtering freq
    for k = 1:no_GPS
        u_hat_lowpass(:,k) = lowpass(u_hat(:,k), fil,Fs,'steepness',0.99);
        vm_hat_lowpass(:,k) = lowpass(vm_hat(:,k), fil,Fs,'steepness',0.99);
        U_m_lowpass(:,k)  = sqrt(u_hat_lowpass(:,k).^2+vm_hat_lowpass(:,k).^2);
        U_GPS_lowpass(:,k) = lowpass(U_GPS(:,k), fil,Fs,'steepness',0.99);
    end
    r_lowpass     = deg2rad(lowpass(data.YawAngVel,fil,Fs,'steepness',0.99));
    %%% convert speed to ship fix coordinate 
    for k = 1: no_GPS
        for j=1:height(data)
            u_shipfix_kf(j,k) = u_hat_kf(j,k) * cos(psi_hat(j)) + vm_hat_kf(j,k) * sin(psi_hat(j));
            vm_shipfix_kf(j,k) = vm_hat_kf(j,k) * cos(psi_hat(j)) - u_hat_kf(j,k) * sin(psi_hat(j));
            beta_kf(j,k) = atan2(vm_shipfix_kf(j,k), u_shipfix_kf(j,k));
        
            u_shipfix_lp(j,k) = u_hat_lowpass(j,k) * cos(psi_hat(j)) + vm_hat_lowpass(j,k) * sin(psi_hat(j));
            vm_shipfix_lp(j,k) = vm_hat_lowpass(j,k) * cos(psi_hat(j)) - u_hat_lowpass(j,k) * sin(psi_hat(j));
            beta_lp(j,k) = atan2(vm_shipfix_lp(j,k), u_shipfix_lp(j,k)); 
        
            % unfilterd velocity data
            u_shipfix_raw(j,k) = u_hat(j,k) * cos(psi_hat(j)) + vm_hat(j,k) * sin(psi_hat(j));
            vm_shipfix_raw(j,k) = vm_hat(j,k) * cos(psi_hat(j)) - u_hat(j,k) * sin(psi_hat(j));
        end
    end
    %%
    
    %%% wind calculation
    
    for n = 1:no_anemo
        if n==1
            gamma_anemo(:,n) = deg2rad(data.AnemoFrontDir);
            UA_anemo(:,n)     = data.AnemoFrontSpeed; 
        elseif n==2
            gamma_anemo(:,n) = deg2rad(data.AnemoMidDir);
            UA_anemo(:,n)     = data.AnemoMidSpeed;
        end
        for j= 1:height(data)

        % compute ship speed on Anemotemer of ship fix coordinate
        u_shipfix_anemo(j,n) = u_shipfix_kf(j, chosenGNSS);
        v_shipfix_anemo(j,n) = vm_shipfix_kf(j,chosenGNSS) + r_lowpass(j) * L_anemo(n);
        U_anemo(j,n)         = sqrt(u_shipfix_anemo(j,n)^2+v_shipfix_anemo(j,n)^2);
        beta_anemo(j,n)     =  atan2(v_shipfix_anemo(j,n), u_shipfix_anemo(j,n));
        
        % compute true wind velo and dir(velocity to ground)
        [UT_anemo(j,n), zeta_wind_anemo(j,n)] = apparent2true(UA_anemo(j,n), U_anemo(j,n), gamma_anemo(j,n), beta_anemo(j,n), psi_hat(j));
        [UA_anemo_check(j,n), gamma_anemo_check(j,n)] = true2apparent(UT_anemo(j,n), U_anemo(j,n), zeta_wind_anemo(j,n), beta_anemo(j,n), psi_hat(j));
        end
    end
    % average true wind speed and dir.
    UT_ave = mean(UT_anemo,2);
    
    x_dir = mean(cos(zeta_wind_anemo),2);
    y_dir = mean(sin(zeta_wind_anemo),2);
    
    zeta_wind = atan2(y_dir, x_dir);
    for j=1:length(zeta_wind)
        if zeta_wind(j) < 0
            zeta_wind(j) = zeta_wind(j) + 2*pi;
        end
        %%% comvert true wind to apparent on midship
        [UA_midship(j), gamma_a_midship(j)] = true2apparent(UT_ave(j), U_m_kf(j), zeta_wind(j), beta_kf(j), psi_hat(j));

    end
   
    %%
    %%%    direction of travel to ground on xhat-yhat coordinate of fore GPS;
    %    to check accuracy of GPS direction resolve
        
    for j =1:height(data)
    dir_of_travel_GPS_diff(j) = atan2( vm_hat_GPS_kf(j),u_hat_GPS_kf(j));
    dir_of_travel_GPS(j) = data.GPSTrueNorth_Fore(j)-theta_berth;
        if dir_of_travel_GPS(j) >180
            dir_of_travel_GPS(j) = dir_of_travel_GPS(j) -360;
        end
    end
    %%
    %%%%%%%%%%%%%%%%%%%%%%
    %%%% make plots %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%
    filetype ='png';
    %% GPS quality
    % distance between GPS sensor
    figure(1)
    plot(data.time,dist_f2a_m,'r')
    hold on
    plot(data.time,dist_f23_m,'g')
    plot(data.time,dist_a23_m,'b')
    xlabel("time [s]")
    ylabel("Distance [m]")
    title('distance between GPS sensors')
    legend('fore GPS to aft GPS, actual=2.167', 'fore GPS to 3rd GPS, actual=2.592',...
        'aft GPS to 3rd GPS, actual=0.425','location','east')
    hold off
    figurename=['./figures/',files(i).name,'_1.',filetype];
    saveas(gca,figurename,filetype) 
   
    %%%  GPS quality
    figure(2)
    plot(data.time, data.Quality_Fore,'r-')
    hold on
    plot(data.time, data.Quality_Aft,'b-.')
    plot(data.time, data.Quality_3,'k--')
    legend('Fore', 'Aft', '3rd','location','southeast')
    xlabel('time (s)')
    ylabel('quality index')
    grid on
%     yticks([0 4 5])
    ylim([0 6])
    hold off
    figurename=['./figures/',files(i).name,'_2.',filetype];
    saveas(gca,figurename,filetype) 
    
    %%% Number of sat.
    figure(3)
    plot(data.time, data.SatelliteNumber_Fore,'r-')
    hold on
    plot(data.time, data.SatelliteNumber_Aft,'b-.')
    plot(data.time, data.SatelliteNumber_3,'k--')
    legend('Fore', 'Aft', '3rd','location','southeast')
    xlabel('time (s)')
    ylabel('Number of Satellites obtained')
    grid on
%     yticks([0 4 5])
    ylim([6 16])
    hold off
    figurename=['./figures/',files(i).name,'_3.',filetype];
    saveas(gca,figurename,filetype) 
%%
    % midship location difference by GPS anttena
    figure(4)
    plot(data.time, x_hat_diff(:,1),'r.')
    hold on
    plot(data.time, x_hat_diff(:,2), 'b.')
    plot(data.time, x_hat_diff(:,3),'k.')
    
    hleg = legend('$$|\hat{x}_{m_F}-\hat{x}_{m_A}|$$','$$|\hat{x}_{m_A}-\hat{x}_{m_3}|$$',...
        '$$|\hat{x}_{m_F}-\hat{x}_{m_3}|$$','location','best');
    set(hleg, 'Interpreter','latex')
    xlabel('time (s)')
    ylabel('position (m)')
    ylim([0 1])
    htit = title('Difference of midship $$\hat{x}$$ location between GPS output');
    set(htit, 'Interpreter','latex')     
    hold off
    figurename=['./figures/',files(i).name,'_4.',filetype];
    saveas(gca,figurename,filetype) 
    
    % midship location difference by GPS anttena
    figure(5)
    plot(data.time, y_hat_diff(:,1),'r.')
    hold on
    plot(data.time, y_hat_diff(:,2),'b.')
    plot(data.time, y_hat_diff(:,3),'k.')
    hleg = legend('$$|\hat{y}_{m_F}-\hat{y}_{m_A}|$$','$$|\hat{y}_{m_A}-\hat{y}_{m_3}|$$',...
        '$$|\hat{y}_{m_F}-\hat{y}_{m_3}|$$','location','best');
    set(hleg, 'Interpreter','latex')
    xlabel('time (s)')
    ylabel('position (m)')
    ylim([0 1])
    htit = title('Difference of midship $$\hat{y}$$ location between GPS output');
    set(htit, 'Interpreter','latex') 
    hold off
    figurename=['./figures/',files(i).name,'_5.',filetype];
    saveas(gca,figurename,filetype) 
  %% heading, rudder, drift  
    figure(6)
    plot(data.time, data.RudderAng,'m')
    hold on
    plot(data.time, rad2deg(psi(:)),'c')
    plot(data.time, rad2deg(beta_kf(:,1)), 'r--')
    plot(data.time, rad2deg(beta_kf(:,2)), 'b--')
    plot(data.time, rad2deg(beta_kf(:,3)), 'k--')
    xlabel("time")
    ylabel("degree")
    legend('\delta','\psi (Gryo)','\beta (fore GPS)','\beta (Aft GPS)',...
        '\beta (3rd GPS)','location','northwest')
    title('\delta vs \psi vs \beta vs direction of travel')
    hold off
    
    figurename=['./figures/',files(i).name,'_6.',filetype];
    saveas(gca,figurename,filetype)  
   
%% satellite image
 
    % satellite image
    figure(7)
    geoplot(lat_fore,lon_fore,'r-')
    hold on
    geoplot(lat_aft,lon_aft,'b-')
    geoplot(lat_3,lon_3,'w-')
    legend('fore','aft','3rd','color','k','textcolor','w')
    geobasemap satellite
    geolimits([34.8225 34.82292],[135.52458 135.52556])
    title('trajectory of GPS sensors')
    hold off
    figurename=['./figures/',files(i).name,'_7.',filetype];
    saveas(gca,figurename,filetype) 
 %% torajectory
    % trajectory of GPS
    figure (8)
    plot(y_hat_GPS(:,1), x_hat_GPS(:,1),'r--')
    hold on
    plot(y_hat_GPS(:,2), x_hat_GPS(:,2),'b--')
    plot(y_hat_GPS(:,3), x_hat_GPS(:,3),'k--')
    xlabel ('$$\hat{y}$$ position (m)','interpreter','latex')
    ylabel ('$$\hat{x}$$ position (m)','interpreter','latex')
    legend('trajectoryn of fore GPS','trajectory of aft GPS','trajectory of 3rd GPS',...
        'location','southwest')
    axis equal
    title('trajectory of GPS antenna')
    hold off
    figurename=['./figures/',files(i).name,'_8.',filetype];
    saveas(gca,figurename,filetype)
    
    % trajectory of midship
    figure (9)
    plot(y_hat(:,1), x_hat(:,1),'r-')
    hold on
    plot(y_hat(:,2), x_hat(:,2),'b-')
    plot(y_hat(:,3), x_hat(:,3),'k-')
    xlabel ('$$\hat{y}$$ position (m)','interpreter','latex')
    ylabel ('$$\hat{x}$$ position (m)','interpreter','latex')
    legend('midship traj. by fore GPS','midship traj. by AFT GPS','midship traj. by 3rd GPS',...
        'location','best')
    axis equal
    title('trajectory of midship converted from anntena loc.')
    hold off
    figurename=['./figures/',files(i).name,'_9.',filetype];
    saveas(gca,figurename,filetype)
    
    % time history of Antenna trajectory
    figure(10)
    plot(data.time,x_hat_GPS(:,1),'r-')
    hold on
    plot(data.time,y_hat_GPS(:,1),'r-.')
    plot(data.time,x_hat_GPS(:,2),'b-')
    plot(data.time,y_hat_GPS(:,2),'b-.')
    plot(data.time,x_hat_GPS(:,3),'k-')
    plot(data.time,y_hat_GPS(:,3),'k-.')
    hleg = legend('$$\hat{x}_{GPS_F}$$','$$\hat{y}_{GPS_F}$$','$$\hat{x}_{GPS_A}$$','$$\hat{y}_{GPS_A}$$',...
    '$$\hat{x}_{GPS_3}$$','$$\hat{y}_{GPS_3}$$','location','best');
    set(hleg, 'Interpreter','latex')
    xlabel('time (s)')
    ylabel('position (m)')
    title('time history of trajectory of GPS antenna')
    hold off
    figurename=['./figures/',files(i).name,'_10.',filetype];
    saveas(gca,figurename,filetype)
    
    % time history of midship trajectory
    figure(11)
    plot(data.time,x_hat(:,1),'r-')
    hold on
    plot(data.time,y_hat(:,1),'r-.')
    plot(data.time,x_hat(:,2),'b-')
    plot(data.time,y_hat(:,2),'b-.')
    plot(data.time,x_hat(:,3),'k-')
    plot(data.time,y_hat(:,3),'k-.')
    hleg = legend('$$\hat{x}_{m_F}$$','$$\hat{y}_{m_F}$$','$$\hat{x}_{m_A}$$','$$\hat{y}_{m_A}$$',...
    '$$\hat{x}_{m_3}$$','$$\hat{y}_{m_3}$$','location','best');
    set(hleg, 'Interpreter','latex')
    xlabel('time (s)')
    ylabel('position (m)')
    title('time history of mid ship trajectory converted form GPS loc.')
    hold off
    figurename=['./figures/',files(i).name,'_11.',filetype];
    saveas(gca,figurename,filetype)
   %%
    % ship speed
    figure(12)
    plot(data.time,data.GPSSpeedKM_Fore/3.6,'r')
    hold on
    plot(data.time,U_GPS_lowpass(:,1),'r--')
    plot(data.time,data.GPSSpeedKM_Aft/3.6,'b')
    plot(data.time,U_GPS_lowpass(:,2),'b--')
    plot(data.time,data.GPSSpeedKM_3/3.6,'k')
    plot(data.time,U_GPS_lowpass(:,3),'k--')
    xlabel("time (s)")
    ylabel("velocity (m/s)")
    ylim([-1 1])
    legend('Fore GPS output','\it{U_{GPS_F, num. diff., lp}}',...
        'Aft GPS output','\it{U_{GPS_A, num. diff., lp}}',...
        '3rd GPS output','\it{U_{GPS_3, num. diff., lp}}','location','best')
    title('Speed of GPS, GPS output(dopper) vs numerical diff.')
    hold off
    figurename=['./figures/',files(i).name,'_12.',filetype];
    saveas(gca,figurename,filetype) 

    % filtered u
    figure (13)
    plot(data.time, u_shipfix_kf(:,1),'r-')
    hold on 
    plot(data.time, u_shipfix_lp(:,1),'r--')
    plot(data.time, u_shipfix_kf(:,2),'b-')
    plot(data.time, u_shipfix_lp(:,2),'b--')
    plot(data.time, u_shipfix_kf(:,3),'k-')
    plot(data.time, u_shipfix_lp(:,3),'k--')
    xlabel('time (s)')
    ylabel('filetered velocity of midship(m/s)')
    htit=title('ship fix velocity $u$');
    hleg=legend('kf, Fore','lp, Fore','kf, aft','lp, aft',...
        'kf, 3rd','lp, 3rd','location','best','fontsize',14);
    set(hleg, 'Interpreter','latex')
    set(htit, 'Interpreter','latex') 
    ylim([-1 1])
    hold off
    figurename=['./figures/',files(i).name,'_13.',filetype];
    saveas(gca,figurename,filetype)
    
    % filtered vm
    figure (14)
    plot(data.time, vm_shipfix_kf(:,1),'r-')
    hold on 
    plot(data.time, vm_shipfix_lp(:,1),'r--')
    plot(data.time, vm_shipfix_kf(:,2),'b-')
    plot(data.time, vm_shipfix_lp(:,2),'b--')
    plot(data.time, vm_shipfix_kf(:,3),'k-')
    plot(data.time, vm_shipfix_lp(:,3),'k--')
    xlabel('time (s)')
    ylabel('filtered velocity of midship(m/s)')
    htit=title('ship fix velocity $$v_m$$');
    hleg=legend('kf, Fore','lp, Fore','kf, aft','lp, aft',...
        'kf, 3rd','lp, 3rd','location','best','fontsize',14);
    set(hleg, 'Interpreter','latex')
    set(htit, 'Interpreter','latex') 
    ylim([-0.3 0.3])
    hold off
    figurename=['./figures/',files(i).name,'_14.',filetype];
    saveas(gca,figurename,filetype)
    
    % filtered r
    figure(15)
    plot(data.time, data.YawAngVel,'b.')
    hold on
    plot(data.time, rad2deg(r_lowpass),'k--')
    xlabel('time (s)')
    ylabel('Ang velocity (deg/s)')
    legend('Gyro output','lowpass filtered')
    hold off
    figurename=['./figures/',files(i).name,'_15.',filetype];
    saveas(gca,figurename,filetype) 
    
    

   %%
    % heading by GNSS
    figure(16)
    plot(data.time, rad2deg(psi_hat_GNSS(:,1)),'r.')
    hold on
    plot(data.time, rad2deg(psi_hat_GNSS(:,2)),'b.')
    plot(data.time, rad2deg(psi_hat_GNSS(:,3)),'g.')
    plot(data.time, rad2deg(psi_hat_gyro),'k--')
    xlabel('time (s)')
    ylabel('heading angle to berth (hat coordinate)(degree)');
    legend( 'fore GPS to aft GPS', 'fore GPS to 3rd GPS',...
            'aft GPS to 3rd GPS', 'Gryo heading','location','best')
    hold off   
    figurename=['./figures/',files(i).name,'_16.',filetype];
    saveas(gca,figurename,filetype) 

  %% wind 
    % measured wind velocity
    figure(17)
    plot(data.time,data.AnemoMidSpeed,'r')
    hold on
    plot(data.time,data.AnemoFrontSpeed,'b')
    plot(data.time, UA_midship,'k--')
    xlabel("time")
    ylabel("apprarent velocity (m/s)")
    legend('Anemo mid','Anemo fore','converted to midship','location','best')
    title('apparent wind velocity on ship')
    hold off
    figurename=['./figures/',files(i).name,'_17.',filetype];
    saveas(gca,figurename,filetype) 

    %relative wind direction on ship
    figure(18)
    plot(data.time, rad2deg(gamma_anemo(:,1)),'r.')
    hold on
    plot(data.time, rad2deg(gamma_anemo(:,2)),'b.')
    plot(data.time, rad2deg(gamma_a_midship), 'k.')
    xlabel("time")
    ylabel("degree")
    legend('Anemo Fore','Anemo mid','converted to midship','location','southwest')
    title('apparent wind direction on ship')
    hold off
    figurename=['./figures/',files(i).name,'_18.',filetype];
    saveas(gca,figurename,filetype) 
    
    % computed wind velocity
    figure(19)
    plot(data.time, UT_anemo(:,1),'r')
    hold on
    plot(data.time,  UT_anemo(:,2),'b')
    plot(data.time, UT_ave,'k')
    xlabel("time [s]")
    ylabel("true wind velocity [m/s]")
    legend('by fore anemo','by mid anemo','average of Anemometers','location','best')
    title('true wind velocity computed from anemometer')
    hold off
    figurename=['./figures/',files(i).name,'_19.',filetype];
    saveas(gca,figurename,filetype) 

    % computed true wind diretion zeta_w
    figure(20)
    plot(data.time, rad2deg(zeta_wind_anemo(:,1)),'r.')
    hold on
    plot(data.time, rad2deg(zeta_wind_anemo(:,2)),'b.')
    plot(data.time, rad2deg(zeta_wind),'k.')
    xlabel("time")
    ylabel("\zeta_w (degree)")
    legend('by Anemo Fore','by Anemo mid','average of Anemometers','location','southwest')
    title('time history of true wind direction \zeta_w')
    hold off
    figurename=['./figures/',files(i).name,'_20.',filetype];
    saveas(gca,figurename,filetype) 
 %% histogram of control
    figure(21)
    histogram(data.RudderAng,70)
    set(hleg, 'Interpreter','latex')
    xlabel('\delta (degree)')
    title('Histogram of rudder angle')
    hold off
    figurename=['./figures/',files(i).name,'_21.',filetype];
    saveas(gca,figurename,filetype) 
    
    figure(22)
    histogram(data.PropellerRPM,50)
    xlabel('prop. rev (rpm)') 
    title('Histogram of Prop rev.')
    hold off
    figurename=['./figures/',files(i).name,'_22.',filetype];
    saveas(gca,figurename,filetype) 
%% for wind conversion check 
% % time history of beta at midship and anemomters
%     figure (19)
%     plot(data.time, beta_anemo(:,1),'r')
%     hold on
%     plot(data.time, beta_anemo(:,2),'b')
%     plot(data.time, beta_kf,'k')
%     plot(data.time, beta_lp,'k--')
%     xlabel("time (sec)")
%     ylabel("\beta (radian)")
%     legend('Anemo fore', 'Anemo Mid', 'Midship Kalman F', 'Midship Lowpass F',...
%     'location','northwest')
%     figurename=['./figures/',files(i).name,'_19.',filetype];
%     saveas(gca,figurename,filetype) 
%     
%  % convertion check :relative wind direction on ship
%     figure(20)
%     plot(data.time, rad2deg(gamma_anemo(:,1)),'r--')
%     hold on
%     plot(data.time, rad2deg(gamma_anemo(:,2)),'b')
%     plot(data.time, rad2deg(gamma_anemo_check(:,1)),'r^','markersize',1.5)
%     plot(data.time, rad2deg(gamma_anemo_check(:,2)),'b^','markersize',1.5)
%     xlabel("time (sec)")
%     ylabel("\gamma_a (degree)")
%     legend('Anemo Fore measured','Anemo Mid measured','Anemo Fore converted',...
%         'Anemo Mid converted','location','southwest')
%     title('apparent wind direction on ship')
%     hold off
%      figurename=['./figures/',files(i).name,'_20.',filetype];
%     saveas(gca,figurename,filetype) 
    %% distance between GPS and actual
    
    figure(23)
    plot(data.time,abs(dist_f2a_m-(L_GPS(1)-L_GPS(2))),'r')
    hold on
    plot(data.time,abs(dist_f23_m-(L_GPS(1)-L_GPS(3))),'g')
    plot(data.time,abs(dist_a23_m-(L_GPS(2)-L_GPS(3))),'b')
    xlabel("time [s]")
    ylabel("difference from actual, absolute[m]")
    title('distance between GPS sensors, difference')
    legend('fore GPS to aft GPS, actual=2.167', 'fore GPS to 3rd GPS, actual=2.592',...
        'aft GPS to 3rd GPS, actual=0.425','location','northeast')
    hold off
    ylim([0 0.5])
    figurename=['./figures/',files(i).name,'_23.',filetype];
    saveas(gca,figurename,filetype) 
    
    %% Drift of Gyro on comparison with GNSS heading
    driftgyro=zeros(height(data),no_GPS);
    driftgyro(:,1) = rad2deg(psi_hat_gyro(:)-psi_hat_GNSS(:,1));
    driftgyro(:,2) = rad2deg(psi_hat_gyro(:)-psi_hat_GNSS(:,2));
    driftgyro(:,3) = rad2deg(psi_hat_gyro(:)-psi_hat_GNSS(:,3));
    
    figure(24)
   plot(data.time, driftgyro(:,3),'g.')
        hold on
    plot(data.time,driftgyro(:,2),'b.')
    plot(data.time, driftgyro(:,1),'r.')
     xlabel('time (s)')
    ylabel('Drift(degree)');
    legend( 'aft GPS to 3rd GPS', 'fore GPS to 3rd GPS',...
            'fore GPS to aft GPS','location','best')
    title('Drift of Gyro on comparison with GNSS heading')
    ylim([-15, 15])
    hold off   
    figurename=['./figures/',files(i).name,'_24.',filetype];
    saveas(gca,figurename,filetype) 
    %% save output files

    % name variables
    n_prop = data.PropellerRPM/60;  % rpm 2 rps
    delta_rudder = deg2rad(data.RudderAng);
    typeGPS = chosenGNSS * ones(height(data),1);
    typeheading = useheading * ones(height(data),1);   
    
    % criate output table and save
    varNames ={'t [s]','x_position_mid [m]','u_velo [m/s]', 'y_position_mid [m]', 'vm_velo [m/s]', 'psi_hat [rad]',...
    'r_angvelo [rad/s]', 'n_prop [rps]', 'delta_rudder [rad]','wind_velo_relative_mid [m/s]',...
    'wind_dir_relative_mid [rad]','wind_velo_true [m/s]','wind_dir_true [rad]',...
    'u_velo_raw [m/s]','vm_velo_raw [m/s]','psi [rad]','psi_hat_GNSS_f2a [rad]','psi_hat_GNSS_f23 [rad]','psi_hat_GNSS_a23 [rad]','r_angvelo_raw [rad/s]','chosen_GPS','chosen_headind_method'};
    savetable = table(data.time, x_hat(:,chosenGNSS), u_shipfix_kf(:,chosenGNSS), y_hat(:,chosenGNSS), vm_shipfix_kf(:,chosenGNSS), psi_hat,...
    r_lowpass, n_prop, delta_rudder, UA_midship, gamma_a_midship, UT_ave,zeta_wind,...
    u_shipfix_raw(:,chosenGNSS), vm_shipfix_raw(:,chosenGNSS), psi, psi_hat_GNSS(:,1),psi_hat_GNSS(:,2),psi_hat_GNSS(:,3),...
    deg2rad(data.YawAngVel),typeGPS, typeheading,'variablenames',varNames); 

    savename = ['./processed_csv/stopping_processed_',files(i).name, '.csv'];
    writetable(savetable,savename)
 end

%% 線形カルマンフィルタのFunction文
function [ xhat_new, P_new, G ] = kf( A, B, Bu, C, Q, R, u, y, xhat, P )
% KF 線形カルマンフィルタの更新式
% [ xhat_new, P_new, G ] = kf( A, B, B1, C, Q, R, u, y, xhat, P )
% 線形カルマンフィルタの推定値更新を行う
% ----  引数  ----
% A,B,h,C: 対象システム
% x(k+1) = Ax(k)+Bv(k)+Bu u(k)
% y(k) = C7x(k)+w(k)
% のシステム行列
% Q,R:雑音v,wの共分散行列, v,wは正規白色雑音で
% E[v(k)]= E[w(k)]=0
% E[v{k}'v{k}]=Q,E[w(k)'w(k)]=R
% であることを想定
% u: 状態更新前時点での制御入力 u(k-1)
% y: 状態更新後時点での制御入力 y(k)
% xhat, P: 更新前の状態推定量 xhat(k-1)・誤差共分散行列P(k-1)
% ----  戻り値  ----
% xhat_new: 更新後の状態推定値xhat(k)
% P_new: 更新後の誤差共分散行列P(k)
% G: カルマンゲイン G(k)
% 
% 列ベクトルに変換
xhat=xhat(:); u=u(:);y=y(:);
% 事前推定値
xhatm = A*xhat + Bu*u; % 状態
Pm = A*P*A'+B*Q*B'; % 共分散行列
% カルマンゲイン行列
G = Pm * C/(C' * Pm *C +R);
% 事後推定値
xhat_new = xhatm + G * (y-C'*xhatm);% 状態
P_new = (eye(size(A)) - G * C') * Pm; % 誤差共分散
end
%%
function [ dist, co ] = Cal_dist_co( ALong, ALat )
%****************************************************************************    
%  compute disntance traveled and course of ship by of 2 latitude and longitude of points using 漸長緯度航法
%  input are DICIMAL DEGREE (ddd.dddd)
%
%   Input variables
%   ALong: Longitude vector of two points [start end](deg.)
%   ALat:  Latitude vector of two point [start end](deg.)
%
%   Output variables
%   dist: distance traveled（min, to convert to meter, *1852）
%   co: course (degree, north = 0)
%
%****************************************************************************
    
    %
    DR = pi / 180.0d0;
    RD = 180.0d0 / pi;
    %  地球を真球として計算する
    %  漸長緯度を求める
    AM1= 7915.7045D0 * log10( tan( ( 45.0D0 + ALat( 1 ) / 2.0D0 ) * DR ) );
    AM2= 7915.7045D0 * log10( tan( ( 45.0D0 + ALat( 2 ) / 2.0D0 ) * DR ) );   
    %  漸長緯差を求める(単位：分)
    AMD = AM2 - AM1;
    %  緯差を求める(単位：度)
    ALatD=ALat(2)-ALat(1);
    %  経差を求める(単位：分)
    ALongD=(ALong(2)-ALong(1)) * 60.0D0;
    %  中分緯度を求める
    ALatM=(ALat(2)+ALat(1)) / 2.0D0;
    %  東西距を求める(単位：分)
    ADep = ALongD*cos(ALatM*DR);
    %  平均中分緯度航行で進路の計算
    if abs(ADep) <= eps && abs(ALatD) <= eps
%     if ADep <= eps && ALatD <= eps
        co = 0;
    else
        co = atan2(ADep/60.0D0, ALatD);
%           co = abs(atan(ADep/ALatD/60.0D0));
    end
    if co > 2*pi
        co = co-2*pi;
    elseif co < 0
        co = co + 2*pi;
    end
%     if (ALatD >= 0.0D0) && (ALongD >= 0.0D0)
%         co = co;
%     elseif (ALatD >= 0.0D0) && (ALongD <= 0.0D0)
%         co=2.0D0*pi-co;
%     elseif (ALatD <= 0.0D0) && (ALongD <= 0.0D0)
%         co=1.0D0*pi+co;
%     elseif (ALatD <= 0.0D0) && (ALongD >= 0.0D0)
%         co = 1.0D0*pi-co;
%     end
    co = co * RD;
    %
    %   90度270度近辺では中分緯度航法で計算    
    if ( co >= 89.0D0) && (co <= 91.0D0)
        %  航程を求める
        dist = sqrt( ( ALatD * 60.0D0 ) * ( ALatD * 60.0D0 ) + ( ADep ) * ( ADep ) );
    elseif (co >= 269.0D0) && (co <= 271.0D0)
        %  航程を求める
        dist = sqrt( ( ALatD * 60.0D0 ) * ( ALatD * 60.0D0 ) + ( ADep ) * ( ADep ) );
    else
        %  漸長緯度航法で進路の計算
        if abs(ALongD) <= eps && abs(AMD) <= eps
           co =0;
        else
             co =  atan2( ALongD , AMD  );
%              co = abs( atan( ALongD / AMD ) );
        end
        
        if co > 2*pi
            co = co-2*pi;
        elseif co < 0
            co = co + 2*pi;
        end
%         if (AMD >= 0.0D0) && (ALongD >= 0.0D0)
%             co = co;
%         elseif (AMD >= 0.0D0) && (ALongD <= 0.0D0)
%             co = 2.0D0 * pi - co;
%         elseif (AMD <= 0.0D0) && ( ALongD <= 0.0D0)
%             co=1.0D0 * pi + co;
%         elseif (AMD <= 0.0D0) && ( ALongD >= 0.0D0)
%             co = 1.0D0 * pi-co;
%         end
        co = co * RD;
        %  航程の計算
        dist = abs( ( ALatD ) * 60.0D0 / cos( co * DR ) );
    end
    %
end

%%
function[UA, gamma_a] = true2apparent(UT, U_ship, zeta_wind, beta, psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute apparent(relative) wind (or current) speed and direction from true wind measured on ship
% !!! angle are RADIANS !!!
%
% [input]
% UT      : norm of true wind speed (speed to ground).
% U_ship  : norm of speed of arbitary point on ship [m/s], given by GNSS
% zeta_wind:true wind direction (direction to ground). 0 at x direction of
%           earth fix coordinate. (depend on earth fix coor. system. NOT north = zero)
% beta    : angle between ship heading and U_ship (drift angle) [rad] [0, 2pi].
% psi     : heading of ship related to earth fix coordinate [rad] [-pi, pi].
%
% [output]
% UA      : apparent wind speed (norm) on midship. [m/s]
% gamma_a : apparent wind direction on midship.[rad] [0, 2pi]. 0 at ship heading direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute chi(agle between ship heading and true wind direction)
    chi = mod(3*pi - zeta_wind + psi, 2*pi);
    
    if beta < 0
        beta = beta + 2*pi;
    end
    
    % compute chi_beta(angle between ship traveling direction and true wind)
    chi_beta = chi+beta ;
    
    % convert chi_beta tp [0, 2pi]
    if chi_beta < 0
        chi_beta = chi_beta + 2*pi;
    elseif chi_beta > 2*pi
        chi_beta = chi_beta - 2*pi;
    end
    
    % compute apparent wind speed
    UA = sqrt(UT^2 + U_ship^2 -2*UT*U_ship*cos(chi_beta));
    
    % avoid numerical error by acos (acos input are [-1 1])
    RV = (UA^2 + U_ship^2 -UT^2)/(2*UA*U_ship);
    if abs(RV) >1
       RV = sign(RV);
    end
    
    % compute gamma_beta(angle between ship traveling direction and apparent wind)
    % use if sentence because acos range is (0 pi) and gamma_beta (0 2pi)
    if chi_beta <= pi
        gamma_beta = acos(RV);
    else
        gamma_beta = 2*pi-acos(RV);
    end
    
    % compute apparent wind direction gamma_a 
    gamma_a = gamma_beta + beta; 
    
    % avoid NaN at wind or shipspeed are zero
    if UT < 1e-5
        gamma_a = beta;
    elseif U_ship < 1e-5
        gamma_a = mod(2*pi-chi, 2*pi);
    end
    
    % convert to [0, 2pi]
    if gamma_a > 2*pi
        gamma_a = gamma_a-2*pi;
    elseif gamma_a < 0
        gamma_a = gamma_a + 2*pi;
    end   
end

function[UT, zeta_wind] = apparent2true(UA, U_ship, gamma_a, beta, psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute true wind(or current) speed and direction from apparent wind measured on ship
% !!! angle are RADIANS !!!
%
% [input]
% UA      : apparent wind speed (norm) measure on arbitary point of ship. measured by Anemonometer [m/s]
% U_ship  : norm of speed of arbitary point on ship [m/s], given by GNSS
% gamma_a : apparent wind direction measured by Anemometer at arbitary
%           point of ship. [rad] [0, 2pi]. 0 at ship heading direction.
% beta    : angle between ship heading and U_ship (drift angle) [rad] [0, 2pi].
% psi     : heading of ship related to earth fix coordinate [rad] [-pi, pi].
%
% [output]
% UT      : norm of true wind speed (speed to ground).
% zeta_wind:true wind direction (direction to ground). 0 at x direction of
%           earth fix coordinate. (depend on earth fix coor. system. NOT north = zero)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if beta < 0
        beta = beta + 2*pi;
    end
    % compute gammma_beta(angle between ship traveling direction and apparent wind direction gamma_a)
    gamma_beta = gamma_a-beta;
    if gamma_beta < 0
        gamma_beta =gamma_beta + 2*pi;
    end
    % compute true wind speed 
    UT = sqrt(UA^2 + U_ship^2 -2*UA*U_ship*cos(gamma_beta));
       
    % avoid numerical error by acos (acos input are [-1 1])
    RV = (UT^2 + U_ship^2 -UA^2)/(2*UT*U_ship);
    if abs(RV) >1
        RV = sign(RV);
    end
    
    % compute chi_beta(angle between ship traveling direction and true wind)
    % use if sentence because acos range is (0 pi) and chi_beta (0 2pi)
    if gamma_beta <= pi
        chi_beta = acos(RV);
    else
        chi_beta = 2*pi-acos(RV);
    end
    
    % convert to true wind direction
    zeta_wind = psi-chi_beta + beta + pi;
    
    % convert to [0, 2pi]
    if zeta_wind > 2*pi
        zeta_wind = zeta_wind - 2*pi;
    elseif zeta_wind < 0
        zeta_wind = zeta_wind + 2*pi;
    end
    % avoid NaN at wind is zero
    if UT < 1e-5
        zeta_wind = 0;
    end
end